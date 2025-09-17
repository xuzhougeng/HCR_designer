import os
import shutil
import gzip
from datetime import datetime
from werkzeug.utils import secure_filename
from zipfile import ZipFile
import csv
import re
import json

from flask import Flask, request, render_template, send_file, jsonify

from scripts.utils import load_alias, load_fasta, generate_unique_id, save_sequence, create_blast_db
from scripts.create_triplet_probe import main as create_triplet_probe
from scripts.create_split_probe import main as create_split_probe
from scripts.bp_suggestion import design_multiple_probes
from scripts.bp_suggestion import query_conflicting_probes
from scripts.bp_suggestion import BridgeProbeSystem
from scripts.bp_suggestion import analyze_input_conflicts
from src.common.sequence_utils import calculate_tm, calculate_gc_content, is_valid_probe, check_poly_n




app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'gz'}

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


# shortest isoform will be used for probe design
FASTA_FILE_PATH = {
    "TAIR10_cdna" : "resources/Athaliana.cdna.fasta.gz",
    "TAIR10_cds" : "resources/Athaliana.cds.fasta.gz"
}

GENE_ALIAS_PATH = {
    "TAIR10" : "resources/gene_aliases_20220930.txt.gz"
}

cds_dict = load_fasta(FASTA_FILE_PATH['TAIR10_cds'])
cdna_dict = load_fasta(FASTA_FILE_PATH['TAIR10_cdna'])
gene_alias = load_alias(GENE_ALIAS_PATH['TAIR10'])

probe_table_file = "resources/probe_table.txt"
bridge_seq_dict = {}
for line in open(probe_table_file):
    bridge_id, sequence = line.strip().split()
    bridge_seq_dict[bridge_id] = sequence

def get_available_genomes():
    """获取可用的参考基因组列表"""
    genome_dir = os.path.join(app.config['UPLOAD_FOLDER'])
    genomes = []
    if os.path.exists(genome_dir):
        for file in os.listdir(genome_dir):
            if file.endswith(('.fa', '.fasta', '.fna')):
                genomes.append({
                    'filename': file,
                    'display_name': file.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
                })
    return genomes


def filter_probe_by_llm(probe_data):
    """
    目前还有问题，需要修改Prompt
    """
    # 构造系统提示词
    import requests

    DEEPSEEK_API_URL = os.getenv("DEEPSEEK_API_URL")
    API_KEY = os.getenv("DEEPSEEK_API_KEY")
    model = os.getenv("DEEPSEEK_MODEL")

    try:
        # 验证输入数据格式
        if not isinstance(probe_data, dict) or 'probe_sets' not in probe_data:
            raise ValueError("输入数据格式错误：需要包含 probe_sets 字段的字典")

        probe_sets = probe_data['probe_sets']
        if not isinstance(probe_sets, list):
            raise ValueError("probe_sets 必须是列表格式")

        # 构建用户提示
        probe_sets_json = json.dumps(probe_sets, indent=2, ensure_ascii=False)
        
        user_prompt = f"""你是一个生物信息学专家，请严格按以下规则，对探针组合进行筛选(一个探针组合包含2到3个probe序列)：
        1. 探针的特异性排序，优先考虑每个探针对中，每个探针错配≤2的组合
        2. 各个探针组合的间距≥15bp
        3. 覆盖度最大化
        最终输出包含每个组合的筛选结果，必须要有两个字典keep和reason（keep字段为True或False，reason字段用中文说明）。如下是JSON格式的输入:\n\n{probe_sets_json}"""

        # API调用
        headers = {
            "Authorization": f"Bearer {API_KEY}",
            "Content-Type": "application/json"
        }
        
        payload = {
            "model": model,
            "messages": [
                {"role": "user", "content": user_prompt}
            ],
            "stream": False
        }
        
        response = requests.post(DEEPSEEK_API_URL, headers=headers, json=payload)
        response.raise_for_status()
        
        # 解析响应
        response_json = response.json()
        if not isinstance(response_json, dict):
            raise ValueError("API返回的响应不是字典格式")
            
        content = response_json.get('choices', [{}])[0].get('message', {}).get('content', '')
        reasoning_content = response_json.get('choices', [{}])[0].get('message', {}).get('reasoning_content', '')
        if not content:
            raise ValueError("API响应中缺少必要的内容")
            
        # 从响应内容中提取JSON部分
        json_start = content.find('[')
        json_end = content.rfind(']') + 1
        if json_start == -1 or json_end == 0:
            raise ValueError("无法从响应内容中提取有效的JSON数据")
            
        json_str = content[json_start:json_end]
        
        # 解析JSON
        result = json.loads(json_str)
        if not isinstance(result, list):
            raise ValueError("解析后的JSON不是列表格式")
        
        # 更新原始数据中的探针结果
        for probe in probe_sets:
            if not isinstance(probe, dict):
                continue
                
            probe_id = str(probe.get("id", ""))
            # 在result中查找匹配的探针
            filtered = next(
                (item for item in result if str(item.get("id", "")) == probe_id),
                None
            )
            
            if filtered:
                probe["keep"] = filtered.get("keep", False)
                probe["reason"] = filtered.get("reason", "未知原因")
            else:
                probe["keep"] = False
                probe["reason"] = "未找到对应的筛选结果"

        return probe_data, reasoning_content

    except requests.RequestException as e:
        print(f"API请求失败: {str(e)}")
        # 处理API请求错误
        for probe in probe_data.get('probe_sets', []):
            if isinstance(probe, dict):
                probe["keep"] = False
                probe["reason"] = f"API请求失败: {str(e)}"
        return probe_data, "API请求失败，无法获取分析过程。"
        
    except json.JSONDecodeError as e:
        print(f"JSON解析错误: {str(e)}")
        # 处理JSON解析错误
        for probe in probe_data.get('probe_sets', []):
            if isinstance(probe, dict):
                probe["keep"] = False
                probe["reason"] = f"JSON解析错误: {str(e)}"
        return probe_data, "JSON解析错误，无法获取分析过程。"
        
    except Exception as e:
        print(f"处理过程出错: {str(e)}")
        # 处理其他错误
        for probe in probe_data.get('probe_sets', []):
            if isinstance(probe, dict):
                probe["keep"] = False
                probe["reason"] = f"处理过程出错: {str(e)}"
        return probe_data, "处理过程出错，无法获取分析过程。"

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/filter', methods=['GET', 'POST'])
def filter():
    return render_template('filter.html')

@app.route("/help", methods=['GET'])
def help():
    return render_template('help.html')

@app.route('/split', methods=['GET', 'POST'])
def split():
    if request.method == 'POST':
        # 生成随机的文件夹名
        unique_id = generate_unique_id()
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'split_results', unique_id)
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None
        error_message = None
        try:
            if 'ref_genome' in request.files:
                file = request.files['ref_genome']
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file.save(file_path)

                    if filename.endswith('.gz'):
                        decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])
                        with gzip.open(file_path, 'rb') as f_in:
                            with open(decompressed_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(file_path)
                        file_path = decompressed_path
                    
                    # Create BLAST database if needed
                    blast_db = create_blast_db(file_path)

            if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
                existing_genome = request.form['existing_ref_genome']
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)
                # Create BLAST database if needed
                blast_db = create_blast_db(file_path)
        except ValueError as e:
            error_message = str(e)
            return render_template('split.html', error_message=error_message)

        # 获取通用参数
        min_length = int(request.form.get('min_length', 15))
        max_length = int(request.form.get('max_length', 20))
        min_gc = float(request.form.get('min_gc', 40.0))
        max_gc = float(request.form.get('max_gc', 60.0))
        min_tm = float(request.form.get('min_tm', 47.0))
        max_tm = float(request.form.get('max_tm', 53.0))
        min_gap = int(request.form.get('min_gap', 2))
        poly_n = int(request.form.get('poly_n', 4))
        kmer_size = int(request.form.get('kmer_size', 8))
        min_kmer_count = int(request.form.get('min_kmer_count', 2))
        min_complementary_length = int(request.form.get('min_complementary_length', 5))
        max_selected = int(request.form.get('max_selected', 5))

        # 单个处理逻辑
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()
        task_name = request.form['name']
        BP_ID = request.form['bp_id']
        
        # 获取probe sequence
        probe_seq = bridge_seq_dict.get(BP_ID, None)
        if probe_seq is None:
            return jsonify({'error': f'BP_ID {BP_ID} not found in probe table'}), 400

        if len(gene_id) > 0:
            sequence_type = request.form['sequenceType']
            if sequence_type == 'cds':
                if gene_id not in cds_dict:
                    return jsonify({'error': f'GeneID {gene_id} not found in CDS data'}), 400
                seq = cds_dict[gene_id]
            elif sequence_type == 'cdna':
                if gene_id  not in cdna_dict:
                    return jsonify({'error': f'GeneID {gene_id} not found in cDNA data'}), 400
                seq = cdna_dict[gene_id]

        if len(seq) == 0:
            return jsonify({'error': 'empty input sequences'}), 400
        
        create_split_probe(
            name=task_name,
            sequence=seq,
            gene_id=gene_id,
            min_length=min_length,
            max_length=max_length,
            gc_min=min_gc,
            gc_max=max_gc,
            tm_min=min_tm,
            tm_max=max_tm,
            min_gap=min_gap,
            min_complementary_length=min_complementary_length,
            poly_n=poly_n,
            kmer_size=kmer_size,
            min_kmer_count=min_kmer_count,
            ref_genome=blast_db,
            bridge_probe_id=BP_ID,
            bridge_probe=probe_seq,
            output_dir=output_dir,
            max_selected=max_selected
        )

        zip_filename = f"split_results_{unique_id}.zip"
        zip_path = os.path.join(app.config['UPLOAD_FOLDER'], 'split_results', zip_filename)
        
        # get files in output_dir
        files = os.listdir(output_dir)
        fasta_name = f"{gene_id}.fasta"
        fasta_file = os.path.join(output_dir, fasta_name)
        save_sequence(gene_id, seq, fasta_file)
        files.append(fasta_name)
        
        # 创建 ZIP 文件
        with ZipFile(zip_path, 'w') as zipf:
            for file in files:
                if os.path.isfile(os.path.join(output_dir, file)):
                    zipf.write(os.path.join(output_dir, file), os.path.basename(file))
        
        # 清理临时文件
        for file in files:
            if os.path.exists(os.path.join(output_dir, file)):
                os.remove(os.path.join(output_dir, file))
        
        # 返回 ZIP 文件
        return send_file(zip_path, as_attachment=True)
    
    # GET请求处理保持不变
    available_genomes = get_available_genomes()
    return render_template('split.html', available_genomes=available_genomes)

@app.route('/split/batch', methods=['GET', 'POST'])
def split_batch():
    if request.method == 'POST':
        # 生成随机的文件夹名
        unique_id = generate_unique_id()
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'split_results', unique_id)
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None
        ref_genome_path = None
        error_message = None
        try:
            if 'ref_genome' in request.files:
                file = request.files['ref_genome']
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file.save(file_path)

                    if filename.endswith('.gz'):
                        decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])
                        with gzip.open(file_path, 'rb') as f_in:
                            with open(decompressed_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(file_path)
                        file_path = decompressed_path
                    
                    # Create BLAST database if needed
                    blast_db = create_blast_db(file_path)
                    ref_genome_path = file_path

            if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
                existing_genome = request.form['existing_ref_genome']
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)
                # Create BLAST database if needed
                blast_db = create_blast_db(file_path)
                ref_genome_path = blast_db
        except ValueError as e:
            error_message = str(e)
            return render_template('split_batch.html', error_message=error_message)

        # 获取通用参数
        min_length = int(request.form.get('min_length', 15))
        max_length = int(request.form.get('max_length', 20))
        min_gc = float(request.form.get('min_gc', 40.0))
        max_gc = float(request.form.get('max_gc', 60.0))
        min_tm = float(request.form.get('min_tm', 47.0))
        max_tm = float(request.form.get('max_tm', 53.0))
        min_gap = int(request.form.get('min_gap', 2))
        poly_n = int(request.form.get('poly_n', 4))
        kmer_size = int(request.form.get('kmer_size', 8))
        min_kmer_count = int(request.form.get('min_kmer_count', 2))
        min_complementary_length = int(request.form.get('min_complementary_length', 5))
        default_max_selected = int(request.form.get('max_selected', 5))  # 获取默认的max_selected值

        # 批量处理
        if 'batch_file' in request.files:
            batch_file = request.files['batch_file']
            if batch_file:
                # 创建一个子目录存放所有批量结果
                batch_dir = os.path.join(output_dir, 'batch_results')
                os.makedirs(batch_dir, exist_ok=True)
                
                # 读取CSV文件
                import csv
                from io import StringIO
                csv_content = batch_file.read().decode('utf-8')
                csv_reader = csv.reader(StringIO(csv_content))
                
                # 如果有参考基因组，先从中提取所有序列
                ref_sequences = {}
                if ref_genome_path:
                    print(ref_genome_path)
                    ref_sequences = load_fasta(ref_genome_path)
                    # upper all keys
                    ref_sequences = {k.upper(): v for k, v in ref_sequences.items()}
                                
                for row in csv_reader:
                    if len(row) >= 3:  # 至少需要3列
                        task_name, bp_id, gene_id = row[:3]
                        gene_id = gene_id.upper()
                        # 获取max_selected值,如果有第4列就使用第4列的值,否则使用默认值
                        max_selected = default_max_selected
                        if len(row) >= 4:
                            try:
                                max_selected = int(row[3])
                            except (ValueError, TypeError):
                                max_selected = default_max_selected
                        
                        # 获取probe sequence
                        probe_seq = bridge_seq_dict.get(bp_id, None)
                        if probe_seq is None:
                            continue
                            
                        # 获取序列，按优先级尝试不同来源
                        sequence = None
                        
                        # 1. 首先检查基因别名
                        if gene_id in gene_alias:
                            gene_id = gene_alias[gene_id]
                        
                        # 2. 从CDS字典获取
                        if gene_id in cds_dict:
                            sequence = cds_dict[gene_id]
                        # 3. 从cDNA字典获取
                        elif gene_id in cdna_dict:
                            sequence = cdna_dict[gene_id]
                        # 4. 从参考基因组获取
                        elif gene_id in ref_sequences:
                            sequence = ref_sequences[gene_id]
                        # 5. 尝试从参考基因组中模糊匹配
                        elif ref_sequences:
                            # 尝试查找包含gene_id的序列ID
                            for ref_id, ref_seq in ref_sequences.items():
                                if gene_id in ref_id:
                                    sequence = ref_seq
                                    gene_id = ref_id  # 更新为完整的ID
                                    break
                            
                        if sequence:
                            # 为每个任务创建单独的输出目录
                            task_dir = os.path.join(batch_dir, task_name)
                            os.makedirs(task_dir, exist_ok=True)
                            
                            # 保存序列到FASTA文件
                            fasta_file = os.path.join(task_dir, f"{gene_id}.fasta")
                            save_sequence(gene_id, sequence, fasta_file)
                            
                            # 运行split probe设计
                            create_split_probe(
                                name=task_name,
                                sequence=sequence,
                                gene_id=gene_id,
                                min_length=min_length,
                                max_length=max_length,
                                gc_min=min_gc,
                                gc_max=max_gc,
                                tm_min=min_tm,
                                tm_max=max_tm,
                                min_gap=min_gap,
                                min_complementary_length=min_complementary_length,
                                poly_n=poly_n,
                                kmer_size=kmer_size,
                                min_kmer_count=min_kmer_count,
                                ref_genome=blast_db,
                                bridge_probe_id=bp_id,
                                bridge_probe=probe_seq,
                                output_dir=task_dir,
                                max_selected=max_selected  # 使用每行指定的max_selected值
                            )
                        else:
                            print(f"无法找到基因 {gene_id} 的序列")
                
                # 创建ZIP文件
                zip_filename = f"split_batch_results_{unique_id}.zip"
                zip_path = os.path.join(app.config['UPLOAD_FOLDER'], 'split_results', zip_filename)
                
                with ZipFile(zip_path, 'w') as zipf:
                    for root, dirs, files in os.walk(batch_dir):
                        for file in files:
                            file_path = os.path.join(root, file)
                            arcname = os.path.relpath(file_path, batch_dir)
                            zipf.write(file_path, arcname)
                
                # 清理临时文件
                shutil.rmtree(output_dir)
                
                # 返回ZIP文件
                return send_file(zip_path, as_attachment=True)

    # GET请求返回批量处理页面
    available_genomes = get_available_genomes()
    return render_template('split_batch.html', available_genomes=available_genomes)

@app.route('/triplet', methods=['GET', 'POST'])
def triplet():
    if request.method == 'POST':
        # 生成随机的文件夹名
        unique_id = generate_unique_id()
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'triplet_results', unique_id)
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None
        error_message = None
        try:
            if 'ref_genome' in request.files:
                file = request.files['ref_genome']
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file.save(file_path)

                    if filename.endswith('.gz'):
                        decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])
                        with gzip.open(file_path, 'rb') as f_in:
                            with open(decompressed_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(file_path)
                        file_path = decompressed_path
                    
                    # Create BLAST database if needed
                    blast_db = create_blast_db(file_path)

            if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
                existing_genome = request.form['existing_ref_genome']
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)
                # Create BLAST database if needed
                blast_db = create_blast_db(file_path)
        except ValueError as e:
            error_message = str(e)
            return render_template('triplet.html', error_message=error_message)

        # 获取通用参数
        min_length = int(request.form.get('min_length', 15))
        max_length = int(request.form.get('max_length', 20))
        min_gc = float(request.form.get('min_gc', 40.0))
        max_gc = float(request.form.get('max_gc', 60.0))
        min_tm = float(request.form.get('min_tm', 47.0))
        max_tm = float(request.form.get('max_tm', 53.0))
        min_gap = int(request.form.get('min_gap', 2))
        poly_n = int(request.form.get('poly_n', 4))
        kmer_size = int(request.form.get('kmer_size', 8))
        min_kmer_count = int(request.form.get('min_kmer_count', 2))
        min_complementary_length = int(request.form.get('min_complementary_length', 5))
        max_selected = int(request.form.get('max_selected', 5))

        # 原有的单个处理逻辑保持不变
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()
        task_name = request.form['name']
        BP_ID = request.form['bp_id']
        
        # 获取probe sequence
        probe_seq = bridge_seq_dict.get(BP_ID, None)
        if probe_seq is None:
            return jsonify({'error': f'BP_ID {BP_ID} not found in probe table'}), 400

        if len(gene_id) > 0:
            sequence_type = request.form['sequenceType']
            if sequence_type == 'cds':
                if gene_id not in cds_dict:
                    return jsonify({'error': f'GeneID {gene_id} not found in CDS data'}), 400
                seq = cds_dict[gene_id]
            elif sequence_type == 'cdna':
                if gene_id  not in cdna_dict:
                    return jsonify({'error': f'GeneID {gene_id} not found in cDNA data'}), 400
                seq = cdna_dict[gene_id]

        if len(seq) == 0:
            return jsonify({'error': 'empty input sequences'}), 400
        
        create_triplet_probe(
            name=task_name,
            sequence=seq,
            gene_id=gene_id,
            min_length=min_length,
            max_length=max_length,
            gc_min=min_gc,
            gc_max=max_gc,
            tm_min=min_tm,
            tm_max=max_tm,
            min_gap=min_gap,
            min_complementary_length=min_complementary_length,
            poly_n=poly_n,
            kmer_size=kmer_size,
            min_kmer_count=min_kmer_count,
            ref_genome=blast_db,
            bridge_probe_id=BP_ID,
            bridge_probe=probe_seq,
            output_dir=output_dir,
            max_selected=max_selected
        )
        

        zip_filename = f"triplet_results_{unique_id}.zip"
        zip_path = os.path.join(app.config['UPLOAD_FOLDER'], 'triplet_results', zip_filename)
        
        # get files in output_dir
        files = os.listdir(output_dir)
        fasta_name = f"{gene_id}.fasta"
        fasta_file = os.path.join(output_dir, fasta_name)
        save_sequence(gene_id, seq, fasta_file)
        files.append(fasta_name)
        # 创建 ZIP 文件
        with ZipFile(zip_path, 'w') as zipf:
            for file in files:
                if os.path.isfile(os.path.join(output_dir, file)):
                    zipf.write(os.path.join(output_dir, file), os.path.basename(file))
        
        # 清理临时文件
        for file in files:
            if os.path.exists(os.path.join(output_dir, file)):
                os.remove(os.path.join(output_dir, file))
        
        # 返回 ZIP 文件
        return send_file(zip_path, as_attachment=True)
    
    # GET请求处理保持不变
    available_genomes = get_available_genomes()
    return render_template('triplet.html', available_genomes=available_genomes)

@app.route('/triplet/batch', methods=['GET', 'POST'])
def triplet_batch():
    if request.method == 'POST':
        # 生成随机的文件夹名
        unique_id = generate_unique_id()
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'triplet_results', unique_id)
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None
        ref_genome_path = None
        error_message = None
        try:
            if 'ref_genome' in request.files:
                file = request.files['ref_genome']
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    file.save(file_path)

                    if filename.endswith('.gz'):
                        decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])
                        with gzip.open(file_path, 'rb') as f_in:
                            with open(decompressed_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(file_path)
                        file_path = decompressed_path
                    
                    # Create BLAST database if needed
                    blast_db = create_blast_db(file_path)
                    ref_genome_path = file_path

            if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
                existing_genome = request.form['existing_ref_genome']
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)
                # Create BLAST database if needed
                blast_db = create_blast_db(file_path)
                ref_genome_path = blast_db
        except ValueError as e:
            error_message = str(e)
            return render_template('triplet_batch.html', error_message=error_message)

        # 获取通用参数
        min_length = int(request.form.get('min_length', 15))
        max_length = int(request.form.get('max_length', 20))
        min_gc = float(request.form.get('min_gc', 40.0))
        max_gc = float(request.form.get('max_gc', 60.0))
        min_tm = float(request.form.get('min_tm', 47.0))
        max_tm = float(request.form.get('max_tm', 53.0))
        min_gap = int(request.form.get('min_gap', 2))
        poly_n = int(request.form.get('poly_n', 4))
        kmer_size = int(request.form.get('kmer_size', 8))
        min_kmer_count = int(request.form.get('min_kmer_count', 2))
        min_complementary_length = int(request.form.get('min_complementary_length', 5))
        default_max_selected = int(request.form.get('max_selected', 5))  # 获取默认的max_selected值

        # 批量处理
        if 'batch_file' in request.files:
            batch_file = request.files['batch_file']
            if batch_file:
                # 创建一个子目录存放所有批量结果
                batch_dir = os.path.join(output_dir, 'batch_results')
                os.makedirs(batch_dir, exist_ok=True)
                
                # 读取CSV文件
                import csv
                from io import StringIO
                csv_content = batch_file.read().decode('utf-8')
                csv_reader = csv.reader(StringIO(csv_content))
                
                # 如果有参考基因组，先从中提取所有序列
                ref_sequences = {}
                if ref_genome_path:
                    print(ref_genome_path)
                    ref_sequences = load_fasta(ref_genome_path)
                    # upper all keys
                    ref_sequences = {k.upper(): v for k, v in ref_sequences.items()}
                
                for row in csv_reader:
                    if len(row) >= 3:  # 至少需要3列
                        task_name, bp_id, gene_id = row[:3]
                        gene_id = gene_id.upper()
                        
                        # 获取max_selected值,如果有第4列就使用第4列的值,否则使用默认值
                        max_selected = default_max_selected
                        if len(row) >= 4:
                            try:
                                max_selected = int(row[3])
                            except (ValueError, TypeError):
                                max_selected = default_max_selected
                        
                        # 获取probe sequence
                        probe_seq = bridge_seq_dict.get(bp_id, None)
                        if probe_seq is None:
                            continue
                            
                        # 获取序列，按优先级尝试不同来源
                        sequence = None
                        
                        # 1. 首先检查基因别名
                        if gene_id in gene_alias:
                            gene_id = gene_alias[gene_id]
                        
                        # 2. 从CDS字典获取
                        if gene_id in cds_dict:
                            sequence = cds_dict[gene_id]
                        # 3. 从cDNA字典获取
                        elif gene_id in cdna_dict:
                            sequence = cdna_dict[gene_id]
                        # 4. 从参考基因组获取
                        elif gene_id in ref_sequences:
                            sequence = ref_sequences[gene_id]
                        # 5. 尝试从参考基因组中模糊匹配
                        elif ref_sequences:
                            # 尝试查找包含gene_id的序列ID
                            for ref_id, ref_seq in ref_sequences.items():
                                if gene_id in ref_id:
                                    sequence = ref_seq
                                    gene_id = ref_id  # 更新为完整的ID
                                    break
                            
                        if sequence:
                            # 为每个任务创建单独的输出目录
                            task_dir = os.path.join(batch_dir, task_name)
                            os.makedirs(task_dir, exist_ok=True)
                            
                            # 保存序列到FASTA文件
                            fasta_file = os.path.join(task_dir, f"{gene_id}.fasta")
                            save_sequence(gene_id, sequence, fasta_file)
                            
                            # 运行triplet probe设计
                            create_triplet_probe(
                                name=task_name,
                                sequence=sequence,
                                gene_id=gene_id,
                                min_length=min_length,
                                max_length=max_length,
                                gc_min=min_gc,
                                gc_max=max_gc,
                                tm_min=min_tm,
                                tm_max=max_tm,
                                min_gap=min_gap,
                                min_complementary_length=min_complementary_length,
                                poly_n=poly_n,
                                kmer_size=kmer_size,
                                min_kmer_count=min_kmer_count,
                                ref_genome=blast_db,
                                bridge_probe_id=bp_id,
                                bridge_probe=probe_seq,
                                output_dir=task_dir,
                                max_selected=max_selected  # 使用每行指定的max_selected值
                            )
                        else:
                            print(f"无法找到基因 {gene_id} 的序列")
                
                # 创建ZIP文件
                zip_filename = f"triplet_batch_results_{unique_id}.zip"
                zip_path = os.path.join(app.config['UPLOAD_FOLDER'], 'triplet_results', zip_filename)
                
                with ZipFile(zip_path, 'w') as zipf:
                    for root, dirs, files in os.walk(batch_dir):
                        for file in files:
                            file_path = os.path.join(root, file)
                            arcname = os.path.relpath(file_path, batch_dir)
                            zipf.write(file_path, arcname)
                
                # 清理临时文件
                shutil.rmtree(output_dir)
                
                # 返回ZIP文件
                return send_file(zip_path, as_attachment=True)

    # GET请求返回批量处理页面
    available_genomes = get_available_genomes()
    return render_template('triplet_batch.html', available_genomes=available_genomes)

@app.route('/bridge', methods=['GET', 'POST'])
def bridge():
    if request.method == 'POST':
        # 获取表单参数
        input_bp_ids = request.form.get('bp_ids', '').strip()
        n = int(request.form.get('number', 1))
        k = int(request.form.get('kmer', 9))
        bp_start = request.form.get('bp_start', '').strip()
        bp_end = request.form.get('bp_end', '').strip()
        
        # 处理BP_ID列表
        bp_id_list = []
        
        # 处理文件上传
        if 'bp_file' in request.files:
            file = request.files['bp_file']
            if file and file.filename:
                content = file.read().decode('utf-8')
                file_bp_ids = [line.strip() for line in content.splitlines() if line.strip()]
                bp_id_list.extend(file_bp_ids)
        
        # 处理文本框输入
        if input_bp_ids:
            text_bp_ids = [bp.strip() for bp in input_bp_ids.split(',') if bp.strip()]
            bp_id_list.extend(text_bp_ids)
        
        # 去除重复的BP_ID
        bp_id_list = list(dict.fromkeys(bp_id_list))
        
        try:
            result_probes = design_multiple_probes(
                probe_table_file="resources/probe_table.txt",
                input_bp_ids=bp_id_list,
                n=n,
                k=k,
                bp_range=(bp_start, bp_end) if bp_start and bp_end else None
            )
            
            return render_template('bridge.html', 
                                result_probes=result_probes,
                                bp_ids=input_bp_ids,
                                number=n,
                                kmer=k,
                                bp_start=bp_start,
                                bp_end=bp_end)
            
        except Exception as e:
            return render_template('bridge.html', 
                                error=str(e),
                                bp_ids=input_bp_ids,
                                number=n,
                                kmer=k,
                                bp_start=bp_start,
                                bp_end=bp_end)
    
    return render_template('bridge.html', number=1, kmer=9)

@app.route('/bridge/query', methods=['GET', 'POST'])
def bridge_query():
    if request.method == 'POST':
        query_type = request.form.get('query_type', 'bp_id')
        k = int(request.form.get('kmer', 9))
        
        try:
            if query_type == 'bp_id':
                bp_id = request.form.get('bp_id', '').strip()
                if not bp_id:
                    raise ValueError("请输入BP_ID")
                
                conflicts = query_conflicting_probes(
                    probe_table_file="resources/probe_table.txt",
                    bp_id=bp_id,
                    k=k
                )
                return render_template('bridge_query.html',
                                       query_type='bp_id',
                                       conflicts=conflicts,
                                       query_bp_id=bp_id,
                                       query_sequence=bridge_seq_dict.get(bp_id),
                                       kmer=k)
            else:
                sequences = request.form.get('sequences', '').strip().split('\n')
                sequences = [seq.strip() for seq in sequences if seq.strip()]
                
                if not sequences:
                    raise ValueError("请输入至少一个序列")
                
                # 直接使用BridgeProbeSystem的方法
                system = BridgeProbeSystem("resources/probe_table.txt", k=k)
                results = system.analyze_sequences_conflicts(sequences)
                
                return render_template('bridge_query.html',
                                       query_type='sequence',
                                       sequence_results=results)
              
        except Exception as e:
            return render_template('bridge_query.html', 
                                  query_type=query_type,
                                  bp_id=request.form.get('bp_id', ''),
                                  kmer=k,
                                  error=str(e))
    
    # 处理GET请求，检查是否有bp_id参数
    bp_id = request.args.get('bp_id', '').strip()
    if bp_id:
        try:
            k = 9  # 默认k-mer长度
            conflicts = query_conflicting_probes(
                probe_table_file="resources/probe_table.txt",
                bp_id=bp_id,
                k=k
            )
            return render_template('bridge_query.html',
                                   query_type='bp_id',
                                   conflicts=conflicts,
                                   query_bp_id=bp_id,
                                   query_sequence=bridge_seq_dict.get(bp_id),
                                   kmer=k)
        except Exception as e:
            return render_template('bridge_query.html',
                                   query_type='bp_id',
                                   bp_id=bp_id,
                                   kmer=9,
                                   error=str(e))
    
    return render_template('bridge_query.html', query_type='bp_id', kmer=9)

@app.route('/bridge/list', methods=['GET'])
def bridge_list():
    # 获取所有探针信息并按BP_ID排序
    probes = []
    for bp_id, seq in bridge_seq_dict.items():
        tm = round(calculate_tm(seq), 1)  # 计算Tm值并四舍五入到1位小数
        gc = round(calculate_gc_content(seq), 1)  # 计算GC含量并四舍五入到1位小数
        probes.append({
            'bp_id': bp_id,
            'sequence': seq,
            'length': len(seq),
            'tm': tm,
            'gc': gc
        })
    
    # 按BP_ID排序
    probes.sort(key=lambda x: x['bp_id'])
    
    return render_template('bridge_list.html', probes=probes)

@app.route('/bridge/sequence', methods=['GET', 'POST'])
def bridge_sequence():
    if request.method == 'POST':
        k = int(request.form.get('kmer', 9))
        sequences = request.form.get('sequences', '').strip().split('\n')
        sequences = [seq.strip() for seq in sequences if seq.strip()]
        
        try:
            if not sequences:
                raise ValueError("请输入至少一个序列")
            
            # 使用BridgeProbeSystem进行序列分析
            system = BridgeProbeSystem("resources/probe_table.txt", k=k)
            results = system.analyze_sequences_conflicts(sequences)
            
            return render_template('bridge_sequence.html', 
                                  sequence_results=results,
                                  sequences='\n'.join(sequences),
                                  kmer=k)
               
        except Exception as e:
            return render_template('bridge_sequence.html', 
                                  error=str(e),
                                  sequences=request.form.get('sequences', ''),
                                  kmer=k)
      
    return render_template('bridge_sequence.html', kmer=9)

@app.route('/bridge/analyze', methods=['GET', 'POST'])
def bridge_analyze():
    if request.method == 'POST':
        k = int(request.form.get('kmer', 9))
        
        # 处理文件上传
        if 'input_file' in request.files:
            file = request.files['input_file']
            if file:
                content = file.read().decode('utf-8')
            else:
                content = request.form.get('input_content', '').strip()
        else:
            content = request.form.get('input_content', '').strip()
            
        if not content:
            return render_template('bridge_analyze.html',
                                error="Please provide input content",
                                kmer=k)
        
        try:
            results = analyze_input_conflicts(
                probe_table_file="resources/probe_table.txt",
                input_content=content,
                k=k
            )
            
            return render_template('bridge_analyze.html',
                                results=results,
                                input_content=content,
                                kmer=k)
            
        except Exception as e:
            return render_template('bridge_analyze.html',
                                error=str(e),
                                input_content=content,
                                kmer=k)
    
    return render_template('bridge_analyze.html', kmer=9)

@app.route('/bridge/validate', methods=['GET', 'POST'])
def bridge_validate():
    # 默认参数值
    default_params = {
        'min_length': 17,
        'max_length': 20,
        'min_gc': 40.0,
        'max_gc': 60.0,
        'min_tm': 47.0,
        'max_tm': 53.0,
        'poly_n': 4
    }
    
    if request.method == 'POST':
        # 获取用户设置的参数
        params = {
            'min_length': int(request.form.get('min_length', default_params['min_length'])),
            'max_length': int(request.form.get('max_length', default_params['max_length'])),
            'min_gc': float(request.form.get('min_gc', default_params['min_gc'])),
            'max_gc': float(request.form.get('max_gc', default_params['max_gc'])),
            'min_tm': float(request.form.get('min_tm', default_params['min_tm'])),
            'max_tm': float(request.form.get('max_tm', default_params['max_tm'])),
            'poly_n': int(request.form.get('poly_n', default_params['poly_n']))
        }
        
        if 'probe_file' not in request.files:
            return render_template('bridge_validate.html', error="请上传文件", **params)
        
        file = request.files['probe_file']
        if file.filename == '':
            return render_template('bridge_validate.html', error="未选择文件", **params)
            
        if not file.filename.endswith('.csv'):
            return render_template('bridge_validate.html', error="请上传CSV文件", **params)
            
        try:
            # 读取CSV文件
            content = file.read().decode('utf-8')
            reader = csv.reader(content.splitlines())
            next(reader)  # 跳过表头
            
            results = []
            for row in reader:
                if len(row) < 2:
                    continue
                    
                probe_id, sequence = row[0].strip(), row[1].strip()
                is_valid = True
                reasons = []
                
                # 检查序列格式
                if not re.match('^[ATCG]+$', sequence):
                    is_valid = False
                    reasons.append("序列包含非法字符")
                
                # 使用用户设置的参数进行验证
                length = len(sequence)
                if length < params['min_length'] or length > params['max_length']:
                    is_valid = False
                    reasons.append(f"长度不符合要求: {length}bp")
                
                gc = calculate_gc_content(sequence)
                if gc < params['min_gc'] or gc > params['max_gc']:
                    is_valid = False
                    reasons.append(f"GC含量不适合: {gc:.1f}%")
                
                try:
                    tm = calculate_tm(sequence)
                    if tm < params['min_tm'] or tm > params['max_tm']:
                        is_valid = False
                        reasons.append(f"Tm值不适合: {tm:.1f}°C")
                except:
                    is_valid = False
                    reasons.append("Tm值计算失败")
                
                if not check_poly_n(sequence, params['poly_n']):
                    is_valid = False
                    reasons.append(f"存在{params['poly_n']}个以上连续碱基")
                
                results.append({
                    'probe_id': probe_id,
                    'sequence': sequence,
                    'is_valid': is_valid,
                    'reasons': '；'.join(reasons) if reasons else '符合要求'
                })
            
            return render_template('bridge_validate.html', results=results, **params)
            
        except Exception as e:
            return render_template('bridge_validate.html', error=str(e), **params)
            
    return render_template('bridge_validate.html', **default_params)

@app.route('/filter', methods=['GET', 'POST'])
def filter_probes():
    if request.method == 'POST':
        # 检查是否有文件上传
        if 'probe_file' not in request.files:
            return render_template('filter.html', error="请上传文件")
            
        file = request.files['probe_file']
        if file.filename == '':
            return render_template('filter.html', error="未选择文件")
            
        if not file.filename.endswith('.json'):
            return render_template('filter.html', error="请上传JSON文件")
            
        try:
            # 读取JSON文件
            content = file.read().decode('utf-8')
            probe_data = json.loads(content)
            
            # 调用过滤函数
            filtered_data, reasoning_content = filter_probe_by_llm(probe_data)
            
            # 返回结果
            return render_template('filter.html', results=filtered_data, reasoning_content=reasoning_content)
            
        except Exception as e:
            return render_template('filter.html', error=str(e))
    
    return render_template('filter.html')

def calculate_tpm(cds_file, counts_file):
    """计算基因表达量的TPM值
    
    Args:
        cds_file: FASTA格式的CDS序列文件
        counts_file: 两列的文本文件，第一列是基因ID，第二列是count值
        
    Returns:
        dict: 基因ID到TPM值的映射字典
    """
    try:
        # 1. 读取CDS序列，获取每个基因的长度
        cds_content = cds_file.read().decode('utf-8')
        cds_dict = {}
        current_id = None
        current_seq = []
        
        for line in cds_content.splitlines():
            if line.startswith('>'):
                if current_id and current_seq:
                    cds_dict[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif line.strip():
                current_seq.append(line.strip())
        
        if current_id and current_seq:
            cds_dict[current_id] = ''.join(current_seq)
            
        gene_lengths = {gene_id: len(sequence) for gene_id, sequence in cds_dict.items()}
        
        # 2. 读取counts文件
        counts_dict = {}
        counts_content = counts_file.read().decode('utf-8')
        for line in counts_content.splitlines():
            if line.strip():  # 跳过空行
                parts = line.strip().split()
                if len(parts) >= 2:
                    gene_id = parts[0]
                    try:
                        count = float(parts[1])
                        counts_dict[gene_id] = count
                    except ValueError:
                        continue  # 跳过无法转换为数字的行
        
        # 3. 计算每个基因的RPK (Reads Per Kilobase)
        rpk_dict = {}
        for gene_id, count in counts_dict.items():
            if gene_id in gene_lengths:
                # RPK = reads / (length_kb)
                rpk = count / (gene_lengths[gene_id] / 1000.0)
                rpk_dict[gene_id] = rpk
        
        # 4. 计算RPK总和
        sum_rpk = sum(rpk_dict.values())
        
        # 5. 计算每个基因的TPM
        tpm_dict = {}
        if sum_rpk > 0:  # 避免除以0
            scaling_factor = 1_000_000.0 / sum_rpk
            for gene_id, rpk in rpk_dict.items():
                tpm = rpk * scaling_factor
                tpm_dict[gene_id] = tpm
        
        return tpm_dict
        
    except Exception as e:
        print(f"TPM计算错误: {str(e)}")
        return {}

def get_recommended_probe_number(tpm):
    """根据TPM值推荐探针数量
    
    TPM范围对应的探针数量：
    0-1: 10个探针
    1-10: 8个探针
    10-30: 5个探针
    30-50: 5个探针
    50-150: 4个探针
    >150: 3个探针
    """
    if tpm < 1:
        return 10
    elif tpm < 10:
        return 8
    elif tpm < 30:
        return 5
    elif tpm < 50:
        return 5
    elif tpm < 150:
        return 4
    else:
        return 3

def format_fasta_sequences(input_text, clean_headers=True, uppercase=True, filter_invalid=True):
    """
    Format FASTA sequences according to specified options

    Args:
        input_text: FASTA sequence text
        clean_headers: Remove text after first space in headers
        uppercase: Convert sequences to uppercase
        filter_invalid: Remove sequences with non-ATCG bases

    Returns:
        tuple: (formatted_text, stats_dict)
    """
    lines = input_text.strip().split('\n')
    formatted_lines = []
    stats = {
        'original_count': 0,
        'formatted_count': 0,
        'removed_count': 0,
        'headers_cleaned': 0,
        'removed_sequences': []
    }

    current_header = None
    current_sequence = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('>'):
            # Process previous sequence if exists
            if current_header is not None:
                stats['original_count'] += 1
                sequence = ''.join(current_sequence)

                # Apply formatting options
                if uppercase:
                    sequence = sequence.upper()

                # Check for invalid bases if filtering is enabled
                if filter_invalid:
                    valid_bases = set('ATCG')
                    invalid_bases = set(sequence) - valid_bases
                    if invalid_bases:
                        stats['removed_count'] += 1
                        stats['removed_sequences'].append((current_header.strip('>'), f"Contains invalid bases: {', '.join(invalid_bases)}"))
                        current_header = None
                        current_sequence = []
                        continue

                # Add formatted sequence
                formatted_lines.append(current_header)
                formatted_lines.append(sequence)
                stats['formatted_count'] += 1

            # Process new header
            current_header = line
            if clean_headers and ' ' in line:
                # Keep only the part before the first space
                current_header = line.split(' ')[0]
                stats['headers_cleaned'] += 1

            current_sequence = []
        else:
            current_sequence.append(line)

    # Process the last sequence
    if current_header is not None:
        stats['original_count'] += 1
        sequence = ''.join(current_sequence)

        if uppercase:
            sequence = sequence.upper()

        if filter_invalid:
            valid_bases = set('ATCG')
            invalid_bases = set(sequence) - valid_bases
            if invalid_bases:
                stats['removed_count'] += 1
                stats['removed_sequences'].append((current_header.strip('>'), f"Contains invalid bases: {', '.join(invalid_bases)}"))
            else:
                formatted_lines.append(current_header)
                formatted_lines.append(sequence)
                stats['formatted_count'] += 1
        else:
            formatted_lines.append(current_header)
            formatted_lines.append(sequence)
            stats['formatted_count'] += 1

    return '\n'.join(formatted_lines), stats

@app.route('/formatting', methods=['GET', 'POST'])
def formatting():
    if request.method == 'POST':
        try:
            # Get processing options
            clean_headers = 'clean_headers' in request.form
            uppercase = 'uppercase' in request.form
            filter_invalid = 'filter_invalid' in request.form

            # Get input from file upload only
            if 'fasta_file' not in request.files:
                return jsonify({'error': "请上传文件"}), 400

            file = request.files['fasta_file']
            if file.filename == '':
                return jsonify({'error': "未选择文件"}), 400

            if not allowed_file(file.filename):
                return jsonify({'error': "请上传有效的FASTA文件 (.fasta, .fa, .gz)"}), 400

            # Read file content directly from memory, don't save to disk
            if file.filename.endswith('.gz'):
                content = gzip.open(file, 'rt').read()
            else:
                content = file.read().decode('utf-8')

            # Format sequences
            formatted_sequences, stats = format_fasta_sequences(
                content, clean_headers, uppercase, filter_invalid
            )

            if not formatted_sequences.strip():
                return jsonify({'error': "处理后没有有效序列"}), 400

            # Create formatted filename
            original_name = file.filename.rsplit('.', 1)[0]  # Remove extension
            if original_name.endswith('.fasta'):
                original_name = original_name[:-6]  # Remove .fasta
            formatted_filename = f"{original_name}_formatted.fasta"

            # Return formatted file for download using send_file with BytesIO
            from io import BytesIO
            output = BytesIO()
            output.write(formatted_sequences.encode('utf-8'))
            output.seek(0)

            return send_file(
                output,
                as_attachment=True,
                download_name=formatted_filename,
                mimetype='text/plain'
            )

        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('formatting.html')

@app.route('/probe_recommend', methods=['GET', 'POST'])
def probe_recommend():
    if request.method == 'POST':
        # 检查是否有所需的文件
        if 'cds_file' not in request.files or 'counts_file' not in request.files:
            return jsonify({'error': '请上传所有必需的文件'}), 400

        cds_file = request.files['cds_file']
        counts_file = request.files['counts_file']

        # 检查文件是否被选择
        if cds_file.filename == '' or counts_file.filename == '':
            return jsonify({'error': '请选择所有必需的文件'}), 400

        try:
            # 保存文件内容到临时变量
            cds_content = cds_file.read()
            counts_content = counts_file.read()
            
            # 重置文件指针
            cds_file.seek(0)
            counts_file.seek(0)
            
            # 计算所有基因的TPM值
            tpm_dict = calculate_tpm(cds_file, counts_file)
            
            # 从TPM字典中获取所有基因ID
            gene_list = list(tpm_dict.keys())

            # 根据TPM值推荐探针数量
            result = []
            for gene_id in gene_list:
                tpm = tpm_dict[gene_id]
                probe_number = get_recommended_probe_number(tpm)
                result.append({
                    'gene_id': gene_id,
                    'tpm': round(tpm, 2),  # 保留两位小数
                    'probe_number': probe_number
                })

            # 生成CSV文件
            unique_id = generate_unique_id()
            csv_filename = f"probe_recommend_{unique_id}.csv"
            csv_path = os.path.join(app.config['UPLOAD_FOLDER'], csv_filename)
            
            with open(csv_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['gene_id', 'tpm', 'probe_number'])
                writer.writeheader()
                writer.writerows(result)

            # 返回CSV文件
            return send_file(csv_path, as_attachment=True)

        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('probe_recommend.html')

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="6789", debug=True)
