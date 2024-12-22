import os
import shutil
import gzip
from datetime import datetime
from werkzeug.utils import secure_filename
from zipfile import ZipFile

from flask import Flask, request, render_template, send_file, jsonify

from scripts.utils import load_alias, load_fasta, generate_unique_id, save_sequence
from scripts.create_triplet_probe import main as create_triplet_probe
from scripts.create_split_probe import main as create_split_probe
from scripts.bp_suggestion import design_multiple_probes
from scripts.bp_suggestion import query_conflicting_probes
from scripts.bp_suggestion import BridgeProbeSystem
from scripts.bp_suggestion import analyze_input_conflicts
from src.common.sequence_utils import calculate_tm, calculate_gc_content

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
    bridge_id, sequence = line.strip().split(" ")
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

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/split', methods=['GET', 'POST'])
def split():
    if request.method == 'POST':
        # 生成随机的文件夹名
        unique_id = generate_unique_id()
        # Create output directory if it doesn't exist
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'split_results', unique_id)
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None

        # Handle file upload for reference genome
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)
                    file_path = decompressed_path
                
                blast_db = file_path

        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

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

        # 检查是否为批量处理
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
                
                for row in csv_reader:
                    if len(row) >= 3:
                        task_name, bp_id, gene_id = row[:3]
                        print(task_name, bp_id, gene_id)
                        gene_id = gene_id.upper()
                        
                        # 获取probe sequence
                        probe_seq = bridge_seq_dict.get(bp_id, None)
                        if probe_seq is None:
                            print(f"Bridge probe ID {bp_id} not found in database")
                            return jsonify({'error': f'Bridge probe ID {bp_id} not found in database'}), 400
                            
                        # 获取序列
                        if gene_id in gene_alias:
                            gene_id = gene_alias[gene_id]
                            
                        sequence = None
                        if gene_id in cds_dict:
                            sequence = cds_dict[gene_id]
                        elif gene_id in cdna_dict:
                            sequence = cdna_dict[gene_id]
                            
                        if sequence:
                            # 为每个任务创建单独的输出目录
                            task_dir = os.path.join(batch_dir, task_name)
                            os.makedirs(task_dir, exist_ok=True)
                            
                            # Add this line to save sequence to FASTA file
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
                                output_dir=task_dir
                            )
                
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
                
        else:
            # 原有的单个处理逻辑
            seq = request.form['seq']
            gene_id = request.form['geneID'].upper()
            task_name = request.form['name']
            BP_ID = request.form['bp_id']
            print(BP_ID + " " + gene_id)
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
                output_dir=output_dir
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
                
                blast_db = file_path

        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

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
                
                for row in csv_reader:
                    if len(row) >= 3:
                        task_name, bp_id, gene_id = row[:3]
                        gene_id = gene_id.upper()
                        
                        # 获取probe sequence
                        probe_seq = bridge_seq_dict.get(bp_id, None)
                        if probe_seq is None:
                            continue
                            
                        # 获取序列
                        if gene_id in gene_alias:
                            gene_id = gene_alias[gene_id]
                            
                        sequence = None
                        if gene_id in cds_dict:
                            sequence = cds_dict[gene_id]
                        elif gene_id in cdna_dict:
                            sequence = cdna_dict[gene_id]
                            
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
                                output_dir=task_dir
                            )
                
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

        # Handle file upload for reference genome
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
                
                blast_db = file_path

        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

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
            output_dir=output_dir
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
                
                blast_db = file_path

        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

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
                
                for row in csv_reader:
                    if len(row) >= 3:
                        task_name, bp_id, gene_id = row[:3]
                        gene_id = gene_id.upper()
                        
                        # 获取probe sequence
                        probe_seq = bridge_seq_dict.get(bp_id, None)
                        if probe_seq is None:
                            continue
                            
                        # 获取序列
                        if gene_id in gene_alias:
                            gene_id = gene_alias[gene_id]
                            
                        sequence = None
                        if gene_id in cds_dict:
                            sequence = cds_dict[gene_id]
                        elif gene_id in cdna_dict:
                            sequence = cdna_dict[gene_id]
                            
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
                                output_dir=task_dir
                            )
                
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
        tm = round(calculate_tm(seq), 1)  # 计算Tm值并四舍五��到1位小数
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

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="6789", debug=True)
