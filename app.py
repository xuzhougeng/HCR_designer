import os
import shutil
import gzip
from datetime import datetime
from Bio import SeqIO
import pandas as pd
from werkzeug.utils import secure_filename
from zipfile import ZipFile
import glob

from flask import Flask, request, render_template, send_file, jsonify

from scripts.splint import create_primer as create_splint_primer
from scripts.hcr import create_primer as create_hcr_primer
from scripts.snail import create_primer as create_snail_primer
from scripts.utils import load_alias, load_fasta, generate_unique_id, create_unique_zip
from scripts.tcr import PrimerDesignConfig, design_prime_set

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


@app.route('/hcr', methods=['GET', 'POST'])
def hcr():
    if request.method == 'POST':
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        # Create output directory if it doesn't exist
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'hcr_results')
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None

        # Handle file upload
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])  # remove '.gz' from filename
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)  # Remove the original .gz file
                    file_path = decompressed_path  # Update file_path to point to the decompressed file
                
                blast_db = file_path

        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        # 获取表单参数
        gc_min = float(request.form.get('min_gc', 40.0))
        gc_max = float(request.form.get('max_gc', 60.0))
        tm_min = float(request.form.get('min_tm', 47.0))
        tm_max = float(request.form.get('max_tm', 53.0))

        # HCR特定参数
        probe_size = int(request.form.get('probe_size', 50))
        polyN = int(request.form.get('polyN', 5))
        initiator_type = request.form.get('initiator_type', 'B1')
        kmer = int(request.form.get('kmer', 8))
        prefix = request.form.get('name', 'probe')  # 使用name作为prefix，如果没有则默认为'probe'

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
        
        # 设计探针 - 直接传递所需参数
        probe_sets = create_hcr_primer(
            seq=seq,
            prefix=prefix,
            probe_size=probe_size,
            polyN=polyN,
            min_gc=gc_min,
            max_gc=gc_max,
            min_tm=tm_min,
            max_tm=tm_max,
            initiator_type=initiator_type,
            kmer=kmer
        )
        
        # 使用固定的文件名
        probe_csv = os.path.join(output_dir, "probes.csv")
        
        # 生成探针结果文件
        #output_primer_sets(probe_sets, probe_csv, blast_db=blast_db)
        
        # 准备要打包的文件列表
        files = [probe_csv]
        
        # 如果存在 BLAST 详细结果文件，也添加到文件列表中
        blast_details = probe_csv.rsplit('.', 1)[0] + '_blast_details.txt'
        if os.path.exists(blast_details):
            files.append(blast_details)

        # 创建 README 文件
        readme_path = os.path.join(output_dir, "readme.txt")
        with open(readme_path, 'w') as f:
            f.write(f"HCR Probe Design Results\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Probe Size: {probe_size} bp\n")
            f.write(f"- Poly N Length: {polyN} bp\n")
            f.write(f"- GC Content: {gc_min}%-{gc_max}%\n")
            f.write(f"- Melting Temperature: {tm_min}°C-{tm_max}°C\n")
            f.write(f"- Initiator Type: {initiator_type}\n")
            f.write(f"- k-mer Size: {kmer}\n")
            if blast_db:
                f.write(f"- Reference Genome: {os.path.basename(blast_db)}\n")
        files.append(readme_path)
        
        # 生成随机的ZIP文件名
        unique_id = generate_unique_id()
        zip_filename = f"hcr_results_{unique_id}.zip"
        zip_path = os.path.join(output_dir, zip_filename)
        
        # 创建 ZIP 文件
        with ZipFile(zip_path, 'w') as zipf:
            for file in files:
                if os.path.isfile(file):
                    zipf.write(file, os.path.basename(file))
        
        # 清理临时文件
        for file in files:
            if os.path.exists(file):
                os.remove(file)
        
        # 返回 ZIP 文件
        return send_file(zip_path, as_attachment=True)
    
    # GET 请求时，获取可用的参考基因组列表并传递给模板
    available_genomes = get_available_genomes()
    return render_template('hcr.html', available_genomes=available_genomes)


@app.route('/splint', methods=['GET', 'POST'])
def splint():
    if request.method == 'POST':
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        # Create output directory if it doesn't exist
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'splint_results')
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None

        # Handle file upload
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])  # remove '.gz' from filename
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)  # Remove the original .gz file
                    file_path = decompressed_path  # Update file_path to point to the decompressed file
                
                blast_db = file_path

        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        # 获取表单参数
        probe_size = int(request.form.get('min_probe_size', 17))  # 使用min_probe_size作为probe_size
        gc_min = float(request.form.get('min_gc', 0.3))  # 保持0-1范围
        gc_max = float(request.form.get('max_gc', 0.7))  # 保持0-1范围
        tm_min = float(request.form.get('min_tm', 45))
        tm_max = float(request.form.get('max_tm', 55))
        polyN = int(request.form.get('polyN', 5))
        kmer = int(request.form.get('kmer', 8))
        fluor = request.form.get('fluor', 'AF488')
        prefix = request.form.get('name', 'probe')  # 使用name作为prefix

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
        
        # 设计探针
        probe_df, blast_df = create_splint_primer(
            seq=seq,
            prefix=prefix,
            probe_size=probe_size,
            polyN=polyN,
            min_gc=gc_min,
            max_gc=gc_max,
            min_tm=tm_min,
            max_tm=tm_max,
            fluor=fluor,
            kmer=kmer,
            background=blast_db
        )
        
        # 使用固定的文件名
        probe_csv = os.path.join(output_dir, "probes.csv")
        blast_csv = os.path.join(output_dir, "blast_results.csv")
        
        # 保存探针结果
        probe_df.to_csv(probe_csv, index=False)
        
        # 准备要打包的文件列表
        files = [probe_csv]
        
        # 如果有blast结果，保存它
        if blast_df is not None:
            blast_df.to_csv(blast_csv, index=False)
            files.append(blast_csv)

        # 创建 README 文件
        readme_path = os.path.join(output_dir, "readme.txt")
        with open(readme_path, 'w') as f:
            f.write(f"Splint-FISH Probe Design Results\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Probe Size: {probe_size} bp\n")
            f.write(f"- Poly N Length: {polyN} bp\n")
            f.write(f"- GC Content: {gc_min*100}%-{gc_max*100}%\n")
            f.write(f"- Melting Temperature: {tm_min}°C-{tm_max}°C\n")
            f.write(f"- Fluorophore: {fluor}\n")
            f.write(f"- k-mer Size: {kmer}\n")
            if blast_db:
                f.write(f"- Reference Genome: {os.path.basename(blast_db)}\n")
        files.append(readme_path)
        
        # 生成随机的ZIP文件名
        unique_id = generate_unique_id()
        zip_filename = f"splint_results_{unique_id}.zip"
        zip_path = os.path.join(output_dir, zip_filename)
        
        # 创建 ZIP 文件
        with ZipFile(zip_path, 'w') as zipf:
            for file in files:
                if os.path.isfile(file):
                    zipf.write(file, os.path.basename(file))
        
        # 清理临时文件
        for file in files:
            if os.path.exists(file):
                os.remove(file)
        
        # 返回 ZIP 文件
        return send_file(zip_path, as_attachment=True)
    
    # GET 请求时，获取可用的参考基因组列表并传递给模板
    available_genomes = get_available_genomes()
    return render_template('splint.html', available_genomes=available_genomes)

@app.route('/snail', methods=['GET', 'POST'])
def snail():
    if request.method == 'GET':
        available_genomes = get_available_genomes()
        return render_template('snail.html', available_genomes=available_genomes)

    if request.method == 'POST':
        # print the form
        print(request.form)
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        probe_size = int(request.form['probe_size'])  # 使用表单中的probe_size参数
        min_gap = int(request.form['min_gap'])
        
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        min_tm = float(request.form['min_tm'])
        max_tm = float(request.form['max_tm'])
        min_comp = int(request.form['min_comp'])

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

        # Create output directory if it doesn't exist
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'snail_results')
        os.makedirs(output_dir, exist_ok=True)

        # background file upload, for blast
        blast_db = None

        # Handle file upload
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])  # remove '.gz' from filename
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)  # Remove the original .gz file
                    file_path = decompressed_path  # Update file_path to point to the decompressed file
                
                blast_db = file_path
        
        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

        # 设计探针
        probe_df = create_snail_primer(
            seq=seq,
            prefix=name,
            probe_size=probe_size,
            polyN=5,  # 使用默认值
            min_gc=min_gc,
            max_gc=max_gc,
            min_tm=min_tm,
            max_tm=max_tm,
            fulor="AF488",  # 使用默认值
            kmer=8  # 使用默认值
        )
        
        # 使用固定的文件名
        probe_csv = os.path.join(output_dir, "probes.csv")
        
        # 保存探针结果
        probe_df.to_csv(probe_csv, index=False)
        
        # 准备要打包的文件列表
        files = [probe_csv]

        # 创建 README 文件
        readme_path = os.path.join(output_dir, "readme.txt")
        with open(readme_path, 'w') as f:
            f.write(f"SNAIL Probe Design Results\n")
            f.write(f"=========================\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Name: {name}\n")
            f.write(f"- Probe Size: {probe_size}\n")
            f.write(f"- Minimum Gap: {min_gap}\n")
            f.write(f"- GC Content Range: {min_gc}% - {max_gc}%\n")
            f.write(f"- Melting Temperature Range: {min_tm}°C - {max_tm}°C\n")
            f.write(f"- Minimum Complementary Length: {min_comp}\n")
            if blast_db:
                f.write(f"- Reference Genome: {os.path.basename(blast_db)}\n")
        
        files.append(readme_path)

        # 创建ZIP文件
        zip_path = create_unique_zip(files, prefix='snail_results')
        
        return send_file(zip_path, as_attachment=True, download_name=f'snail_results_{name}.zip')
    
    return render_template('snail.html')


@app.route('/tcr', methods=['GET', 'POST'])
def tcr():
    if request.method == 'POST':
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        # Create output directory if it doesn't exist
        output_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'tcr_results')
        os.makedirs(output_dir, exist_ok=True)

        # 处理参考基因组/BLAST数据库
        blast_db = None

        # Handle file upload
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join(app.config['UPLOAD_FOLDER'], filename[:-3])  # remove '.gz' from filename
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)  # Remove the original .gz file
                    file_path = decompressed_path  # Update file_path to point to the decompressed file
                
                blast_db = file_path

        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            blast_db = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        # 获取表单参数
        min_length = int(request.form.get('min_length', 15))
        max_length = int(request.form.get('max_length', 20))
        gc_min = float(request.form.get('min_gc', 40.0))
        gc_max = float(request.form.get('max_gc', 60.0))
        tm_min = float(request.form.get('min_tm', 47.0))
        tm_max = float(request.form.get('max_tm', 53.0))
        min_gap = int(request.form.get('min_gap', 2))
        min_complementary_length = int(request.form.get('min_complementary_length', 5))

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
        
        # 创建配置对象
        config = PrimerDesignConfig(
            min_length=min_length,
            max_length=max_length,
            gc_min=gc_min,
            gc_max=gc_max,
            tm_min=tm_min,
            tm_max=tm_max,
            min_gap=min_gap,
            min_complementary_length=min_complementary_length
        )
        
        # 设计探针
        probe_sets = design_prime_set(seq, config)
        
        # 使用固定的文件名
        probe_csv = os.path.join(output_dir, "probes.csv")
        
        # 生成探针结果文件
        # output_primer_sets(probe_sets, probe_csv, blast_db=blast_db)
        # read the probe_csv
        probe_df = pd.read_csv(probe_csv)
        # creat bed file from  probe_df
        bed_file = probe_csv.rsplit('.', 1)[0] + '.bed'
        # feature: probe_df['Type'] ,  probe_df['Set']
        bed_df = pd.DataFrame(
            {
                'Chr': name,
                'Start': probe_df['Start'],
                'End': probe_df['End'],
                'Feature': probe_df['Type'] + '_Set' + probe_df['Set'].astype(str)
            }
        )
        bed_df.to_csv(bed_file, sep='\t', index=False, header=False)

        # 准备要打包的文件列表
        files = [probe_csv, bed_file]
        
        # 如果存在 BLAST 详细结果文件，也添加到文件列表中
        blast_details = probe_csv.rsplit('.', 1)[0] + '_blast_details.txt'
        if os.path.exists(blast_details):
            files.append(blast_details)

        # 创建 README 文件
        readme_path = os.path.join(output_dir, "readme.txt")
        with open(readme_path, 'w') as f:
            f.write(f"TCR Probe Design Results\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Probe Length: {min_length}-{max_length} bp\n")
            f.write(f"- GC Content: {gc_min}%-{gc_max}%\n")
            f.write(f"- Melting Temperature: {tm_min}°C-{tm_max}°C\n")
            f.write(f"- Minimum Gap: {min_gap} bp\n")
            f.write(f"- Minimum Complementary Length: {min_complementary_length} bp\n")
            if blast_db:
                f.write(f"- Reference Genome: {os.path.basename(blast_db)}\n")
        files.append(readme_path)
        
        # 生成随机的ZIP文件名
        unique_id = generate_unique_id()
        zip_filename = f"tcr_results_{unique_id}.zip"
        zip_path = os.path.join(output_dir, zip_filename)
        
        # 创建 ZIP 文件
        with ZipFile(zip_path, 'w') as zipf:
            for file in files:
                if os.path.isfile(file):
                    zipf.write(file, os.path.basename(file))
        
        # 清理临时文件
        for file in files:
            if os.path.exists(file):
                os.remove(file)
        
        # 返回 ZIP 文件
        return send_file(zip_path, as_attachment=True)
    
    # GET 请求时，获取可用的参考基因组列表并传递给模板
    available_genomes = get_available_genomes()
    return render_template('tcr.html', available_genomes=available_genomes)

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="6789")
    #app.run(debug=True)
