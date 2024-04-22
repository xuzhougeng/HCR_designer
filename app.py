
import os
import shutil
import gzip

from Bio import SeqIO
import pandas as pd
from werkzeug.utils import secure_filename


from flask import Flask, request, render_template, send_file, jsonify

from scripts.splint import create_primer as create_splint_primer
from scripts.hcr import create_primer as create_hcr_primer
from scripts.snail import create_primer as create_snail_primer
from scripts.utils import load_alias, load_fasta, generate_unique_id, create_unique_zip

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


@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')


@app.route('/hcr', methods=['GET', 'POST'])
def hcr():
    if request.method == 'POST':
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        probe_size = int(request.form['probe_size']) # single probe size, total probe size will be probe_size * 2
        initiator_type = str(request.form['initiator'])
        polyN = int(request.form['polyN'])
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        min_tm = float(request.form['min_tm'])
        max_tm = float(request.form['max_tm'])

        kmer = int(request.form['kmer'])


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
        
        output = 'output.csv'
        probe_df = create_hcr_primer(seq, name, probe_size * 2, polyN, min_gc, max_gc, min_tm, max_tm, initiator_type, kmer)
        probe_df.to_csv(output)
        return send_file(output, as_attachment=True)

    return render_template('hcr.html')

@app.route('/splint', methods=['GET', 'POST'])
def splint():
    if request.method == 'POST':
        # print the form
        print(request.form)
        name = request.form['name']
        seq = request.form['seq'].upper()
        gene_id = request.form['geneID'].upper()

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        min_probe_size = int(request.form['min_probe_size']) # single probe size, total probe size will be probe_size * 2 defaul 17
        max_probe_size = int(request.form['max_probe_size']) # single probe size, total probe size will be probe_size * 2 defaul 17

        polyN = int(request.form['polyN'])
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        min_tm = float(request.form['min_tm'])
        max_tm = float(request.form['max_tm'])

        kmer = int(request.form['kmer'])

        fluor_type = request.form['fluor']
        
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

        # background file upload, for blast
        background = None

        # Handle file upload
        if 'ref_genome' in request.files:
            file = request.files['ref_genome']
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join( app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)

                # Check if the file is a gzip file and decompress it
                if filename.endswith('.gz'):
                    decompressed_path = os.path.join( app.config['UPLOAD_FOLDER'] , filename[:-3])  # remove '.gz' from filename
                    with gzip.open(file_path, 'rb') as f_in:
                        with open(decompressed_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(file_path)  # Remove the original .gz file
                    file_path = decompressed_path  # Update file_path to point to the decompressed file
                
                background = file_path
        
        # Alternatively, check if an existing genome was selected
        if 'existing_ref_genome' in request.form and request.form['existing_ref_genome'] != '':
            existing_genome = request.form['existing_ref_genome']
            background = os.path.join(app.config['UPLOAD_FOLDER'], existing_genome)  # Path to the selected existing file

        probe_df_dict = {}
        blast_df_dict = {}

        for probe_size in range(min_probe_size, max_probe_size + 1):

            key_name = f'p_{probe_size}'
            res = create_splint_primer(seq, name,probe_size,  polyN, min_gc, max_gc, min_tm, max_tm, fluor_type, kmer, background)
            probe_df_dict[key_name] = res[0]
            blast_df_dict[key_name] = res[1]


        probe_df = pd.concat(probe_df_dict, axis=0)
        blast_df = pd.concat(blast_df_dict, axis=0)

        # Generate unique filenames for this instance of data processing
        unique_id = generate_unique_id()
        output_csv = f"output_{unique_id}.csv"
        blast_csv = f"blast_{unique_id}.csv"

        # Saving dataframes to unique csv files
        probe_df.to_csv(output_csv)
        blast_df.to_csv(blast_csv)
        # List of files to be zipped
        files = [output_csv, blast_csv]
        
        zip_path = create_unique_zip(files, "results")
        return send_file(zip_path, as_attachment=True) 
    
    return render_template('splint.html')

@app.route('/snail', methods=['GET', 'POST'])
def snail():
    if request.method == 'POST':
        # print the form
        print(request.form)
        name = request.form['name']
        seq = request.form['seq']
        gene_id = request.form['geneID'].upper()

        if gene_id in gene_alias:
            gene_id = gene_alias[gene_id]

        min_probe_size = int(request.form['min_probe_size']) # single probe size, total probe size will be probe_size * 2 defaul 17
        max_probe_size = int(request.form['max_probe_size']) # single probe size, total probe size will be probe_size * 2 defaul 17

        polyN = int(request.form['polyN'])
        min_gc = float(request.form['min_gc'])
        max_gc = float(request.form['max_gc'])
        min_tm = float(request.form['min_tm'])
        max_tm = float(request.form['max_tm'])

        kmer = int(request.form['kmer'])

        fluor_type = request.form['fluor']
        
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
        
        # blast = request.form['blast'] == "y"
        # blastdb = None
        # if blast:
        #     #blastdb = request.form['blastdb']
        #     #print(blastdb)
        #     #if blastdb == 'Arabidopsis':
        #     blastdb = 'db/Athaliana.fa';

        output = 'output.csv'

        probe_df_dict = {}
        for probe_size in range(min_probe_size, max_probe_size + 1):

            key_name = f'p_{probe_size}'
            probe_df_dict[key_name] =  create_snail_primer(seq, name,probe_size,  polyN, min_gc, max_gc, min_tm, max_tm, fluor_type, kmer)

        probe_df = pd.concat(probe_df_dict, axis=0)
        probe_df.to_csv(output)
        return send_file(output, as_attachment=True)
    
    return render_template('snail.html')

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="6789")
    #app.run(debug=True)
