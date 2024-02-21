from flask import Flask, request, render_template, send_file, jsonify
from scripts.splint import create_primer as create_splint_primer
from scripts.hcr import create_primer as create_hcr_primer
from scripts.snail import create_primer as create_snail_primer
from Bio import SeqIO
import gzip
import pandas as pd
app = Flask(__name__)

# This dictionary will store your sequences
# The keys are sequence IDs, and the values are the sequences themselves
# SEUQENCE_ID_FILE_PATH = "sequence_id.txt"

# gene_id_list = []
# with open(SEUQENCE_ID_FILE_PATH, "r") as f:
#     for line in f:
#         gene_id_list.append(line.strip())

def load_fasta(file_path):

    sequences = {}
    with gzip.open(file_path, 'rt') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            sequences[record.id] = str(record.seq)
        
    return sequences

def load_alias(file_path):
    gene_alias = {}
    with gzip.open(file_path, 'r') as f:
        for line in f:
            line = line.decode('latin-1').strip()
            parts = line.strip().split("\t")
            locus_name = parts[0]
            symbol = parts[1]
            symbol = symbol.upper()
            gene_alias[symbol] = locus_name.upper()
    return gene_alias


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
        
        blast = request.form['blast'] == "y"
        blastdb = None
        if blast:
            #blastdb = request.form['blastdb']
            #print(blastdb)
            #if blastdb == 'Arabidopsis':
            blastdb = 'db/Athaliana.fa'

        output = 'output.csv'
        probe_df = create_hcr_primer(seq, name, probe_size * 2, polyN, min_gc, max_gc, min_tm, max_tm, initiator_type)
        probe_df.to_csv(output)
        return send_file(output, as_attachment=True)

    return render_template('hcr.html')

@app.route('/splint', methods=['GET', 'POST'])
def splint():
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
            probe_df_dict[key_name] =  create_splint_primer(seq, name,probe_size,  polyN, min_gc, max_gc, min_tm, max_tm, fluor_type)

        probe_df = pd.concat(probe_df_dict, axis=0)
        probe_df.to_csv(output)
        return send_file(output, as_attachment=True)
    
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
            probe_df_dict[key_name] =  create_snail_primer(seq, name,probe_size,  polyN, min_gc, max_gc, min_tm, max_tm, fluor_type)

        probe_df = pd.concat(probe_df_dict, axis=0)
        probe_df.to_csv(output)
        return send_file(output, as_attachment=True)
    
    return render_template('snail.html')

if __name__ == '__main__':
    app.run(host="0.0.0.0", port="6789")
    #app.run(debug=True)
