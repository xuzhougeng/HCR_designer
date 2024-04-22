import gzip
import io
import uuid
import os


import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as bn
from Bio.Blast.Applications import NcbimakeblastdbCommandline as makeblastdb
from zipfile import ZipFile

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

def generate_unique_filenames():
    unique_id = uuid.uuid4().hex
    output_filename = f"output_{unique_id}.csv"
    p1_blast_filename = f"p1_blast_{unique_id}.csv"
    p2_blast_filename = f"p2_blast_{unique_id}.csv"
    return output_filename, p1_blast_filename, p2_blast_filename


def create_unique_zip(files, base_dir='results'):
    unique_id = uuid.uuid4().hex
    zip_filename = f"output_{unique_id}.zip"
    zip_path = os.path.join(base_dir, zip_filename)
    
    with ZipFile(zip_path, 'w') as zipf:
        for file in files:
            if os.path.isfile(file):
                zipf.write(file)
                os.remove(file)  # Optionally remove the file after adding it to the zip
            else:
                print(f"Warning: The file {file} does not exist and will not be included in the zip.")
                
    return zip_path


def create_blastn_db(fasta_file):
    """
    :param fasta_file: fasta文件
    """

    # check if the blast db already exists
    db_name = os.path.abspath(fasta_file)
    if os.path.exists(f"{db_name}.nhr"):
        return db_name

    try:
        cline = makeblastdb(
            #cmd = "/opt/biosoft/ncbi-blast-2.10.1+/bin/makeblastdb",
            dbtype="nucl",
            input_file=fasta_file,
            title=db_name,
            out=db_name,
        )
        cline()
    except:
        return False
    
    return db_name

def blastn(seq_id, seq, blastn_db, blastn_evalue=10, blastn_word_size=7):
    """
    :param seq: 欲blast的序列
    :param blastn_db: blastn数据库
    :param blastn_evalue: evalue
    :param blastn_word_size: word_size
    :param blastn_num_threads: 线程数
    """

    # create temp file for the sequence using tempfile
    import tempfile
    seq_file = tempfile.mktemp()
    with open(seq_file, "w") as f:
        f.write(f">{seq_id}\n{seq}\n")

    cline = bn(
            #cmd = "/opt/biosoft/ncbi-blast-2.10.1+/bin/blastn",
            query= seq_file,
            db= blastn_db,
            outfmt=6,
            task="blastn-short",
            evalue=blastn_evalue,
            word_size=blastn_word_size,
            #num_threads=blastn_num_threads,

        )  # this uses biopython's blastn formatting function and creates a commandline compatible command
    (
        stdout,
        _,
    ) = (
        cline()
    )  # cline() calls the string as a command and passes it to the command line, outputting the blast results to one variable and errors to the other

    ## From results of blast creating a numpy array (and Pandas database)
    try:
        blastresult = pd.read_csv(io.StringIO(stdout), delimiter="\t", header=None)
        blastresult.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        # 增加qlen, plen
        blastresult["qlen"] = np.abs(blastresult["qend"] - blastresult["qstart"] ) + 1
        blastresult["plen"] = np.abs(blastresult["send"] - blastresult["sstart"] ) + 1
    except pd.errors.EmptyDataError:
        blastresult = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","qlen","plen" ])
    
    return blastresult