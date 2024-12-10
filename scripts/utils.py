import gzip
import uuid
from Bio import SeqIO

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

def save_sequence(gene, sequence, file_path):
    """
    保存序列到FASTA文件
    
    Parameters:
    -----------
    gene : str
        基因ID
    sequence : str
        序列
    file_path : str
        输出文件路径
    """
    with open(file_path, 'w') as f:
        f.write(f">{gene}\n")
        # 每行输出60个碱基
        for i in range(0, len(sequence), 60):
            f.write(f"{sequence[i:i+60]}\n")

def generate_unique_id():
    unique_id = uuid.uuid4().hex
    return unique_id


