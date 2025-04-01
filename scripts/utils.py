import gzip
import uuid
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import os

def load_fasta(file_path):
    opener = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with opener(file_path, mode) as fasta_file:
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, 'fasta')}
    
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

def create_blast_db(fasta_path):
    """Create BLAST database if it doesn't exist
    
    Parameters:
    -----------
    fasta_path : str
        Path to FASTA file
        
    Returns:
    --------
    str: Path to BLAST database
        
    Raises:
    -------
    ValueError: If sequence validation fails
    """
    # Check if BLAST database files exist
    db_files = [f"{fasta_path}.nhr", f"{fasta_path}.nin", f"{fasta_path}.nsq"]
    if not all(os.path.exists(f) for f in db_files):
        # Validate sequences before creating database
        sequences = load_fasta(fasta_path)
        
        # Check for duplicate IDs
        if len(sequences) != len(set(sequences.keys())):
            raise ValueError("Duplicate sequence IDs found in input file")
            
        # Check sequence format
        for seq_id, seq in sequences.items():
            # Check if sequence is uppercase
            if not seq.isupper():
                raise ValueError(f"Sequence {seq_id} contains lowercase characters")
                
            # Check if sequence contains only valid DNA characters
            if not all(c in 'ATCG' for c in seq):
                raise ValueError(f"Sequence {seq_id} contains invalid DNA characters")
        
        # Create BLAST database
        cline = NcbimakeblastdbCommandline(
            dbtype="nucl",
            input_file=fasta_path,
            out=fasta_path
        )
        cline()
    return fasta_path


