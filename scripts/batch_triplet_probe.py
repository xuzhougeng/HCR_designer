from scripts.utils import load_fasta
import os
from datetime import datetime
import multiprocessing as mp
from functools import partial
import zipfile
import argparse
from .create_triplet_probe import main as tcr_main

# column two, BRIDGE ID and sequence
probe_table_file = "./resources/probe_table.txt"
bridge_seq_dict = {}
for line in open(probe_table_file):
    bridge_id, sequence = line.strip().split(" ")
    bridge_seq_dict[bridge_id] = sequence

# Load cDNA sequences from TAIR10
FASTA_FILE_PATH = {
    "TAIR10_cdna" : "./resources/Athaliana.cdna.fasta.gz",
}

# BLAST database path
BLAST_DB = "./uploads/Athaliana_cdna.fasta"
#BLAST_DB = None

cdna_dict = load_fasta(FASTA_FILE_PATH['TAIR10_cdna'])

gene_list = """AT3G25730
AT1G68840
AT1G13260
AT1G25560
AT1G01010
AT1G01260
AT1G01720
AT1G02230
AT1G03040
AT1G03970
AT1G05710
AT1G05805
AT1G06070
AT1G06850
AT1G07530
AT1G07900
AT1G07980
AT1G08000
AT1G08010
AT1G08970
AT1G10120
AT1G10170
AT1G10240
AT1G10585"""

def save_cdna_sequence(gene, sequence, output_dir):
    """
    保存cDNA序列到FASTA文件
    
    Parameters:
    -----------
    gene : str
        基因ID
    sequence : str
        cDNA序列
    output_dir : str
        输出目录
    """
    fasta_file = os.path.join(output_dir, f"{gene}_cdna.fasta")
    with open(fasta_file, 'w') as f:
        f.write(f">{gene}\n")
        # 每行输出60个碱基
        for i in range(0, len(sequence), 60):
            f.write(f"{sequence[i:i+60]}\n")

def process_gene(gene, bridge_seq_dict, gene_to_bp_id, cdna_dict, blast_db):
    """处理单个基因的函数"""
    # 获取bridge probe ID和序列
    bp_id = f"BP_{int(gene_to_bp_id[gene]):04d}"
    bridge_probe = bridge_seq_dict.get(bp_id)
    
    if not bridge_probe:
        print(f"Warning: No bridge probe found for {bp_id}, skipping {gene}")
        return
        
    # 获取基因序列
    sequence = cdna_dict.get(gene)
    if not sequence:
        print(f"Warning: No sequence found for gene {gene}, skipping")
        return
    
    # 创建输出目录
    output_dir = os.path.join("./output", gene)
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存cDNA序列到FASTA文件
    save_cdna_sequence(gene, sequence, output_dir)
    
    # 设置输出文件路径
    output_file = os.path.join(output_dir, f"{gene}_probes.txt")
    
    print(f"Processing gene {gene} with bridge probe {bp_id}")
    
    # 运行TCR探针设计
    print(f"Running TCR probe design for gene {gene} with bridge probe {bp_id}")
    tcr_main(
        sequence=sequence,  
        name=gene,
        gene_id=bp_id,
        ref_genome=blast_db,
        output_dir=output_dir
    )

def zip_output(output_dir="output"):
    """将output目录打包为zip文件"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_filename = f"output_{timestamp}.zip"
    
    with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, output_dir)
                zipf.write(file_path, arcname)
    
    return zip_filename

if __name__ == "__main__":
    # 将基因列表字符串转换为列表
    genes = [gene.strip() for gene in gene_list.split('\n') if gene.strip()]
    
    # 创建基因到BP_ID的映射
    gene_to_bp_id = {gene: i+1001 for i, gene in enumerate(genes)}
    
    # 创建进程池
    num_processes = min(mp.cpu_count(), len(genes))  # 使用CPU核心数或基因数量中的较小值
    pool = mp.Pool(num_processes)
    
    # 准备偏函数，固定其他参数
    process_gene_partial = partial(
        process_gene,
        bridge_seq_dict=bridge_seq_dict,
        gene_to_bp_id=gene_to_bp_id,
        cdna_dict=cdna_dict,
        blast_db=BLAST_DB
    )
    
    try:
        # 使用进程池处理基因
        pool.map(process_gene_partial, genes)
    finally:
        pool.close()
        pool.join()
    
    # 打包输出文件
    zip_file = zip_output()
    print(f"\nCompleted processing all genes")
    print(f"Results have been saved to: {zip_file}")
