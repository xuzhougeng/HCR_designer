from Bio import SeqIO
import gzip

import sys

fasta_file = sys.argv[1]

# 假设fasta_file是一个打开的FASTA文件
with gzip.open(fasta_file, "rt") as fh:
    records = list(SeqIO.parse(fh, "fasta"))

# 创建一个字典来存储每个基因的最短isoform
shortest_isoforms = {}

for record in records:
    gene_id = record.id.split(".")[0]
    if gene_id not in shortest_isoforms or len(record.seq) < len(shortest_isoforms[gene_id].seq):
        shortest_isoforms[gene_id] = record

# 输出最短isoforms
for gene_id, record in shortest_isoforms.items():
    print(f">{gene_id}\n{record.seq}")
