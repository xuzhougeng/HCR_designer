from src.triple_probe import TripleProbeDesigner, ProbeSet, TripleProbeConfig,ProbeOutputHandler
import os
import random
from datetime import datetime
import logging
import argparse
# set logging level info
logging.basicConfig(level=logging.INFO)

# 定义碱基互补对应关系
complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def create_config(**kwargs):
    """创建三探针配置对象
    
    Args:
        **kwargs: 配置参数，可选参数见TripleProbeConfig类定义
        
    Returns:
        TripleProbeConfig: 配置对象
    """
    return TripleProbeConfig(**kwargs)

def save_triplet_probes(triplet_probes, output_file, task_name, BP_ID, delimiter=','):
    """
    将三合一探针保存为CSV格式
    
    Parameters:
    -----------
    triplet_probes : list
        三合一探针列表，每个元素包含 (L, M, R) 序列
    output_file : str
        输出文件名
    task_name : str
        任务名称，用作探针ID的前缀
    delimiter : str, optional
        CSV文件分隔符 (default: ',')
    """
    # 从output_file获取基础文件名（不包含扩展名）
    base_name = output_file.rsplit('.', 1)[0]
    triplet_file = f"{base_name}_triplet.csv"
    
    with open(triplet_file, 'w') as f:
        # 写入表头
        headers = ['Probe ID', 'Type', 'Sequence']
        f.write(delimiter.join(headers) + '\n')
        
        # 写入探针序列
        for i, (L, M, R) in enumerate(triplet_probes, 1):
            # 写入L探针
            f.write(f"{task_name}-{i}-L_{BP_ID}{delimiter}{L}\n")
            # 写入M探针
            f.write(f"{task_name}-{i}-M{delimiter}{M}\n")
            # 写入R探针
            f.write(f"{task_name}-{i}-R{delimiter}{R}\n")


def generate_triplet_probe(sequence: str, 
    task_name:str, 
    BP_ID:str, 
    bridge_probe: str, 
    config: TripleProbeConfig):
    """
    为序列创建TripleProbe探针, 用于HCR

    Parameters:
    -----------
    sequence : str
        目标序列
    bridge_probe : str
        桥接探针序列, 长度必须为19个碱基
    task_name : str
        任务名称，用作探针ID的前缀
    BP_ID : str
        探针ID
    config : TripleProbeConfig
        配置对象

    Returns:
    --------
    list: HCR探针列表，每个元素包含：
        L: L + N + brigde_probe[0:16] + (bridge_probe[17:19]的互补序列) + AAGATA
        M: ACATTA + M
        R: R + TAATGTTATCTT
    """
    if not bridge_probe or len(bridge_probe) != 19:
        raise ValueError("桥接探针必须为19个碱基长度")

    # create probe designer
    designer = TripleProbeDesigner(config)
    
    # 设计探针
    probe_sets = designer.design_probes(sequence)
    
    # 输出探针信息, 包括详细结果和BED格式
    output_handler = ProbeOutputHandler(config.output_dir)
    output_handler.save_probe_sets(probe_sets, "triplet_probe", "triplet_probe", config.blast_db)

    # 获取杂交探针(HCR probe)
    # 随机碱基选择
    random_base = random.choice(['A', 'T', 'C', 'G'])
    
    # 获取bridge_probe最后两个碱基的互补序列
    bridge_end_complement = ''.join(complement[base] for base in bridge_probe[17:19])
    
    triplet_probes = []
    for probe_set in probe_sets:
        # 使用probe_sets中的序列信息
        L = probe_set.left_probe[0]  # 获取左引物序列
        M = probe_set.middle_probe[0]  # 获取中间引物序列
        R = probe_set.right_probe[0]  # 获取右引物序列
        
        # 构建三个探针
        L_probe = L + random_base + bridge_probe[0:17] + bridge_end_complement + "AAGATA"
        M_probe = "ACATTA" + M
        R_probe = R + "TAATGTTATCTT"
        
        triplet_probes.append((L_probe, M_probe, R_probe))

    # 创建 README 文件
    readme_path = os.path.join(config.output_dir, "readme.txt")
    with open(readme_path, 'w') as f:
        f.write(f"TCR Probe Design Results\n")
        f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Parameters:\n")
        f.write(f"- Probe Length: {config.min_length}-{config.max_length} bp\n")
        f.write(f"- GC Content: {config.gc_min}%-{config.gc_max}%\n")
        f.write(f"- Melting Temperature: {config.tm_min}°C-{config.tm_max}°C\n")
        f.write(f"- Minimum Gap: {config.min_gap} bp\n")
        f.write(f"- Minimum Complementary Length: {config.min_complementary_length} bp\n")
        if config.blast_db:
            f.write(f"- Reference Genome: {os.path.basename(config.blast_db)}\n")

    # save triplet probes to csv
    output_file = os.path.join(config.output_dir, "triplet_probe.csv")

    save_triplet_probes(triplet_probes, output_file, task_name, BP_ID, delimiter=',')
    return probe_sets

def main(sequence, name, gene_id, ref_genome, **kwargs):
    sequence = sequence.replace("\n", "")
    bridge_probe = 19 * "A"
    task_name = name
    BP_ID = gene_id
    blast_db = ref_genome
    #blast_db = "upload/Athaliana.cdna.fasta"
    config = create_config(**kwargs)
    print(f"Running TCR probe design for gene {name} with bridge probe {gene_id}")
    generate_triplet_probe(sequence, task_name, BP_ID, bridge_probe, config)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='HCR引物设计工具')
    # 必需参数
    parser.add_argument('--name', required=True, help='任务名称，将用作探针ID的前缀')
    parser.add_argument('--sequence', required=True, help='目标序列')
    parser.add_argument('--gene-id', required=True, help='基因ID')
    
    # 可选参数
    parser.add_argument('--min-length', type=int, default=15, help='最小引物长度 (default: 15)')
    parser.add_argument('--max-length', type=int, default=20, help='最大引物长度 (default: 20)')
    parser.add_argument('--min-gc', type=float, default=40.0, help='最小GC含量百分比 (default: 40.0)')
    parser.add_argument('--max-gc', type=float, default=60.0, help='最大GC含量百分比 (default: 60.0)')
    parser.add_argument('--min-tm', type=float, default=47.0, help='最小熔解温度 (default: 47.0)')
    parser.add_argument('--max-tm', type=float, default=53.0, help='最大熔解温度 (default: 53.0)')
    parser.add_argument('--min-gap', type=int, default=2, help='探针之间最小间距 (default: 2)')
    parser.add_argument('--ref-genome', help='参考基因组路径（用于BLAST分析）')
    parser.add_argument('--output-dir', default='output', help='输出文件夹 (default: output)')
    parser.add_argument('--bridge-probe', help='桥接探针序列，用于生成三合一探针')
    args = parser.parse_args()
    main(
        sequence=args.sequence,
        name=args.name,
        gene_id=args.gene_id,
        ref_genome=args.ref_genome,
        **vars(args)
    )