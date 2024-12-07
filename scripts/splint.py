import argparse
import logging
import os
import tempfile

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast.Applications import NcbiblastnCommandline

from .probe_generator import create_probes

# 设置日志记录器
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class PrimerDesignConfig:
    """Splint探针设计的配置类"""
    def __init__(
        self,
        probe_size: int = 17,      # 探针大小
        polyN: int = 5,            # 连续碱基的最大允许数量
        min_gc: float = 0.3,       # 最小GC含量
        max_gc: float = 0.7,       # 最大GC含量
        min_tm: float = 45.0,      # 最小熔点温度
        max_tm: float = 55.0,      # 最大熔点温度
        kmer: int = 8              # k-mer大小
        
    ):
        self.probe_size = probe_size
        self.polyN = polyN
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.kmer = kmer

        # 探针的固定序列
        self.anchor_seq = "TCTCGTTTTCTGAACTTAGC"
        self.oligo_5p = "CGGTATCAAG"
        self.oligo_3p = "CTGTTTAAGA"

    def to_dict(self):
        """将配置转换为字典格式"""
        return {
            'probe_size': self.probe_size,
            'polyN': self.polyN,
            'min_gc': self.min_gc,
            'max_gc': self.max_gc,
            'min_tm': self.min_tm,
            'max_tm': self.max_tm,
            'kmer': self.kmer
        }

# TODO: 可能存在bug
def save_splint_probes_to_csv(probes, output_file, task_name, delimiter=','):
    """
    将Splint探针保存为CSV格式
    
    Parameters:
    -----------
    probes : list
        探针列表，每个元素包含探针序列信息
    output_file : str
        输出文件名
    task_name : str
        任务名称，用作探针ID的前缀
    delimiter : str, optional
        CSV文件分隔符 (default: ',')
    """
    # 从output_file获取基础文件名（不包含扩展名）
    base_name = output_file.rsplit('.', 1)[0]
    probe_file = f"{base_name}_probes.csv"
    
    with open(probe_file, 'w') as f:
        # 写入表头
        headers = ['Probe ID', 'Sequence']
        f.write(delimiter.join(headers) + '\n')
        
        # 写入探针序列
        for i, probe in enumerate(probes, 1):
            probe_id = f"{task_name}-{i}"
            f.write(f"{probe_id}{delimiter}{probe}\n")

# 可能存在bug
def save_probe_positions_bed(probes, positions, output_file, task_name):
    """
    将探针位置信息保存为BED格式
    
    Parameters:
    -----------
    probes : list
        探针列表
    positions : list
        探针位置列表，每个元素为(start, end)元组
    output_file : str
        输出文件名
    task_name : str
        任务名称，用作染色体名
    """
    # 从output_file获取基础文件名（不包含扩展名）
    base_name = output_file.rsplit('.', 1)[0]
    bed_file = f"{base_name}_positions.bed"
    
    with open(bed_file, 'w') as f:
        for i, ((start, end), probe) in enumerate(zip(positions, probes), 1):
            f.write(f"{task_name}\t{start}\t{end}\tProbe-{i}\n")

def analyze_blast_results(primer, db_path):
    """
    分析BLAST结果
    """
    """
    对单个引物进行BLAST分析并统计不同错配数量的匹配数
    只考虑完全长度匹配的情况，gap也计入错配数

    Parameters:
    -----------
    primer : str
        引物序列
    db_path : str
        BLAST数据库路径

    Returns:
    --------
    tuple: (dict, list)
        - dict: 包含不同错配数量的统计结果
        - list: 详细的匹配信息列表，每个元素为(subject_id, mismatches, gaps)
    """
    primer_length = len(primer)
    results = {}
    detailed_matches = []  # 存储详细的匹配信息
    
    # 创建临时文件存储查询序列
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp:
        temp.write(f">query\n{primer}\n")
        query_file = temp.name

    try:
        # 运行BLAST
        output = NcbiblastnCommandline(
            query=query_file,
            db=db_path,
            outfmt="6 qseqid sseqid qstart qend sstart send nident length mismatch gapopen qseq sseq",
            word_size=10,
            task="blastn-short",
            dust="no",
            perc_identity=80
        )()[0]

        # 处理BLAST结果
        for line in output.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 12:
                continue

            # 解析BLAST结果字段
            subject_id = fields[1]     # 目标序列ID
            length = int(fields[7])    # 匹配长度
            mismatches = int(fields[8])# 错配数
            gaps = int(fields[9])      # gap数
            
            # 只考虑完全长度匹配的情况
            if length == primer_length:
                # 计算总错配数（包括gaps）
                total_mismatches = mismatches + gaps
                
                # 更新统计结果
                if total_mismatches not in results:
                    results[total_mismatches] = 0
                results[total_mismatches] += 1
                
                # 保存详细匹配信息
                detailed_matches.append((subject_id, mismatches, gaps))

    finally:
        # 清理临时文件
        os.unlink(query_file)

    # 如果没有任何匹配结果，返回空结果
    if not results:
        return {0: 0}, []

    return results, detailed_matches

def get_detailed_blast_results(primer, db_path):
    """
    获取详细的BLAST比对结果
    
    Parameters:
    -----------
    primer : str
        引物序列
    db_path : str
        BLAST数据库路径
    
    Returns:
    --------
    list: 包含详细比对信息的记录列表
    """
    # 创建临时文件用于BLAST输入
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write(f">query\n{primer}\n")
        query_file = f.name

    # 设置BLAST参数
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=db_path,
        outfmt="6 qseqid sseqid qstart qend sstart send nident length mismatch gapopen qseq sseq",
        word_size=4,
        task="blastn-short",
        dust="no",
        perc_identity=80
    )

    # 运行BLAST
    stdout, stderr = blastn_cline()
    
    # 删除临时文件
    os.unlink(query_file)
    
    # 解析BLAST结果
    detailed_results = []
    
    for line in stdout.split('\n'):
        if not line.strip():
            continue
            
        fields = line.split('\t')
        if len(fields) < 12:
            continue
            
        subject_id = fields[1]
        mismatches = int(fields[8])
        query_seq = fields[10]
        subject_seq = fields[11]
        
        if mismatches <= 4:  # 只记录4个或更少错配的结果
            detailed_results.append({
                'subject_id': subject_id,
                'mismatches': mismatches,
                'alignment': f"Query:  {query_seq}\nSubject: {subject_seq}"
            })
    
    return detailed_results


def design_prime_set(sequence, task_name, config, blast_db):
    """
    Parameters:
    -----------
    sequence : str
        目标序列
    config : PrimerDesignConfig
        引物设计参数配置
        
    Returns:
    --------
    list: 符合条件的探针组合列表，每个元素为包含2个引物信息的元组，
          每个引物信息包含：(序列, 起始位置, 结束位置, Tm值, GC含量)
    """        
    probes = create_probes(sequence, 
                           probe_size= config.probe_size * 2, 
                           inner_gap=0, 
                           polyN=config.polyN, 
                           min_gc=config.min_gc, 
                           max_gc=config.max_gc, 
                           min_tm=config.min_tm, 
                           max_tm=config.max_tm, 
                           k=config.kmer)
    # 将探针转换为所需格式
    primer_sets = []
    for pos, probe in probes.items():
        # 将探针分成两半
        half_size = len(probe) // 2
        left_probe = probe[:half_size]
        right_probe = probe[half_size:]
        
        # 计算每个探针的Tm值和GC含量
        try:
            left_tm = mt.Tm_NN(left_probe, nn_table=mt.DNA_NN4)
            right_tm = mt.Tm_NN(right_probe, nn_table=mt.DNA_NN4)
            
            left_gc = (left_probe.count('G') + left_probe.count('C')) / len(left_probe)
            right_gc = (right_probe.count('G') + right_probe.count('C')) / len(right_probe)
            
            # 添加到结果列表
            primer_sets.append(
                ((left_probe, pos, pos + half_size, left_tm, left_gc),
                 (right_probe, pos + half_size, pos + len(probe), right_tm, right_gc))
            )
        except Exception as e:
            logging.error(f"计算探针参数时出错: {e}")
            continue
            
    logging.info(f'生成了 {len(primer_sets)} 个探针组合')
    return primer_sets
    
    
        


def main(sequence=None, config=None, blast_db=None, output_file="splint_probes.txt", 
         bridge_probe=None, BP_ID=None, task_name=None):
    """
    主函数
    
    Parameters:
    -----------
    sequence : str, optional
        目标序列，如果为None则使用默认序列
    config : SplintConfig, optional
        探针设计参数配置，如果为None则使用默认配置
    blast_db : str, optional
        BLAST数据库路径
    output_file : str, optional
        输出文件名 (default: "splint_probes.txt")
    task_name : str, optional
        任务名称，用于生成探针ID
    """
    if sequence is None:
        sequence = "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    
    if config is None:
        config = PrimerDesignConfig()
    
    if task_name is None:
        task_name = "test"
    
    try:
        # 设计探针
        primer_sets = design_prime_set(sequence, config)
        logger.info(f"设计了 {len(primer_sets)} 个引物组合")
        
        # 保存引物组合
        save_primer_sets_with_blast(primer_sets, output_file, blast_db)
        
        # 保存引物位置信息为BED格式
        save_primer_positions_bed(primer_sets, output_file, task_name)

        if bridge_probe:
            try:
                triplet_probes = create_splint_probes(primer_sets, bridge_probe)
                logger.info(f"生成了 {len(triplet_probes)} 个Splint探针组合")
                save_splint_probes_to_csv(triplet_probes, output_file, task_name, BP_ID)
            except ValueError as e:
                logger.error(f"生成Splint探针时出错: {str(e)}")
                
    except Exception as e:
        logger.error(f"运行过程中出错: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Splint探针设计工具')
    # 必需参数
    parser.add_argument('--name', required=True, help='任务名称，将用作探针ID的前缀')
    parser.add_argument('--seq', required=True, help='目标序列')
    
    # 可选参数
    parser.add_argument('--probe-size', type=int, default=17, help='探针大小 (default: 17)')
    parser.add_argument('--poly-n', type=int, default=5, help='连续碱基的最大允许数量 (default: 5)')
    parser.add_argument('--min-gc', type=float, default=0.3, help='最小GC含量 (default: 0.3)')
    parser.add_argument('--max-gc', type=float, default=0.7, help='最大GC含量 (default: 0.7)')
    parser.add_argument('--min-tm', type=float, default=45.0, help='最小熔点温度 (default: 45.0)')
    parser.add_argument('--max-tm', type=float, default=55.0, help='最大熔点温度 (default: 55.0)')
    parser.add_argument('--fluor', default='AF488', help='荧光标记 (default: AF488)')
    parser.add_argument('--kmer', type=int, default=8, help='k-mer大小 (default: 8)')
    parser.add_argument('--ref-genome', help='参考基因组路径（用于BLAST分析）')
    parser.add_argument('--output', default='splint_probes.txt', help='输出文件名 (default: splint_probes.txt)')
    
    args = parser.parse_args()
    
    # 创建配置对象
    config = SplintConfig(
        probe_size=args.probe_size,
        polyN=args.poly_n,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        fluor=args.fluor,
        kmer=args.kmer
    )
    
    # 运行主函数
    main(sequence=args.seq,
         config=config,
         blast_db=args.ref_genome,
         output_file=args.output,
         task_name=args.name)



# 兼容性代码, 暂时是为了app.py的兼容性
# 17 bp +  2nt gap + 17 bp 
def create_primer(seq, prefix, probe_size=17, polyN=5, min_gc=0.3, max_gc=0.7, min_tm=45, max_tm=55, fluor: str = "AF488", kmer:int = 8, background=None):
    
    """设计splint的探针序列

    输入数据为 cds或者cdna的序列
    输出数据为splint的探针序列

    第一步: 基于csd序列设计探针 ( 分开的单一探针的Tm值, SplintPLP在50度左右 > 45 < 55)
    第二步: 获取探针
    第三步: 对探针反省互补
    """
    from .utils import create_blastn_db,blastn,distance_stat

    filler_seq = "TCTCGTTTTCTGAACTTAGC"

    fluor_probe = {
        "AF488": "TCGCGCTTGGTATAATCGCT",
        "Cy3": "AGTAGCCGTGACTATCGACT",
        "Texas": "TGCGTCTATTTAGTGGAGCC",
        "Cy5": "CCTCAATGCTGCTGCTGTACTAC"    
    }

    fluor_probe_name = {
        "AF488": "P1",
        "Cy3": "P2",
        "Texas": "P3",
        "Cy5": "P4"
    }


    # prober_size = left + right for SPLINT is 34
    probes = create_probes(seq, probe_size= probe_size * 2, inner_gap=0,  polyN=polyN, min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm, k=kmer)

    color_seq = fluor_probe[fluor]
    probe_name_suffix = fluor_probe_name[fluor]

    oligo_5p = "CGGTATCAAG"
    oligo_3p = "CTGTTTAAGA"

    probes_pos = []
    probes_list = []
    P1_name_list = []
    P1_list = []
    P1_tm_list = []
    P2_name_list = []
    P2_list = []
    P2_tm_list = []

    count = 0

    # BLAST result of P1 and P2
    blast_list = []
    # BLAST result stat, dist = 0, 1, 2, 3, 4
    blast_stat = []


    if background is not None:
        dbname = create_blastn_db(background)

    for pos,seq in probes.items():
        count += 1
        probes_pos.append( int(pos) + 1)
        probes_list.append(seq)

        probe = Seq(seq) 
        probe_5p = probe[:probe_size].reverse_complement() 
        probe_3p = probe[-probe_size:].reverse_complement()
        
        primer_5p_tm = mt.Tm_NN(probe_5p)
        primer_3p_tm = mt.Tm_NN(probe_3p)

        primer_5p = f"{probe_5p}{filler_seq}{oligo_5p}"
        primer_3p = f"{oligo_3p}{color_seq}{probe_3p}"

        P1_list.append(primer_5p)
        P2_list.append(primer_3p)

        P1_name = f"{prefix}-{count}-5{probe_name_suffix}"
        P1_name_list.append(P1_name)
        P2_name = f"{prefix}-{count}-3{probe_name_suffix}"
        P2_name_list.append(P2_name)

        if background is not None and dbname:
            probe_name = f"{prefix}-{count}"
            blast_res = blastn(probe_name, str(probe), background)
            blast_res['qlen'] = len(probe)
            blast_list.append(blast_res)
            blast_stat.append(distance_stat(blast_res))

        P1_tm_list.append( primer_5p_tm )
        P2_tm_list.append( primer_3p_tm )

    probe_df = pd.DataFrame({
        "probe_pos" : probes_pos,
        "probe_seq" : probes_list,
        "P1_name": P1_name_list,
        "P1": P1_list,
        "P1_Tm": P1_tm_list,
        "P2_name": P2_name_list,
        "P2": P2_list,
        "P2_Tm": P2_tm_list,
        }
    )
    if len(blast_list) > 0:
        probe_df["blast_stat"] = blast_stat

    blast_df = None
    
    if len(blast_list) > 0:
        blast_df = pd.concat(blast_list)
    
    return probe_df, blast_df
    

