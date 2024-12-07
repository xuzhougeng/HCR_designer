"""
TCR引物设计工具

该工具用于设计TCR探针,主要功能包括:

1. 引物设计流程:
- 从左到右扫描目标序列,寻找符合条件的引物组合
- 每个组合包含左引物(L)、中间引物(M)和右引物(R)
- 引物需满足GC含量、Tm值、长度等条件
- 避免引物间的二聚体和发卡结构

2. 主要参数:
- 引物长度: 15-20bp
- GC含量: 40-60%
- Tm值: 47-53°C
- 引物间最小间距: 2bp
- 连续相同碱基数: ≤4

3. 质量控制:
- BLAST比对检查特异性
- 检查引物互补性
- 检查发卡结构
- 检查连续碱基

4. 输出结果:
- 引物序列及位置信息
- BLAST分析结果
- 引物参数(GC含量、Tm值等)
- BED格式的位置文件
- 可选生成三合一探针

5. 使用方法:
python tcr.py --name <任务名> --seq <目标序列> --gene-id <基因ID> [可选参数]


"""

import logging
from Bio.SeqUtils import GC as gc
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO
import subprocess
import tempfile
import argparse
import random
import os

# 配置日志
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

## Helper Functions
def calculate_gc_content(sequence):
    """
    计算序列的GC含量
    
    Parameters:
    -----------
    sequence : str
        DNA序列
    
    Returns:
    --------
    float: GC含量百分比
    """
    return gc(Seq(sequence))

def check_poly_n(seq, n=4):
    """检查序列中是否存在连续的N个相同碱基"""
    for base in ['A', 'T', 'G', 'C']:
        if base * n in seq:
            return False
    return True

def is_complementary(seq1, seq2):
    """检查两个序列是否互补"""
    if len(seq1) != len(seq2):
        return False
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for b1, b2 in zip(seq1, seq2):
        if b2 != pairs.get(b1):
            return False
    return True

def has_hairpin(seq, min_stem_length=4):
    """检查序列是否可能形成发卡结构"""
    seq_length = len(seq)
    for i in range(seq_length - min_stem_length):
        for j in range(i + min_stem_length, seq_length):
            # 获取潜在的茎部序列
            stem1 = seq[i:i+min_stem_length]
            stem2 = seq[j-min_stem_length+1:j+1]
            # 检查是否互补
            if is_complementary(stem1, stem2[::-1]):
                return True
    return False

def check_complementarity(primer1, primer2, min_complementary_length=4):
    """检查两个引物序列之间的互补性"""
    len1, len2 = len(primer1), len(primer2)
    
    # 检查局部互补性
    for i in range(len1 - min_complementary_length + 1):
        for j in range(len2 - min_complementary_length + 1):
            if is_complementary(primer1[i:i+min_complementary_length], 
                              primer2[j:j+min_complementary_length][::-1]):
                return {'has_complementarity': True, 'end_complementarity': False}
    
    # 检查3'端互补性
    end_length = min(5, min_complementary_length)
    end_complementarity = is_complementary(primer1[-end_length:], primer2[-end_length:][::-1])
    
    return {
        'has_complementarity': False,
        'end_complementarity': end_complementarity
    }

def has_dimer_issues(primers, min_complementary_length=4):
    """
    检查一组引物是否存在二聚体问题
    
    Parameters:
    -----------
    primers : list
        引物序列列表
    min_complementary_length : int, optional
        判定为互补的最小连续碱基长度 (default: 4)
        
    Returns:
    --------
    bool: 如果存在二聚体问题则返回True，否则返回False
    
    Example:
        >>> has_dimer_issues(["ATGC", "GCAT"])
        True
    """
    if not primers:
        return False
        
    # 检查所有可能的引物对
    for i in range(len(primers)):
        # 检查自身互补
        if len(primers[i]) >= min_complementary_length * 2:
            result = check_complementarity(primers[i][:len(primers[i])//2], 
                                        primers[i][len(primers[i])//2:], 
                                        min_complementary_length)
            if result['has_complementarity']:
                return True
            
        # 检查与其他引物的互补
        for j in range(i + 1, len(primers)):
            result = check_complementarity(primers[i], primers[j], min_complementary_length)
            if result['has_complementarity']:
                return True
                
    return False

def is_valid_primer(primer_seq, min_length=17, max_length=20, gc_min=40.0, gc_max=60.0, tm_min=47.0, tm_max=53.0, poly_n=4):
    """
    检查引物是否符合要求
    
    Parameters:
    -----------
    primer_seq : str
        引物序列
    min_length : int, optional
        最小引物长度 (default: 17)
    max_length : int, optional
        最大引物长度 (default: 20)
    gc_min : float, optional
        最小GC含量百分比 (default: 40.0)
    gc_max : float, optional
        最大GC含量百分比 (default: 60.0)
    tm_min : float, optional
        最小熔解温度 (default: 47.0)
    tm_max : float, optional
        最大熔解温度 (default: 53.0)
    poly_n : int, optional
        连续相同碱基的最大允许数量 (default: 4)
        
    Returns:
    --------
    bool: 如果引物符合要求则返回True，否则返回False
    
    Example:
        >>> is_valid_primer("ATGC")
        True
        >>> is_valid_primer("ATGC", min_length=10)
        False
    """
    length = len(primer_seq)
    if length < min_length or length > max_length:
        return False

    gc_content = calculate_gc_content(primer_seq)
    if gc_content < gc_min or gc_content > gc_max:
        return False

    try:
        tm = mt.Tm_NN(primer_seq, nn_table=mt.DNA_NN4)
        if tm < tm_min or tm > tm_max:
            return False
    except Exception as e:
        logger.error(f"计算熔点温度时出错: {e}")
        return False

    # 避免连续N个相同的碱基
    seq_str = str(primer_seq)
    bases = ['A', 'T', 'G', 'C']
    for base in bases:
        # 检查是否有大于或等于poly_n个连续碱基
        base_count = 0
        for current_base in seq_str:
            if current_base == base:
                base_count += 1
                if base_count >= poly_n:  # 大于或等于poly_n个连续碱基
                    return False
            else:
                base_count = 0

    return True

def analyze_blast_results(primer, db_path):
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

class PrimerDesignConfig:
    """引物设计参数配置类"""
    def __init__(self,
                 min_length=15,
                 max_length=20,
                 gc_min=40.0,
                 gc_max=60.0,
                 tm_min=47.0,
                 tm_max=53.0,
                 min_gap=15,  # 探针之间最小间距
                 r_m_gap=2,  # R探针和M探针之间的固定间隔
                 min_complementary_length=5,
                 poly_n=4):
        """
        初始化引物设计参数配置
        
        Parameters:
        -----------
        min_length : int, optional
            最小引物长度 (default: 15)
        max_length : int, optional
            最大引物长度 (default: 20)
        gc_min : float, optional
            最小GC含量百分比 (default: 40.0)
        gc_max : float, optional
            最大GC含量百分比 (default: 60.0)
        tm_min : float, optional
            最小熔解温度 (default: 47.0)
        tm_max : float, optional
            最大熔解温度 (default: 53.0)
        min_gap : int, optional
            引物之间的最小间隔 (default: 2)
        r_m_gap : int, optional
            R探针和M探针之间的固定间隔 (default: 2)
        min_complementary_length : int, optional
            判定为互补的最小连续碱基长度 (default: 4)
        poly_n : int, optional
            连续相同碱基的最大允许数量 (default: 4)
        """
        self.min_length = min_length
        self.max_length = max_length
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.min_gap = min_gap  # 引物之间的最小间隔
        self.r_m_gap = r_m_gap  # R探针和M探针之间的固定间隔
        self.min_complementary_length = min_complementary_length
        self.poly_n = poly_n

def find_valid_primers(sequence, start_pos, end_pos, config, used_positions):
    """
    在指定区域内寻找符合条件的引物
    
    Parameters:
    -----------
    sequence : str
        目标序列
    start_pos : int
        起始位置
    end_pos : int
        终止位置
    config : PrimerDesignConfig
        引物设计参数配置
    used_positions : set
        已使用的位置集合
        
    Returns:
    --------
    list: 符合条件的引物列表，每个元素为(start, end, sequence)元组
    
    Example:
        >>> find_valid_primers("ATGC", 0, 4, PrimerDesignConfig(), set())
        [(0, 4, 'ATGC')]
    """
    valid_primers = []
    
    for start in range(start_pos, end_pos - config.min_length + 1):
        for length in range(config.min_length, config.max_length + 1):
            end = start + length
            if end > end_pos:
                continue
                
            # 检查是否与已使用区域重叠
            if any(pos in used_positions for pos in range(start, end)):
                continue
                
            primer_seq = sequence[start:end]
            if is_valid_primer(primer_seq,
                             config.min_length,
                             config.max_length,
                             config.gc_min,
                             config.gc_max,
                             config.tm_min,
                             config.tm_max,
                             poly_n=config.poly_n):
                tm = mt.Tm_NN(primer_seq, nn_table=mt.DNA_NN4)
                gc = calculate_gc_content(primer_seq)
                valid_primers.append((start, end, primer_seq, tm, gc))
                
    return valid_primers

def design_prime_set(sequence, config):
    """
    使用贪婪算法设计探针组合。
    
    算法步骤：
    1. 从左到右扫描序列，找到第一个符合条件的左引物
    2. 在左引物后找到第一个符合条件的中间引物
    3. 在中间引物后找到第一个符合条件的右引物
    4. 如果找到完整的组合，保存并继续从下一个位置寻找
    
    Parameters:
    -----------
    sequence : str
        目标序列
    config : PrimerDesignConfig
        引物设计参数配置
        
    Returns:
    --------
    list: 符合条件的探针组合列表，每个元素为包含三个引物信息的元组，
          每个引物信息包含：(序列, 起始位置, 结束位置, Tm值, GC含量)
    """
    sequence_length = len(sequence)
    primer_sets = []
    used_positions = set()
    
    def find_next_valid_primer(start_pos, end_pos, used_pos, must_start_at=None):
        """在指定范围内找到第一个有效的引物
        
        Parameters:
        -----------
        start_pos : int
            起始位置
        end_pos : int
            结束位置
        used_pos : set
            已使用的位置集合
        must_start_at : int, optional
            必须从这个位置开始（用于确保Middle紧接Left）
        """
        # 如果指定了必须的起始位置，就只检查从这个位置开始的引物
        if must_start_at is not None:
            start_pos = must_start_at
            
        for pos in range(start_pos, end_pos - config.min_length + 1):
            # 如果指定了必须的起始位置，但当前位置不是该位置，则跳过
            if must_start_at is not None and pos != must_start_at:
                break
                
            # 检查位置是否已被使用
            if any(p in used_pos for p in range(pos, pos + config.max_length)):
                continue
                
            # 尝试不同长度的引物
            for length in range(config.min_length, min(config.max_length + 1, end_pos - pos + 1)):
                primer = sequence[pos:pos + length]
                try:
                    tm = mt.Tm_NN(primer, nn_table=mt.DNA_NN4)
                    gc = calculate_gc_content(primer)
                    if (gc >= config.gc_min and gc <= config.gc_max and 
                        tm >= config.tm_min and tm <= config.tm_max and
                        check_poly_n(primer, config.poly_n) and not has_hairpin(primer)):
                        # 在返回之前进行反向互补
                        rev_comp_primer = str(Seq(primer).reverse_complement())
                        logger.debug(f"Found valid primer: {rev_comp_primer}, Tm: {tm:.1f}, GC: {gc:.1f}")
                        return (rev_comp_primer, pos, pos + length, tm, gc)
                except Exception as e:
                    logger.error(f"计算引物参数时出错: {e}")
                    continue
            return None
    
    def check_primer_compatibility(primer1, primer2, primer3=None):
        """检查引物之间是否存在互补性问题"""
        primers = [p[0] for p in [primer1, primer2, primer3] if p is not None]
        return not has_dimer_issues(primers, config.min_complementary_length)
    
    # 从左到右扫描序列
    pos = 0
    while pos < sequence_length - 3 * config.min_length:
        logger.debug(f"Searching for primer set starting at position {pos}")
        
        # 寻找左引物
        left_result = find_next_valid_primer(pos, sequence_length, used_positions)
        if not left_result:
            pos += 1
            continue
            
        left_primer, left_start, left_end, left_tm, left_gc = left_result
        
        # 寻找中间引物（必须紧接着左引物）
        mid_result = find_next_valid_primer(left_end, sequence_length, used_positions, must_start_at=left_end)
        if not mid_result:
            pos += 1
            continue
            
        mid_primer, mid_start, mid_end, mid_tm, mid_gc = mid_result
        
        # 检查左引物和中间引物的兼容性
        if not check_primer_compatibility(left_result, mid_result):
            pos += 1
            continue
        
        # 寻找右引物（与中间引物之间保持固定的gap）
        right_start = mid_end + config.r_m_gap  # 使用配置中的固定gap值
        if right_start >= sequence_length - config.min_length:
            pos += 1
            continue
            
        right_result = find_next_valid_primer(right_start, sequence_length, used_positions)
        if not right_result:
            pos += 1
            continue
            
        right_primer, right_start, right_end, right_tm, right_gc = right_result
        
        # 检查三个引物的兼容性
        if check_primer_compatibility(left_result, mid_result, right_result):
            logger.info(f"Found valid primer set:")
            logger.info(f"Left primer: {left_primer}, Tm: {left_tm:.1f}")
            logger.info(f"Middle primer: {mid_primer}, Tm: {mid_tm:.1f}")
            logger.info(f"Right primer: {right_primer}, Tm: {right_tm:.1f}")
            primer_sets.append((left_result, mid_result, right_result))
            
            # 更新已使用的位置
            for i in range(left_start, left_end):
                used_positions.add(i)
            for i in range(mid_start, mid_end):
                used_positions.add(i)
            for i in range(right_start, right_end):
                used_positions.add(i)
            
            # 从左引物后开始继续搜索
            pos = left_end
        else:
            pos += 1
    
    return primer_sets

def create_triplet_probe(primer_sets, bridge_probe):
    """
    创建三合一探针, 用于TCR

    Parameters:
    -----------
    primer_sets : list
        探针组合列表，每个元素包含三个探针信息 (L, M, R)
    bridge_probe : str
        桥接探针序列，长度必须为19个碱基

    Returns:
    --------
    list: 转换后的三合一探针列表，每个元素包含：
        L: L + N + brigde_probe[0:16] + (bridge_probe[17:19]的互补序列) + AAGATA
        M: ACATTA + M
        R: R + TAATGTTATCTT
    """
    if not bridge_probe or len(bridge_probe) != 19:
        raise ValueError("桥接探针必须为19个碱基长度")

    # 定义碱基互补对应关系
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # 随机碱基选择
    random_base = random.choice(['A', 'T', 'C', 'G'])
    
    # 获取bridge_probe最后两个碱基的互补序列
    bridge_end_complement = ''.join(complement[base] for base in bridge_probe[17:19])
    
    triplet_probes = []
    for left, middle, right in primer_sets:
        # 使用primer_sets中的序列信息
        L = left[0]  # 获取左引物序列
        M = middle[0]  # 获取中间引物序列
        R = right[0]  # 获取右引物序列
        
        # 构建三个探针
        L_probe = L + random_base + bridge_probe[0:17] + bridge_end_complement + "AAGATA"
        M_probe = "ACATTA" + M
        R_probe = R + "TAATGTTATCTT"
        
        triplet_probes.append((L_probe, M_probe, R_probe))
    
    return triplet_probes

def save_primer_sets_with_blast(primer_sets, output_file, blast_db=None):
    """
    将引物组合及其BLAST分析结果保存到文件
    
    Parameters:
    -----------
    primer_sets : list
        引物组合列表
    output_file : str
        输出文件名
    blast_db : str, optional
        BLAST数据库路径
    """
    with open(output_file, 'w') as f:
        for i, (left, middle, right) in enumerate(primer_sets, 1):
            f.write(f"\n引物组合 {i}:\n")
            f.write("-" * 80 + "\n")
            
            for primer_type, primer in [("Left", left), ("Middle", middle), ("Right", right)]:
                sequence = primer[0]
                f.write(f"\n{primer_type} Primer: {sequence}\n")
                f.write(f"Length: {len(sequence)}bp\n")
                f.write(f"GC Content: {calculate_gc_content(sequence):.1f}%\n")
                f.write(f"Tm: {mt.Tm_NN(sequence, nn_table=mt.DNA_NN4):.1f}°C\n")
                
                if blast_db:
                    mismatch_counts, detailed_matches = analyze_blast_results(sequence, blast_db)
                    f.write("\nBLAST Analysis:\n")
                    f.write("Mismatches\tCount\tDescription\n")
                    f.write("-" * 50 + "\n")
                    for mismatches, count in sorted(mismatch_counts.items()):
                        description = "Perfect match" if mismatches == 0 else f"{mismatches} mismatch(es)"
                        f.write(f"{mismatches}\t\t{count}\t\t{description}\n")
                    
                    # 输出详细的匹配信息
                    if detailed_matches:
                        f.write("\nDetailed Matches:\n")
                        f.write("Subject ID\tMismatches\tGaps\tTotal Mismatches\n")
                        f.write("-" * 70 + "\n")
                        for subject_id, mismatches, gaps in detailed_matches:
                            total = mismatches + gaps
                            f.write(f"{subject_id}\t{mismatches}\t{gaps}\t{total}\n")
                
                f.write("\n")

def save_triplet_probes_to_csv(triplet_probes, output_file, task_name, BP_ID, delimiter=','):
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

def analyze_blast_results(primer, db_path):
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

def create_triplet_probe(primer_sets, bridge_probe):
    """
    创建三合一探针, 用于TCR

    Parameters:
    -----------
    primer_sets : list
        探针组合列表，每个元素包含三个探针信息 (L, M, R)
    bridge_probe : str
        桥接探针序列，长度必须为19个碱基

    Returns:
    --------
    list: 转换后的三合一探针列表，每个元素包含：
        L: L + N + brigde_probe[0:16] + (bridge_probe[17:19]的互补序列) + AAGATA
        M: ACATTA + M
        R: R + TAATGTTATCTT
    """
    if not bridge_probe or len(bridge_probe) != 19:
        raise ValueError("桥接探针必须为19个碱基长度")

    # 定义碱基互补对应关系
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # 随机碱基选择
    random_base = random.choice(['A', 'T', 'C', 'G'])
    
    # 获取bridge_probe最后两个碱基的互补序列
    bridge_end_complement = ''.join(complement[base] for base in bridge_probe[17:19])
    
    triplet_probes = []
    for left, middle, right in primer_sets:
        # 使用primer_sets中的序列信息
        L = left[0]  # 获取左引物序列
        M = middle[0]  # 获取中间引物序列
        R = right[0]  # 获取右引物序列
        
        # 构建三个探针
        L_probe = L + random_base + bridge_probe[0:17] + bridge_end_complement + "AAGATA"
        M_probe = "ACATTA" + M
        R_probe = R + "TAATGTTATCTT"
        
        triplet_probes.append((L_probe, M_probe, R_probe))
    
    return triplet_probes

def save_primer_positions_bed(primer_sets, output_file, task_name):
    """
    将引物位置信息保存为BED格式
    
    Parameters:
    -----------
    primer_sets : list
        引物组合列表，每个元素包含三个引物信息 (L, M, R)，每个引物信息包含 (序列, 起始位置, 结束位置, ...)
    output_file : str
        输出文件名
    task_name : str
        任务名称，用作染色体名
    """
    # 从output_file获取基础文件名（不包含扩展名）
    base_name = output_file.rsplit('.', 1)[0]
    bed_file = f"{base_name}_positions.bed"
    
    with open(bed_file, 'w') as f:
        for i, (left, middle, right) in enumerate(primer_sets, 1):
            # 获取每个引物的位置信息
            _, l_start, l_end, _, _ = left
            _, m_start, m_end, _, _ = middle
            _, r_start, r_end, _, _ = right
            
            # 写入左引物位置
            f.write(f"{task_name}\t{l_start}\t{l_end}\tL-{i}\n")
            # 写入中间引物位置
            f.write(f"{task_name}\t{m_start}\t{m_end}\tM-{i}\n")
            # 写入右引物位置
            f.write(f"{task_name}\t{r_start}\t{r_end}\tR-{i}\n")

def main(sequence=None, config=None, blast_db=None, output_file="primer_set.txt",
         bridge_probe=None, BP_ID=None, task_name=None):
    """
    主函数
    
    Parameters:
    -----------
    sequence : str, optional
        目标序列，如果为None则使用默认序列
    config : PrimerDesignConfig, optional
        引物设计参数配置，如果为None则使用默认配置
    blast_db : str, optional
        BLAST数据库路径
    output_file : str, optional
        输出文件名 (default: "primer_set.txt")
    bridge_probe : str, optional
        桥接探针序列，用于生成三合一探针
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
        # 设计引物组合
        primer_sets = design_prime_set(sequence, config)
        logger.info(f"设计了 {len(primer_sets)} 个引物组合")
        
        # 保存引物组合
        save_primer_sets_with_blast(primer_sets, output_file, blast_db)
        
        # 保存引物位置信息为BED格式
        save_primer_positions_bed(primer_sets, output_file, task_name)
        
        if bridge_probe:
            try:
                triplet_probes = create_triplet_probe(primer_sets, bridge_probe)
                logger.info(f"生成了 {len(triplet_probes)} 个三合一探针组合")
                save_triplet_probes_to_csv(triplet_probes, output_file, task_name, BP_ID)
            except ValueError as e:
                logger.error(f"生成三合一探针时出错: {str(e)}")
                
    except Exception as e:
        logger.error(f"运行过程中出错: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='HCR引物设计工具')
    # 必需参数
    parser.add_argument('--name', required=True, help='任务名称，将用作探针ID的前缀')
    parser.add_argument('--seq', required=True, help='目标序列')
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
    parser.add_argument('--output', default='primer_set.txt', help='输出文件名 (default: primer_set.txt)')
    parser.add_argument('--bridge-probe', help='桥接探针序列，用于生成三合一探针')
    
    args = parser.parse_args()
    
    # 创建配置对象
    config = PrimerDesignConfig(
        min_length=args.min_length,
        max_length=args.max_length,
        gc_min=args.min_gc,
        gc_max=args.max_gc,
        tm_min=args.min_tm,
        tm_max=args.max_tm,
        min_gap=args.min_gap
    )
    
    # 运行主函数
    main(sequence=args.seq, 
         config=config, 
         blast_db=args.ref_genome, 
         output_file=args.output, 
         bridge_probe=args.bridge_probe,
         task_name=args.name)