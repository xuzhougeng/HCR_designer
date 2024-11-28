"""
HCR (Hybridization Chain Reaction) 引物设计工具

该模块实现了一个用于HCR引物设计的自动化工具。HCR是一种信号放大技术，
通过设计特定的引物序列，可以在目标序列存在的情况下触发链式反应，从而实现信号放大。

工作流程:
1. 输入处理
   - 接收目标DNA序列
   - 配置引物设计参数（长度、GC含量、Tm值等）

2. 引物候选位置搜索
   - 扫描目标序列的所有可能位置
   - 对每个位置生成不同长度的候选引物
   - 初步筛选符合基本要求的引物（长度、GC含量等）

3. 引物质量评估
   - 检查每个候选引物的以下特性：
     * GC含量是否在合适范围
     * 熔解温度（Tm）是否适合
     * 是否存在连续重复碱基
     * 是否存在自身互补性
   - 评估引物间的互补性：
     * 全长互补性检查
     * 3'端互补性检查
     * 部分序列互补性检查

4. 探针组合设计
   - 为每个目标位置设计三个引物（左、中、右）
   - 确保引物之间：
     * 位置合适（保持最小间隔）
     * 不存在显著互补性
     * 满足所有设计参数要求

5. 结果输出
   - 生成满足所有条件的探针组合
   - 将结果写入文件，包含：
     * 左引物序列
     * 中间引物序列
     * 右引物序列

主要功能:
- 引物有效性检查（is_valid_primer）
- 碱基互补性检查（is_complementary）
- 序列互补性分析（check_complementarity）
- 二聚体问题检测（has_dimer_issues）
- 有效引物搜索（find_valid_primers）
- 探针组合设计（design_probe_set）
- 结果输出（output_probe）

使用方法:
1. 直接运行脚本，使用默认参数：
   python tcr2.py

2. 在其他代码中导入使用：
   from tcr2 import design_probe_set, PrimerDesignConfig
   config = PrimerDesignConfig()
   probes = design_probe_set(sequence, config)

注意事项:
1. 输入序列应为有效的DNA序列（ATCG）
2. 参数配置应根据实际需求调整
3. 结果质量取决于参数设置的合理性
4. 建议进行实验验证

依赖:
- Bio.Seq：用于DNA序列操作
- Bio.SeqUtils：用于序列分析工具
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

# 配置日志
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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

def is_complementary(seq1, seq2):
    """检查两个序列是否互补"""
    if len(seq1) != len(seq2):
        return False
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for b1, b2 in zip(seq1, seq2):
        if b2 != pairs.get(b1):
            return False
    return True

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

def check_poly_n(seq, n=4):
    """检查序列中是否存在连续的N个相同碱基"""
    for base in ['A', 'T', 'C', 'G']:
        if base * n in seq:
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

def design_probe_set(sequence, config):
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
    probe_sets = []
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
                        logger.debug(f"Found valid primer: {primer}, Tm: {tm:.1f}, GC: {gc:.1f}")
                        return (primer, pos, pos + length, tm, gc)
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
        logger.debug(f"Searching for probe set starting at position {pos}")
        
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
            logger.info(f"Found valid probe set:")
            logger.info(f"Left primer: {left_primer}, Tm: {left_tm:.1f}")
            logger.info(f"Middle primer: {mid_primer}, Tm: {mid_tm:.1f}")
            logger.info(f"Right primer: {right_primer}, Tm: {right_tm:.1f}")
            probe_sets.append((left_result, mid_result, right_result))
            
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
    
    return probe_sets

def output_probe(probe_sets, output_file, blast_db=None, delimiter=','):
    """
    将探针组合输出到文件
    
    Parameters:
    -----------
    probe_sets : list
        探针组合列表
    output_file : str
        输出文件名
    blast_db : str, optional
        BLAST数据库路径
    delimiter : str, optional
        CSV文件分隔符 (default: ',')
    """
    # 主输出文件
    with open(output_file, 'w') as f:
        # 写入表头
        headers = ['Set', 'Type', 'Sequence', 'Length', 'Start', 'End', 'Tm(°C)', 'GC(%)',
                  'Perfect Match', '1bp Mismatch', '2bp Mismatch', '3bp Mismatch', '4bp Mismatch']
        f.write(delimiter.join(headers) + '\n')
        
        for i, (left, middle, right) in enumerate(probe_sets, 1):
            for probe_type, probe in [("Left", left), ("Middle", middle), ("Right", right)]:
                sequence, start, end = probe[0], probe[1], probe[2]
                length = end - start
                tm = probe[3]
                gc = probe[4]
                
                # 获取BLAST结果
                blast_counts = {i: 0 for i in range(5)}  # 默认值
                if blast_db:
                    blast_counts = analyze_blast_results(sequence, blast_db)
                
                # 格式化输出行
                row = [
                    str(i),
                    probe_type,
                    sequence,
                    str(length),
                    str(start),
                    str(end),
                    f"{tm:.1f}",
                    f"{gc:.1f}",
                    str(blast_counts[0]),
                    str(blast_counts[1]),
                    str(blast_counts[2]),
                    str(blast_counts[3]),
                    str(blast_counts[4])
                ]
                f.write(delimiter.join(row) + '\n')

    # 如果有BLAST数据库，创建详细的比对记录文件
    if blast_db:
        detail_file = output_file.rsplit('.', 1)[0] + '_blast_details.txt'
        with open(detail_file, 'w') as f:
            f.write("Detailed BLAST Analysis Report\n")
            f.write("=" * 80 + "\n\n")
            
            for i, (left, middle, right) in enumerate(probe_sets, 1):
                f.write(f"Probe Set {i}\n")
                f.write("-" * 40 + "\n")
                
                for probe_type, probe in [("Left", left), ("Middle", middle), ("Right", right)]:
                    sequence = probe[0]  # 获取序列
                    f.write(f"\n{probe_type} Primer: {sequence}\n")
                    f.write("Length: {}bp\n".format(len(sequence)))
                    
                    # 获取详细的BLAST结果
                    blast_records = get_detailed_blast_results(sequence, blast_db)
                    
                    if blast_records:
                        f.write("\nSignificant Matches:\n")
                        f.write("Mismatches\tSubject ID\tAlignment\n")
                        f.write("-" * 60 + "\n")
                        for record in blast_records:
                            f.write(f"{record['mismatches']}\t{record['subject_id']}\t{record['alignment']}\n")
                    else:
                        f.write("No significant matches found\n")
                    
                    f.write("\n" + "-" * 40 + "\n")
                
                f.write("\n")

def analyze_blast_results(primer, db_path):
    """
    对单个引物进行BLAST分析并统计不同错配数量的匹配数
    
    Parameters:
    -----------
    primer : str
        引物序列
    db_path : str
        BLAST数据库路径
    
    Returns:
    --------
    dict: 包含不同错配数量的统计结果
    """
    # 创建临时文件存储序列
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(f">primer\n{primer}\n")
        query_file = temp_file.name

    primer_length = len(primer)
    
    # 设置BLAST参数
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=db_path,
        task="blastn-short",
        word_size=11,
        evalue=10,
        outfmt=5,
        dust='no',
        soft_masking='false',
        perc_identity=80
    )
    
    # 运行BLAST
    stdout, stderr = blastn_cline()
    
    # 解析结果
    blast_records = NCBIXML.parse(StringIO(stdout))
    
    # 统计不同错配数量的匹配数
    mismatch_counts = {i: 0 for i in range(5)}  # 0-4个错配
    
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # 只考虑长度完全相同的匹配
                if len(hsp.query) == primer_length and len(hsp.sbjct) == primer_length:
                    # 计算错配数量（包括gaps）
                    mismatches = sum(1 for q, s in zip(hsp.query, hsp.sbjct) if q != s) + hsp.gaps
                    if mismatches <= 4:  # 只统计0-4个错配的情况
                        mismatch_counts[mismatches] += 1
    
    # 清理临时文件
    subprocess.run(['rm', query_file])
    
    return mismatch_counts

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
    # 创建临时文件存储序列
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(f">primer\n{primer}\n")
        query_file = temp_file.name

    primer_length = len(primer)
    
    # 设置BLAST参数
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=db_path,
        task="blastn-short",
        word_size=11,
        evalue=1000,
        outfmt=5,
        dust='no',
        soft_masking='false',
        perc_identity=80
    )
    
    # 运行BLAST
    stdout, stderr = blastn_cline()
    
    # 解析结果
    blast_records = NCBIXML.parse(StringIO(stdout))
    detailed_results = []
    
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # 只考虑长度完全相同的匹配
                if len(hsp.query) == primer_length and len(hsp.sbjct) == primer_length:
                    # 计算错配数量（包括gaps）
                    mismatches = sum(1 for q, s in zip(hsp.query, hsp.sbjct) if q != s) + hsp.gaps
                    if mismatches <= 4:
                        # 构建对齐显示，包括gap位置的标记
                        match_str = ''.join('|' if q == s else ' ' for q, s in zip(hsp.query, hsp.sbjct))
                        alignment_display = (
                            f"\nQuery:  {hsp.query}\n"
                            f"        {match_str}\n"
                            f"Sbjct:  {hsp.sbjct}\n"
                            f"Gaps: {hsp.gaps}, Mismatches: {mismatches-hsp.gaps}\n"
                        )
                        
                        detailed_results.append({
                            'subject_id': alignment.title,
                            'mismatches': mismatches,
                            'gaps': hsp.gaps,
                            'alignment': alignment_display,
                            'score': hsp.score,
                            'e_value': hsp.expect
                        })
    
    # 按错配数量和得分排序
    detailed_results.sort(key=lambda x: (x['mismatches'], -x['score']))
    
    # 清理临时文件
    subprocess.run(['rm', query_file])
    
    return detailed_results

def blast_probe_sets(probe_sets, db_path):
    """
    对所有探针组合进行BLAST分析
    
    Parameters:
    -----------
    probe_sets : list
        探针组合列表
    db_path : str
        BLAST数据库路径
    """
    logger.info("\nBLAST Analysis Results:")
    logger.info("-" * 80)
    
    for i, (left, middle, right) in enumerate(probe_sets, 1):
        logger.info(f"\nProbe Set {i}:")
        for probe_type, probe in [("Left", left), ("Middle", middle), ("Right", right)]:
            sequence = probe[0]  # 获取序列
            results = analyze_blast_results(sequence, db_path)
            
            logger.info(f"\n{probe_type} Primer ({sequence}, Length: {len(sequence)}bp):")
            logger.info("Mismatches\tCount\tDescription")
            logger.info("-" * 50)
            for mismatches, count in results.items():
                description = "Perfect match" if mismatches == 0 else f"{mismatches} mismatch(es)"
                logger.info(f"{mismatches}\t\t{count}\t\t{description}")



def main(sequence=None, config=None, blast_db=None):
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
    
    Example:
        >>> main("ATGC", PrimerDesignConfig())
    """
    if sequence is None:
        raise ValueError("必须提供目标序列")
    
    if config is None:
        config = PrimerDesignConfig()  # 使用默认配置
        
    probe_sets = design_probe_set(sequence, config)
    
    if len(probe_sets) == 0:
        logger.info("未找到符合条件的探针组合")
    else:
        logger.info(f"找到 {len(probe_sets)} 个符合条件的探针组合：")
        output_probe(probe_sets, "probe_set.txt", blast_db)
        
        # 如果提供了BLAST数据库，进行BLAST分析
        if blast_db:
            blast_probe_sets(probe_sets, blast_db)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='HCR引物设计工具')
    parser.add_argument('--blast_db', help='BLAST数据库路径')
    args = parser.parse_args()
    
    main(blast_db=args.blast_db)