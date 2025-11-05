from typing import Dict, List
import logging
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from collections import Counter

logger = logging.getLogger(__name__)


def calculate_kmer_count(sequence: str, kmer_size: int = 4, min_count: int = 2) -> set:
    """计算序列中k-mer的计数"""
    def reverse_complement(seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complement[base] for base in reversed(seq)])
    
    seq_rc = reverse_complement(sequence)
    
    # 生成k-mer计数
    seq_kmers = Counter([sequence[i:i+kmer_size] for i in range(len(sequence) - kmer_size + 1)])
    rc_kmers = Counter([seq_rc[i:i+kmer_size] for i in range(len(seq_rc) - kmer_size + 1)])
    
    # 合并重复k-mer
    repeat_kmers = set()
    for kmer, count in seq_kmers.items():
        if count >= min_count:
            repeat_kmers.add(kmer)
    for kmer, count in rc_kmers.items():
        if count >= min_count:
            repeat_kmers.add(kmer)

    return repeat_kmers


def has_low_complexity(sequence: str, kmer_list: List[str]) -> bool:
    """
    分析sequence中是否存在低复杂度序列

    如果sequence中存在kmer_list中的k-mer, 则认为sequence包含低复杂度序列
    
    Parameters:
    -----------
    sequence: str
        序列
    kmer_list: List[str]
        k-mer列表, 记录重复的k-mer
    """
    for kmer in kmer_list:
        if kmer in sequence:
            return True
    return False

def calculate_tm(sequence: str) -> float:
    """计算序列的Tm值"""
    return mt.Tm_NN(sequence, nn_table=mt.DNA_NN4)

def calculate_gc_content(sequence: str) -> float:
    """计算序列的GC含量"""
    return gc_fraction(Seq(sequence)) * 100

def check_poly_n(sequence: str, poly_n: int = 4) -> bool:
    """检查序列中是否存在连续的N个相同碱基"""
    seq_str = str(sequence)
    bases = ['A', 'T', 'G', 'C']
    for base in bases:
        # 检查是否有大于或等于poly_n个连续碱基
        if seq_str.find(base * poly_n) != -1:
            return False
    return True


def is_complementary(seq1, seq2):
    """检查两个序列是否互补
    
    ATGC 和 TACG 是互补的

    Parameters:
    -----------
    seq1 : str
        序列1
    seq2 : str
        序列2

    Returns:
    --------
    bool: 如果序列互补则返回True，否则返回False
    """
    if len(seq1) != len(seq2):
        return False
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for b1, b2 in zip(seq1, seq2):
        if b2 != pairs.get(b1):
            return False
    return True

def has_hairpin(seq: str, min_stem_length: int = 4) -> bool:
    """检查序列是否可能形成发卡结构"""
    seq_length = len(seq)
    # Only check if sequence is long enough to form hairpin with given stem length
    if seq_length < 2 * min_stem_length:
        return False
        
    # Modified loop ranges to properly check all possible stem positions
    for i in range(seq_length - 2 * min_stem_length + 1):
        for j in range(i + min_stem_length, seq_length - min_stem_length + 1):
            stem1 = seq[i:i+min_stem_length]
            stem2 = seq[j:j+min_stem_length]
            stem2_reversed = stem2[::-1]  # Reverse the second stem
            if is_complementary(stem1, stem2_reversed):
                return True
    return False

def check_complementarity(seq1: str, seq2: str, min_complementary_length: int = 4) -> Dict[str, bool]:
    """检查两个序列序列之间的互补性"""
    len1, len2 = len(seq1), len(seq2)
    
    # 检查局部互补性
    for i in range(len1 - min_complementary_length + 1):
        for j in range(len2 - min_complementary_length + 1):
            if is_complementary(seq1[i:i+min_complementary_length], 
                              seq2[j:j+min_complementary_length][::-1]):
                return {'has_complementarity': True, 'end_complementarity': False}
    
    # 检查3'端互补性
    end_length = min(5, min_complementary_length)
    end_complementarity = is_complementary(seq1[-end_length:], seq2[-end_length:][::-1])
    
    return {
        'has_complementarity': False,
        'end_complementarity': end_complementarity
    }

def has_dimer_issues(seqs: List[str], min_complementary_length: int = 4) -> bool:
    """
    检查一组引物是否存在二聚体问题
    
    Parameters:
    -----------
    seqs : list
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
    if not seqs:
        return False
        
    # 检查所有可能的引物对
    for i in range(len(seqs)):
        # 检查自身互补
        if len(seqs[i]) >= min_complementary_length * 2:
            result = check_complementarity(seqs[i][:len(seqs[i])//2], 
                                        seqs[i][len(seqs[i])//2:], 
                                        min_complementary_length)
            if result['has_complementarity']:
                return True
            
        # 检查与其他引物的互补
        for j in range(i + 1, len(seqs)):
            result = check_complementarity(seqs[i], seqs[j], min_complementary_length)
            if result['has_complementarity']:
                return True
                
    return False


def is_valid_probe(sequence: str,
        min_length=17, max_length=20, 
        gc_min=40.0, gc_max=60.0, 
        tm_min=47.0, tm_max=53.0, 
        poly_n=4):
    """
    检查引物是否符合要求
    
    Parameters:
    -----------
    sequence : str
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
    length = len(sequence)
    if length < min_length or length > max_length:
        logger.debug(f"长度不在范围内: {length}")
        return False

    gc_content = calculate_gc_content(sequence)
    if gc_content < gc_min or gc_content > gc_max:
        logger.debug(f"GC含量不在范围内: {gc_content}")
        return False

    try:
        tm = calculate_tm(sequence)
        if tm < tm_min or tm > tm_max:
            logger.debug(f"熔点温度不在范围内: {tm}")
            return False
    except Exception as e:
        logger.error(f"计算熔点温度时出错: {e}")
        return False

    # 避免连续N个相同的碱基
    if not check_poly_n(sequence, poly_n):
        logger.debug(f"连续N个相同的碱基: {sequence}")
        return False


    return True


def check_left_probe_gc_junction(left_probe_seq: str) -> bool:
    """
    检查Left探针连接处是否为G/C碱基（SplintR连接酶抑制位点）

    在使用SplintR连接酶的双探针系统中，供体侧（Left探针3'端）连接处
    如果是dG/C或dC/G碱基对，会部分抑制SplintR酶活性，降低连接效率。
    因此应优先选择A/T碱基的探针以获得更好的连接效果。

    技术说明：
    - SplintR连接酶（PBCV-1 DNA连接酶）可催化以互补RNA链为模板的单链DNA连接
    - 该酶对连接处碱基对有选择性：dC/G和dG/C配对会降低酶活性
    - Left探针序列经过reverse_complement转换后，序列第一个碱基对应连接处

    Parameters:
    -----------
    left_probe_seq : str
        Left探针序列（已经过reverse_complement转换）

    Returns:
    --------
    bool: 如果连接处为G或C碱基则返回True（表示会抑制SplintR酶活性，不利于连接）

    Example:
        >>> check_left_probe_gc_junction("GGCATGCATG")  # G开头，不利于连接
        True
        >>> check_left_probe_gc_junction("CGCATGCATG")  # C开头，不利于连接
        True
        >>> check_left_probe_gc_junction("AGCATGCATG")  # A开头，有利于连接
        False
        >>> check_left_probe_gc_junction("TGCATGCATG")  # T开头，有利于连接
        False
    """
    if len(left_probe_seq) < 1:
        logger.warning(f"Left探针序列为空: {left_probe_seq}")
        return False

    # 检查第一个碱基是否是G或C（会抑制SplintR酶活性）
    return left_probe_seq[0] in 'GC'
