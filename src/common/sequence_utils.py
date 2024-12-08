from typing import Dict, List
import logging
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)

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
