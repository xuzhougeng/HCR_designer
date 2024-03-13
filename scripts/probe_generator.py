from typing import List, Dict

from Bio.SeqUtils import MeltingTemp as mt
import logging
from collections import Counter

# Configure the logging to display INFO messages
logging.basicConfig(level=logging.INFO)


def clean_sequence(sequence):
    """
    清除输入的cDNA序列中的非cDNA字符（即非ATGC字符）。

    :param cdna_sequence: 输入的cDNA序列（字符串）。
    :return: 清除非cDNA字符后的序列。
    """
    # 仅保留A, T, G, C字符
    valid_nucleotides = {"A", "T", "G", "C", "a", "t", "c", "g"}
    cleaned_sequence = ''.join([nucleotide for nucleotide in sequence if nucleotide in valid_nucleotides])

    return cleaned_sequence

def generate_all_probes(seq, probe_size):
    """生成所有的探针
    :param seq: 原始序列
    :param length: 探针长度
    """

    #探针字典，记录起始位置和探针序列
    probes = {}
    for i in range(len(seq) - probe_size + 1):
        probes[i] = seq[i: i + probe_size]

    return probes

def generate_odd_even_probes(seq, gap=0, odd_probe_size=25, even_probe_size=25):
    odd_probe = seq[:odd_probe_size]
    even_probe = seq[odd_probe_size + gap: odd_probe_size + gap + even_probe_size]
    
    return odd_probe, even_probe

def filter_probe_by_polyN(probes, polyN=5, gap=0, odd_probe_size=25, even_probe_size=25) -> List[str]:
    """根据多聚体N过滤探针
    :param probes: 探针字典
    :param polyN: 多聚体N的长度
    """
    filtered_probes = {}
    polyA = 'A' * polyN
    polyT = 'T' * polyN
    polyC = 'C' * polyN
    polyG = 'G' * polyN

    # 过滤探针字典
    for pos, probe in probes.items():
        probe = str(probe).upper()
        odd_probe, even_probe = generate_odd_even_probes(probe, gap, odd_probe_size, even_probe_size)
        
        # 当前探针中多聚体N的长度大于polyN时，则跳过
        if polyA in odd_probe or polyT in odd_probe or polyC in odd_probe or polyG in odd_probe:
            continue
        if polyA in even_probe or polyT in even_probe or polyC in even_probe or polyG in even_probe:
            continue
    
        filtered_probes[pos] = probe

    return filtered_probes

def filter_probe_by_gc(probes,gc_min=0.4, gc_max=0.5, gap=0, odd_probe_size=25, even_probe_size=25 ) -> List[str]:
    """根据GC含量过滤探针
    :param probes: 探针列表
    :param gc_min: 最小GC含量
    :param gc_max: 最大GC含量
    """
    filtered_probes = {}
    for pos, probe in probes.items():
        probe = str(probe).upper()
        odd_probe, even_probe = generate_odd_even_probes(probe, gap, odd_probe_size, even_probe_size)
        
        # 当GC含量小于gc_min或者大于gc_max时，则跳过
        if (odd_probe.count('G') + odd_probe.count('C')) / len(odd_probe) < gc_min or (odd_probe.count('G') + odd_probe.count('C')) / len(odd_probe) > gc_max:
            continue
        if (even_probe.count('G') + even_probe.count('C')) / len(even_probe) < gc_min or (even_probe.count('G') + even_probe.count('C')) / len(even_probe) > gc_max:
            continue

        filtered_probes[pos] = probe

    return filtered_probes

def filter_probe_py_tm(probes, min_tm=45, max_tm=55, gap=0, odd_probe_size=25, even_probe_size=25 ) -> List[str]:
    filtered_probes = {}
    for pos, probe in probes.items():
        probe = str(probe).upper()
        odd_probe, even_probe = generate_odd_even_probes(probe, gap, odd_probe_size, even_probe_size)
        odd_tm = mt.Tm_NN(odd_probe)
        even_tm = mt.Tm_NN(even_probe)

        if odd_tm < min_tm or odd_tm > max_tm:
            continue
        if even_tm < min_tm or even_tm > max_tm:
            continue

        filtered_probes[pos] = probe

    return filtered_probes

def filter_probe_by_kmer(probes, seq, k=10):
    
    def reverse_complement(seq):
        """计算序列的反向互补序列"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complement[base] for base in reversed(seq)])
    
    # 生成目标序列及其反向互补序列的k-mer集合
    seq_rc = reverse_complement(seq)

    seq_k_mer_counts = Counter([seq[i:i+k] for i in range(len(seq) - k + 1)])
    seq_rc_k_mer_counts = Counter([seq_rc[i:i+k] for i in range(len(seq_rc) - k + 1)])

    # 过滤出现次数大于1的k-mer及其次数
    seq_k_mer_counts_filtered = {k_mer: count for k_mer, count in seq_k_mer_counts.items() if count > 1}
    seq_rc_k_mer_counts_filtered = {k_mer: count for k_mer, count in seq_rc_k_mer_counts.items() if count > 1}

    # merge seq
    black_list = set(list(seq_k_mer_counts_filtered.keys()) + list(seq_rc_k_mer_counts_filtered.keys()))

    # 初始化过滤后的探针字典
    filtered_probes = {}
    
    # 遍历探针集
    for pos, probe in probes.items():
        probe = str(probe).upper()
        
        # 生成探针的k-mer集合
        probe_k_mers = {probe[i:i+k] for i in range(len(probe) - k + 1)}
        
        # 如果探针的k-mer集合与黑名单没有交集，则保留该探针
        if not probe_k_mers.intersection(black_list):
            filtered_probes[pos] = probe
    
    return filtered_probes



def select_maxium_probes(probes, min_gap = 2, method = "quick") -> Dict[str,str]:
    """基于probes的位置选择探针
    每个探针的位置包含[start,end], 需要从中选择互不重叠的最多的区间

    
    :param probes: 探针列表
    :param min_gap: 探针序列之间的最小间隔 
    :param method: 选择方法, 'quick'为 贪婪算法 , slow 对应 动态规划算法
    """

    pos_table = [ ]
    for pos, probe in probes.items():
        pos_table.append([pos, pos + len(probe) + min_gap])

    #print(pos_table)

    if method == 'quick':
        pos_table = greedy(pos_table)
    elif method == 'slow':
        pos_table = dp(pos_table)
    else:
        # 抛出异常
        raise Exception('method must be quick or slow')
    
    # 根据pos_tabel的start过滤probes
    filtered_probes = {}
    
    pos_start = [_[0] for _ in pos_table ]

    for pos, probe in probes.items():
        if pos in pos_start:
            filtered_probes[pos] = probe
    
    return filtered_probes

def greedy(position_list):
    """Using greedy algorithm to select the maximum number of non-overlapping probes
    :param position_list: a list of probe positions, like [[1, 12], [5, 17], [20, 32]]
    """
    if not position_list:
        return []

    # Sort the positions based on the end point
    position_list.sort(key=lambda x: x[1])
    selected_probes = [position_list[0]]

    for probe in position_list:
        # If the current probe's start position is greater than the last selected probe's end position, 
        # we need to select this probe.
        if probe[0] > selected_probes[-1][1]:
            selected_probes.append(probe)
    
    return selected_probes

def binary_search(probe_list, index):
    left, right = 0, index - 1
    while left <= right:
        mid = (left + right) // 2
        if probe_list[mid][1] < probe_list[index][0]:
            if probe_list[mid + 1][1] < probe_list[index][0]:
                left = mid + 1
            else:
                return mid
        else:
            right = mid - 1
    return -1

def dp(position_list):
    """Using dynamic programming algorithm to select the maximum number of non-overlapping probes
    :param position_list: a list of probe positions, like [[1, 12], [5, 17], [20, 32]]
    """
    if not position_list:
        return []

    position_list.sort(key=lambda x: x[1])

    n = len(position_list)
    dp = [0 for _ in range(n)]
    dp[0] = 1

    for i in range(1, n):
        included = 1 + dp[binary_search(position_list, i)] if binary_search(position_list, i) != -1 else 1
        excluded = dp[i - 1]
        dp[i] = max(included, excluded)

    i = n - 1
    selected_probes = []
    while i >= 0:
        if i == 0 or dp[i] != dp[i-1]:
            selected_probes.append(position_list[i])
            i = binary_search(position_list, i)
        else:
            i -= 1

    return selected_probes[::-1]


def create_probes(seq:str, probe_size:int, inner_gap=0, min_gap=2, 
                  polyN:int = 5, 
                  min_gc:float=0.3, max_gc:float=0.7, 
                  min_tm:int=45, max_tm:int=55, k:int = 8, blastdb = None):
    """主函数
    seq: 输入序列, 可以是cDNA序列, 也可以是cds序列, 但是不能包含非ATGC字符

    :param probe_size: 探针长度, 默认为50, 包括inter_gap

    :param inner_gap: 探针内部引物的间隔
    :param min_gap: 探针之间的最小间隔

    """

    # probe_size = left_probe_size + gap_size + right_probe_size
    seq = clean_sequence(seq)

    odd_probe_size = (probe_size - inner_gap) // 2
    even_probe_size = probe_size - inner_gap - odd_probe_size

    # 生成所有的探针
    probes = generate_all_probes(seq, probe_size)
    logging.info('当前探针数目：%d', len(probes))

    # 基于多聚体N过滤探针
    probes = filter_probe_by_polyN(probes, polyN, inner_gap, odd_probe_size, even_probe_size)
    logging.info('多聚体N过滤, 当前探针数目：%d', len(probes))

    # 基于GC含量过滤探针
    probes = filter_probe_by_gc(probes,  min_gc , max_gc, inner_gap, odd_probe_size, even_probe_size)
    logging.info('GC含量过滤, 当前探针数目：%d', len(probes))

    #  基于Tm值过滤探针
    probes = filter_probe_py_tm(probes,  min_tm, max_tm, inner_gap, odd_probe_size, even_probe_size)
    logging.info('Tm值过滤, 当前探针数目：%d', len(probes))

    # 过滤低复杂度的探针
    probes = filter_probe_by_kmer(probes, seq, k=k)
    logging.info('低复杂度, 当前探针数目：%d', len(probes))

    # 筛选互不重叠的探针
    probes = select_maxium_probes(probes, min_gap=min_gap, method= 'quick')
    logging.info('互不重叠, 当前探针数目：%d', len(probes))

    # 基于BLASTN的结果进行过滤
    # TODO: BLAST考虑的是设计探针的基因， 如果比对的是其他基因组，代码应该不一样
    if blastdb:
        from .filter import blastn, filter_probe_by_blastn
        seq_blast_table = blastn(str(seq), blastdb)
        # 筛选seq_balst_table: 100% identity, 且qlen == slen, 长度大于探针
        filtered_seq_blast_table = seq_blast_table[(seq_blast_table['pident'] == 100) & (seq_blast_table['qlen'] == seq_blast_table['plen'] ) & (seq_blast_table['qlen'] >  probe_size)]
        if filtered_seq_blast_table['sseqid'].nunique() == 1:
            
            seq_chr = filtered_seq_blast_table['sseqid'].unique()[0]
            seq_start = min(filtered_seq_blast_table['sstart'].min(), filtered_seq_blast_table['send'].min())
            seq_end = max(filtered_seq_blast_table['sstart'].max(), filtered_seq_blast_table['send'].max())
            
            # 过滤探针对
            filter_probe_by_blastn(probes, blastdb, seq_chr, seq_start, seq_end)
        else:
            print("函数没写好, 暂时就不过滤了")
    
    print('当前探针数目：', len(probes))
    return probes
