from typing import List, Dict
import numpy as np
import pandas as pd
import io
import tempfile
from multiprocessing import Pool

from Bio.Blast.Applications import NcbiblastnCommandline as bn
from Bio.SeqUtils import MeltingTemp as mt

import logging

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

def blastn(seq, blastn_db, blastn_evalue=1e-1, blastn_word_size=7, blastn_num_threads=4):
    """
    :param seq: 欲blast的序列
    :param blastn_db: blastn数据库
    :param blastn_evalue: evalue
    :param blastn_word_size: word_size
    :param blastn_num_threads: 线程数
    """

    tmp_fasta = ""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        # 如果seq是字符串
        if isinstance(seq, str):
            f.write('>{}\n{}\n'.format("seq", seq))
        # 如果seq是list
        elif isinstance(seq, list):
            for i, s in enumerate(seq):
                f.write('>{}\n{}\n'.format(i, s))
        # 如果seq是字典
        elif isinstance(seq, dict):
            for k, v in seq.items():
                f.write('>{}\n{}\n'.format(k, v))
        else:
            raise TypeError("seq must be str, list or dict")
        
        tmp_fasta = f.name
    
    cline = bn(
            #cmd = "/opt/biosoft/ncbi-blast-2.10.1+/bin/blastn",
            query= tmp_fasta,
            subject= blastn_db,
            outfmt=6,
            task="blastn-short",
            evalue=blastn_evalue,
            word_size=blastn_word_size,
            num_threads=blastn_num_threads,

        )  # this uses biopython's blastn formatting function and creates a commandline compatible command
    (
        stdout,
        _,
    ) = (
        cline()
    )  # cline() calls the string as a command and passes it to the command line, outputting the blast results to one variable and errors to the other

    ## From results of blast creating a numpy array (and Pandas database)
    try:
        blastresult = pd.read_csv(io.StringIO(stdout), delimiter="\t", header=None)
        blastresult.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        # 增加qlen, plen
        blastresult["qlen"] = np.abs(blastresult["qend"] - blastresult["qstart"] ) + 1
        blastresult["plen"] = np.abs(blastresult["send"] - blastresult["sstart"] ) + 1
    except pd.errors.EmptyDataError:
        blastresult = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore","qlen","plen" ])
    
    return blastresult

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



def filter_blastn_result(blastn:pd.DataFrame, seq_chr:str, seq_start:int, seq_end:int):
    """ 删除blastn 在seq_chr:seq_start-seq_end区间的结果
    :param blastn: blastn结果
    :param seq_chr: 原序列在基因组的染色体编号
    :param seq_start: 原序列在基因组的起始位置
    :param seq_end: 原序列在基因组的终止位置
    """

    if blastn.empty:
        return blastn

    # 反选在区间的比对结果
    blastn = blastn[ ~ (blastn['sseqid'] == seq_chr) & (blastn['sstart'] >= seq_start) & (blastn['send'] <= seq_end) ]

    return blastn

def worker(args):
    pos, probe, blastdb, seq_chr, seq_start, seq_end = args
    odd_probe, even_probe = generate_odd_even_probes(probe)

    odd_probe_blastn = blastn(odd_probe, blastdb)
    even_probe_blastn = blastn(even_probe, blastdb)

    odd_probe_blastn = filter_blastn_result(odd_probe_blastn, seq_chr, seq_start, seq_end)
    even_probe_blastn = filter_blastn_result(even_probe_blastn, seq_chr, seq_start, seq_end)

    if odd_probe_blastn.empty and even_probe_blastn.empty:
        return (pos, probe)
    else:
        print(f"odd_probe_blastn: {odd_probe_blastn}")
        print(f"even_probe_blastn: {even_probe_blastn}")
        return ()
    
def filter_probe_by_blastn(probes: Dict, blastdb:str, seq_chr=None, seq_start=None, seq_end=None) -> Dict[int,str]:
    """基于BLASTN结果过滤探针
    为了验证探针的特异性, 我们需要先排除探针在目标基因组的原始区间上(seq_chr:seq_start-seq_end)的比对结果, 然后再检查该探针是否在其他区间存在比对结果

    :param probes: 探针字典 如{0:'aaaaaaa', }
    :param blastdb: BLASTN数据库
    :param seq_chr: 原序列在基因组的染色体编号
    :param seq_start: 原序列在基因组的起始位置
    :param seq_end: 原序列在基因组的结束位置
    """

    pool = Pool(processes=10) # Creates a pool of process, controls worksers, 10 workers in this case
    # The second argument is a list of arguments for 'worker'
    results = pool.map(worker, [(pos, probe, blastdb, seq_chr, seq_start, seq_end) for pos, probe in probes.items()])  
    pool.close()  # Close pool
    pool.join()  # Wait for all workers to finish
    filtered_probes = {pos: probe for result in results if len(result) == 2 for pos, probe in [result]}
    
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
                  min_tm:int=45, max_tm:int=55, blastdb = None):
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
    logging.info('当前探针数目：%d', len(probes))

    # 基于GC含量过滤探针
    probes = filter_probe_by_gc(probes,  min_gc , max_gc, inner_gap, odd_probe_size, even_probe_size)
    logging.info('当前探针数目：%d', len(probes))

    #  基于Tm值过滤探针
    probes = filter_probe_py_tm(probes,  min_tm, max_tm, inner_gap, odd_probe_size, even_probe_size)
    logging.info('当前探针数目：%d', len(probes))

    # 筛选互不重叠的探针
    probes = select_maxium_probes(probes, min_gap=min_gap, method= 'quick')
    logging.info('当前探针数目：%d', len(probes))

    # 基于BLASTN的结果进行过滤
    # TODO: BLAST考虑的是设计探针的基因， 如果比对的是其他基因组，代码应该不一样
    if blastdb:
        
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
