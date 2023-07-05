from Bio.Seq import Seq
from typing import List, Dict
from Bio.Blast.Applications import NcbiblastnCommandline as bn
import numpy as np
import pandas as pd
import io
import tempfile
from multiprocessing import Pool



HCR_ODD_PROBE_SIZE = 25
HCR_EVEN_PROBE_SIZE = 25


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

def get_initiator_sequence(ampl):
    """
    选择Initiator序列: 杂交的时候, 每一个欲检测的靶标序列对应一种Initiator, 包含
    - upspc: 选择的upstream spacer
    - dnspc: 选择的downstream spacer
    - up: 选择的upstream sequence
    - dn: 选择的downstream sequence

    :param ampl: 选择的amplification
    :return: 一个list, 包含upspc, dnspc, up, dn
    """

    if ampl == "B1":
        upspc = "aa"
        dnspc = "ta"
        up = "GAGGAGGGCAGCAAACGG"
        dn = "GAAGAGTCTTCCTTTACG"
    elif ampl == "B2":
        upspc = "aa"
        dnspc = "aa"
        up = "CCTCGTAAATCCTCATCA"
        dn = "ATCATCCAGTAAACCGCC"
    elif ampl == "B3":
        upspc = "tt"
        dnspc = "tt"
        up = "GTCCCTGCCTCTATATCT"
        dn = "CCACTCAACTTTAACCCG"
    elif ampl == "B4":
        upspc = "aa"
        dnspc = "at"
        up = "CCTCAACCTACCTCCAAC"
        dn = "TCTCACCATATTCGCTTC"
    elif ampl == "B5":
        upspc = "aa"
        dnspc = "aa"
        up = "CTCACTCCCAATCTCTAT"
        dn = "CTACCCTACAAATCCAAT"
    elif ampl == "B7":
        upspc = "ww"
        dnspc = "ww"
        up = "CTTCAACCTCCACCTACC"
        dn = "TCCAATCCCTACCCTCAC"
    elif ampl == "B9":
        upspc = "ww"
        dnspc = "ww"
        up = "CACGTATCTACTCCACTC"
        dn = "TCAGCACACTCCCAACCC"
    elif ampl == "B10":
        upspc = "ww"
        dnspc = "ww"
        up = "CCTCAAGATACTCCTCTA"
        dn = "CCTACTCGACTACCCTAG"
    elif ampl == "B11":
        upspc = "ww"
        dnspc = "ww"
        up = "CGCTTAGATATCACTCCT"
        dn = "ACGTCGACCACACTCATC"
    elif ampl == "B13":
        upspc = "ww"
        dnspc = "ww"
        up = "AGGTAACGCCTTCCTGCT"
        dn = "TTATGCTCAACATACAAC"
    elif ampl == "B14":
        upspc = "ww"
        dnspc = "ww"
        up = "AATGTCAATAGCGAGCGA"
        dn = "CCCTATATTTCTGCACAG"
    elif ampl == "B15":
        upspc = "ww"
        dnspc = "ww"
        up = "CAGATTAACACACCACAA"
        dn = "GGTATCTCGAACACTCTC"
    elif ampl == "B17":
        upspc = "ww"
        dnspc = "ww"
        up = "CGATTGTTTGTTGTGGAC"
        dn = "GCATGCTAATCGGATGAG"
    else:
        print("Please try again")
    return [upspc, dnspc, up, dn]

def generate_all_probers(seq, probe_size):
    """生成所有的探针
    :param seq: 原始序列
    :param length: 探针长度
    """

    #探针字典，记录起始位置和探针序列
    probers = {}
    for i in range(len(seq) - probe_size + 1):
        probers[i] = seq[i: i + probe_size]

    return probers

def generate_odd_even_probers(seq, gap=0):
    old_probe = seq[:HCR_ODD_PROBE_SIZE]
    even_probe = seq[HCR_ODD_PROBE_SIZE + gap: HCR_ODD_PROBE_SIZE + gap + HCR_EVEN_PROBE_SIZE]
    #print(seq)
    #print(old_probe, even_probe)
    
    return old_probe, even_probe

def filter_probe_by_polyN(probers, polyN=5) -> List[str]:
    """根据多聚体N过滤探针
    :param probers: 探针字典
    :param polyN: 多聚体N的长度
    """
    filtered_probers = {}
    polyA = 'A' * polyN
    polyT = 'T' * polyN
    polyC = 'C' * polyN
    polyG = 'G' * polyN

    # 过滤探针字典
    for pos, prober in probers.items():
        prober = str(prober).upper()
        old_probe, even_probe = generate_odd_even_probers(prober)
        
        # 当前探针中多聚体N的长度大于polyN时，则跳过
        if polyA in old_probe or polyT in old_probe or polyC in old_probe or polyG in old_probe:
            continue
        if polyA in even_probe or polyT in even_probe or polyC in even_probe or polyG in even_probe:
            continue
    
        filtered_probers[pos] = prober

    return filtered_probers

def filter_probe_by_gc(probers, gc_min=0.4, gc_max=0.5) -> List[str]:
    """根据GC含量过滤探针
    :param probers: 探针列表
    :param gc_min: 最小GC含量
    :param gc_max: 最大GC含量
    """
    filtered_probers = {}
    for pos, prober in probers.items():
        prober = str(prober).upper()
        old_probe, even_probe = generate_odd_even_probers(prober)
        
        # 当GC含量小于gc_min或者大于gc_max时，则跳过
        if (old_probe.count('G') + old_probe.count('C')) / len(old_probe) < gc_min or (old_probe.count('G') + old_probe.count('C')) / len(old_probe) > gc_max:
            continue
        if (even_probe.count('G') + even_probe.count('C')) / len(even_probe) < gc_min or (even_probe.count('G') + even_probe.count('C')) / len(even_probe) > gc_max:
            continue

        filtered_probers[pos] = prober

    return filtered_probers

def export(prefix, probers, initiator_type, output="out.csv"):
    # 获取initiator
    upspc, dnspc, up, dn = get_initiator_sequence(initiator_type)
    
    probers_pos = []
    probers_list = []
    P1_name_list = []
    P1_list = []
    P2_name_list = []
    P2_list = []

    middle = initiator_type.replace("B", "I")
    
    count = 1 

    for pos, prober in probers.items():
        probers_pos.append(int(pos)+1)
        probers_list.append(prober)

        P1 = Seq(prober[0:HCR_ODD_PROBE_SIZE]).reverse_complement()
        P2 = Seq(prober[-HCR_ODD_PROBE_SIZE:]).reverse_complement()

        P1_name_list.append(
            f"{prefix}-{middle}-{count}"
        )
        count += 1

        P1_list.append(
            up + upspc + str(P1)
        )

        P2_name_list.append(
            f"{prefix}-{middle}-{count}"
        )
        count += 1

        P2_list.append(
            str(P2) + dnspc + dn
        )
    
    prober_df = pd.DataFrame(
        {
            "probe_pos": probers_pos,
            "prober_seq": probers_list,
            "P1_name": P1_name_list,
            "P1": P1_list,
            "P2_name": P2_name_list,
            "P2": P2_list
        }
    )
    prober_df.to_csv(output, index=False)

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
    pos, prober, blastdb, seq_chr, seq_start, seq_end = args
    odd_probe, even_probe = generate_odd_even_probers(prober)

    odd_probe_blastn = blastn(odd_probe, blastdb)
    even_probe_blastn = blastn(even_probe, blastdb)

    odd_probe_blastn = filter_blastn_result(odd_probe_blastn, seq_chr, seq_start, seq_end)
    even_probe_blastn = filter_blastn_result(even_probe_blastn, seq_chr, seq_start, seq_end)

    if odd_probe_blastn.empty and even_probe_blastn.empty:
        return (pos, prober)
    else:
        print(f"odd_probe_blastn: {odd_probe_blastn}")
        print(f"even_probe_blastn: {even_probe_blastn}")
        return ()
    
def filter_prober_by_blastn(probers: Dict, blastdb:str, seq_chr=None, seq_start=None, seq_end=None) -> Dict[int,str]:
    """基于BLASTN结果过滤探针
    为了验证探针的特异性, 我们需要先排除探针在目标基因组的原始区间上(seq_chr:seq_start-seq_end)的比对结果, 然后再检查该探针是否在其他区间存在比对结果

    :param probers: 探针字典 如{0:'aaaaaaa', }
    :param blastdb: BLASTN数据库
    :param seq_chr: 原序列在基因组的染色体编号
    :param seq_start: 原序列在基因组的起始位置
    :param seq_end: 原序列在基因组的结束位置
    """



    pool = Pool(processes=10) # Creates a pool of process, controls worksers, 10 workers in this case
    # The second argument is a list of arguments for 'worker'
    results = pool.map(worker, [(pos, prober, blastdb, seq_chr, seq_start, seq_end) for pos, prober in probers.items()])  
    pool.close()  # Close pool
    pool.join()  # Wait for all workers to finish
    filtered_probers = {pos: prober for result in results if len(result) == 2 for pos, prober in [result]}
    
    return filtered_probers


def select_maxium_probers(probers, min_gap = 2, method = "quick") -> List[str]:
    """基于probers的位置选择探针
    每个探针的位置包含[start,end], 需要从中选择互不重叠的最多的区间
    
    :param probers: 探针列表
    :param method: 选择方法, 'quick'为 贪婪算法 , slow 对应 动态规划算法
    """

    pos_table = [ ]
    for pos, prober in probers.items():
        pos_table.append([pos, pos + len(prober) + min_gap])

    #print(pos_table)

    if method == 'quick':
        pos_table = greedy(pos_table)
    elif method == 'slow':
        pos_table = dp(pos_table)
    else:
        # 抛出异常
        raise Exception('method must be quick or slow')
    
    # 根据pos_tabel的start过滤probers
    filtered_probers = {}
    
    pos_start = [_[0] for _ in pos_table ]

    for pos, prober in probers.items():
        if pos in pos_start:
            filtered_probers[pos] = prober
    
    return filtered_probers

def greedy(position_list):
    """Using greedy algorithm to select the maximum number of non-overlapping probers
    :param position_list: a list of prober positions, like [[1, 12], [5, 17], [20, 32]]
    """
    if not position_list:
        return []

    # Sort the positions based on the end point
    position_list.sort(key=lambda x: x[1])
    selected_probers = [position_list[0]]

    for prober in position_list:
        # If the current prober's start position is greater than the last selected prober's end position, 
        # we need to select this prober.
        if prober[0] > selected_probers[-1][1]:
            selected_probers.append(prober)
    
    return selected_probers

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
    """Using dynamic programming algorithm to select the maximum number of non-overlapping probers
    :param position_list: a list of prober positions, like [[1, 12], [5, 17], [20, 32]]
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
    selected_probers = []
    while i >= 0:
        if i == 0 or dp[i] != dp[i-1]:
            selected_probers.append(position_list[i])
            i = binary_search(position_list, i)
        else:
            i -= 1

    return selected_probers[::-1]

def parse_args():
    """解析参数
    """
    import argparse
    parser = argparse.ArgumentParser(description='Generate all probers')
    parser.add_argument('-n', '--name', type=str, help='name of the sequence')
    parser.add_argument('-l', '--length', type=int, default=50, help='prober length')
    # initiator_type
    parser.add_argument('-t', '--initiator_type', type=str, default='B1', help='initiator type')
    parser.add_argument('-p', '--polyN', type=int, default=5, help='polyN length')
    parser.add_argument('-b', '--blastdb', type=str, help='blastdb')
    
    # GC content
    parser.add_argument('--min_gc', type=float, default=0.3, help='min GC content')
    parser.add_argument('--max_gc', type=float, default=0.5, help='max GC content')


    parser.add_argument('-i', '--input', type=str, help='input sequence')
    parser.add_argument('-o', '--output', type=str, help='output file, in csv format')
    
    args = parser.parse_args()
    return args

def main(name, seq, probe_size,initiator_type, polyN, min_gc, max_gc, output, blastdb = None):
    """主函数
    """

    # 生成所有的探针
    probers = generate_all_probers(seq, probe_size)
    print('当前探针数目：', len(probers))

    # 基于多聚体N过滤探针
    probers = filter_probe_by_polyN(probers, polyN=polyN)
    print('当前探针数目：', len(probers))

    # 基于GC含量过滤探针
    probers = filter_probe_by_gc(probers, gc_min= min_gc , gc_max= max_gc)
    print('当前探针数目：', len(probers))

    # 筛选互不重叠的探针
    probers = select_maxium_probers(probers, min_gap=2, method= 'quick')
    print('当前探针数目：', len(probers))

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
            filter_prober_by_blastn(probers, blastdb, seq_chr, seq_start, seq_end)
        else:
            print("函数没写好, 暂时就不过滤了")
    
    print('当前探针数目：', len(probers))


    export(prefix = name, probers=probers, initiator_type = initiator_type , output= output )

if __name__ == '__main__':
    args = parse_args()
    name = args.name
    seq = Seq(open(args.input).read().strip())
    probe_size = args.length # min 50
    initiator_type = args.initiator_type
    polyN = args.polyN
    blastdb = args.blastdb
    min_gc = args.min_gc
    max_gc = args.max_gc
    output = args.output
    main(name, seq, probe_size, initiator_type, polyN, min_gc, max_gc, output, blastdb )