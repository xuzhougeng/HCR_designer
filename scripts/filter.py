
import numpy as np
import pandas as pd
import tempfile
from Bio.Blast.Applications import NcbiblastnCommandline as bn
import io
# from Bio import Align
from multiprocessing import Pool
from typing import Dict

from .probe_generator import generate_odd_even_probes

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



# # 创建一个PairwiseAligner对象 
# aligner = Align.PairwiseAligner()
# # 设置对齐模式为局部对齐
# aligner.mode = 'local'

# # 可以调整对齐的参数，比如匹配分数、不匹配的惩罚、gap的惩罚等
# aligner.match_score = 1
# aligner.mismatch_score = -1
# aligner.open_gap_score = -1
# aligner.extend_gap_score = -0.5