from .probe_generator import create_probes
from Bio.Seq import Seq
import pandas  as pd
from Bio.SeqUtils import MeltingTemp as mt

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


# 17 bp +  2nt gap + 17 bp 
def create_primer(seq, prefix, probe_size=17, polyN=5, min_gc=0.3, max_gc=0.7, min_tm=45, max_tm=55, fulor: str = "AF488"):
    """设计splint的探针序列
    输入数据为 cds或者cdna的序列
    输出数据为splint的探针序列

    第一步: 基于csd序列设计探针 ( 分开的单一探针的Tm值, SplintPLP在50度左右 > 45 < 55)
    第二步: 获取探针
    第三步: 对探针反省互补


    """
    # prober_size = left + right for SPLINT is 34
    probes = create_probes(seq, probe_size= probe_size * 2, inner_gap=0,  polyN=polyN, min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm)

    color_seq = fluor_probe[fulor]
    probe_name_suffix = fluor_probe_name[fulor]

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

        P1_name_list.append(f"{prefix}-{count}-5{probe_name_suffix}")
        P2_name_list.append(f"{prefix}-{count}-3{probe_name_suffix}")

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

    return probe_df
    


