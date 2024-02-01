# left + 2 gap + right

# 60读

from .probe_generator import create_probes
from Bio.Seq import Seq
import pandas  as pd
from Bio.SeqUtils import MeltingTemp as mt

fluor_probe = {
    "AF488": "TCGCGCTTGGTATAATCGCT",
    "Cy3": "TCGCGCTTGGTATAATCGCT",
    "Texas": "TGCGTCTATTTAGTGGAGCC",
    "Cy5": "CCTCAATGCTGCTGCTGTACTAC"    
}

# 220 bp +  2nt gap + 20 bp 
def create_primer(seq, prefix, polyN=5, min_gc=0.3, max_gc=0.7, min_tm=55, max_tm=65, fulor: str = "AF488"):
    """设计splint的探针序列
    输入数据为 cds或者cdna的序列
    输出数据为splint的探针序列

    第一步: 基于csd序列设计探针 ( 分开的单一探针的Tm值, SNIL在60度左右 > 55 < 65)
    第二步: 获取探针
    第三步: 对探针反省互补


    """
    # prober_size = left + 2 + right =  42
    probes = create_probes(seq, probe_size=42, inner_gap=2,  polyN=polyN, min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm)

    color_oligo = fluor_probe[fulor]

    oligo_5p_start = "ACATTA"
    oligo_5p_end = "AAGATA"
    
    oligo_3p = "TAATGTTATCTT"

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
        probe_5p = probe[:20].reverse_complement() 
        probe_3p = probe[-20:].reverse_complement()
        
        primer_5p_tm = mt.Tm_NN(probe_5p)
        primer_3p_tm = mt.Tm_NN(probe_3p)

        # build primer
        primer_5p = f"{oligo_5p_start}{probe_5p}{color_oligo}{oligo_5p_end}"
        primer_3p = f"{probe_3p}{oligo_3p}"

        P1_list.append(primer_5p)
        P2_list.append(primer_3p)

        P1_name_list.append(f"{prefix}-{count}-5p1")
        P2_name_list.append(f"{prefix}-{count}-3p1")

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
    


