from Bio.Seq import Seq
import pandas as pd

from scripts.probe_generator import create_probes
from Bio.SeqUtils import MeltingTemp as mt


HCR_ODD_PROBE_SIZE = 25
HCR_EVEN_PROBE_SIZE = 25

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

def create_primer(seq, prefix, probe_size=50, polyN=5, 
                  min_gc=0.3, max_gc=0.7, 
                  min_tm=45, max_tm=55, 
                  initiator_type: str = "B1"):
    
    inner_gap = probe_size - HCR_EVEN_PROBE_SIZE - HCR_ODD_PROBE_SIZE
    
    probes = create_probes(seq, probe_size, inner_gap=inner_gap, polyN=polyN, min_gc=min_gc, max_gc=max_gc, min_tm=min_tm, max_tm=max_tm)
    
    # 获取initiator
    upspc, dnspc, up, dn = get_initiator_sequence(initiator_type)

    probes_pos = []
    probes_list = []
    P1_name_list = []
    P1_list = []
    P1_tm_list = []
    P2_name_list = []
    P2_list = []
    P2_tm_list = []

    middle = initiator_type.replace("B", "I")
    
    count = 1 

    for pos, probe in probes.items():
        probes_pos.append(int(pos)+1)
        probes_list.append(probe)

        P1 = Seq(probe[0:HCR_ODD_PROBE_SIZE]).reverse_complement()
        P2 = Seq(probe[-HCR_ODD_PROBE_SIZE:]).reverse_complement()
        primer_5p_tm = mt.Tm_NN(P1)
        primer_3p_tm = mt.Tm_NN(P2)


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
    
        P1_tm_list.append( primer_5p_tm )
        P2_tm_list.append( primer_3p_tm )

    probe_df = pd.DataFrame(
        {
            "probe_pos": probes_pos,
            "probe_seq": probes_list,
            "P1_name": P1_name_list,
            "P1": P1_list,
            'P1_Tm': P1_tm_list,
            "P2_name": P2_name_list,
            "P2": P2_list,
            'P2_Tm': P2_tm_list
        }
    )

    return probe_df

