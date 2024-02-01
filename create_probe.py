from scripts import splint


test_seq = """ATGGAGCCGCCACAGCATCAGCATCATCATCATCAAGCCGACCAAGAAAGCGGCAACAAC
AACAACAACAAGTCCGGCTCTGGTGGTTACACGTGTCGCCAGACCAGCACGAGGTGGACA
CCGACGACGGAGCAAATCAAAATCCTCAAAGAACTTTACTACAACAATGCAATCCGGTCA
CCAACAGCCGATCAGATCCAGAAGATCACTGCAAGGCTGAGACAGTTCGGAAAGATTGAG
GGCAAGAACGTCTTTTACTGGTTCCAGAACCATAAGGCTCGTGAGCGTCAGAAGAAGAGA
TTCAACGGAACAAACATGACCACACCATCTTCATCACCCAACTCGGTTATGATGGCGGCT
AACGATCATTATCATCCTCTACTTCACCATCATCACGGTGTTCCCATGCAGAGACCTGCT
AATTCCGTCAACGTTAAACTTAACCAAGACCATCATCTCTATCATCATAACAAGCCATAT
CCCAGCTTCAATAACGGGAATTTAAATCATGCAAGCTCAGGTACTGAATGTGGTGTTGTT
AATGCTTCTAATGGCTACATGAGTAGCCATGTCTATGGATCTATGGAACAAGACTGTTCT
ATGAATTACAACAACGTAGGTGGAGGATGGGCAAACATGGATCATCATTACTCATCTGCA
CCTTACAACTTCTTCGATAGAGCAAAGCCTCTGTTTGGTCTAGAAGGTCATCAAGAAGAA
GAAGAATGTGGTGGCGATGCTTATCTGGAACATCGACGTACGCTTCCTCTCTTCCCTATG
CACGGTGAAGATCACATCAACGGTGGTAGTGGTGCCATCTGGAAGTATGGCCAATCGGAA
GTTCGCCCTTGCGCTTCTCTTGAGCTACGTCTGAACTAG"""

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


if __name__ == '__main__':
    # args = parse_args()
    # name = args.name

    # # if args.input is file, read the file
    # seq = args.input
    # probe_size = args.length # min 50
    # initiator_type = args.initiator_type
    # polyN = args.polyN
    # blastdb = args.blastdb
    # min_gc = args.min_gc
    # max_gc = args.max_gc
    # output = args.output
    # probers = create_probers(seq, probe_size, polyN, min_gc, max_gc, blastdb )

    # export(prefix = name, probers=probers, initiator_type = initiator_type , output= output )
    splint.create_primer(test_seq)