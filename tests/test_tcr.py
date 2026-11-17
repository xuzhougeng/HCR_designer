import unittest
# scripts/tcr.py no longer exists; its primer check lives on as is_valid_probe.
from src.common.sequence_utils import is_valid_probe as is_valid_primer
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import random


def GC(sequence):
    return gc_fraction(sequence) * 100

class TestIsValidPrimer(unittest.TestCase):
    def check_primer_properties(self, seq):
        """辅助函数，打印引物的各项属性"""
        length = len(seq)
        gc = GC(seq)
        tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN4)
        print(f"Length: {length}, GC: {gc:.1f}%, Tm: {tm:.1f}°C")
        return length, gc, tm

    def generate_valid_primer(self, length=17):
        """生成一个符合所有条件的引物序列"""
        bases = ['A', 'T', 'G', 'C']
        while True:
            # 生成一个随机序列，确保GC含量适中
            seq = ''
            gc_count = 0
            for i in range(length):
                if gc_count < length * 0.4 and i >= length * 0.6:  # 需要更多GC
                    base = random.choice(['G', 'C'])
                    gc_count += 1
                elif gc_count > length * 0.6 and i < length * 0.4:  # 需要更多AT
                    base = random.choice(['A', 'T'])
                else:
                    base = random.choice(bases)
                    if base in ['G', 'C']:
                        gc_count += 1
                seq += base
            
            # 检查是否有连续重复
            has_repeats = False
            for base in bases:
                if base * 4 in seq:
                    has_repeats = True
                    break
            
            if has_repeats:
                continue
                
            # 检查GC含量和Tm值
            gc = GC(seq)
            tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN4)
            if 40 <= gc <= 60 and 47 <= tm <= 53:
                return seq

    def test_valid_primer(self):
        # 一个完全符合条件的引物
        valid_primer = self.generate_valid_primer(18)  # 18bp
        length, gc, tm = self.check_primer_properties(valid_primer)
        self.assertTrue(is_valid_primer(valid_primer))

    def test_length_constraints(self):
        # 测试长度限制
        too_short = "ATGC"  # 4bp
        too_long = "ATGCTAGCTAGCTGATCGATGC"  # 22bp
        just_right = self.generate_valid_primer(17)  # 17bp
        
        self.assertFalse(is_valid_primer(too_short))
        self.assertFalse(is_valid_primer(too_long))
        length, gc, tm = self.check_primer_properties(just_right)
        self.assertTrue(is_valid_primer(just_right))

    def test_gc_content(self):
        # 测试GC含量限制
        low_gc = "ATAATATATATATATAT"  # GC% = 0
        high_gc = "GCGCGCGCGCGCGCGCG"  # GC% = 100
        good_gc = self.generate_valid_primer(17)  # 适中的GC含量
        
        self.assertFalse(is_valid_primer(low_gc))
        self.assertFalse(is_valid_primer(high_gc))
        length, gc, tm = self.check_primer_properties(good_gc)
        self.assertTrue(is_valid_primer(good_gc))

    def test_melting_temperature(self):
        # 测试熔解温度限制
        low_tm = "ATATATATATATATAT"  # Tm < 47°C
        high_tm = "GCGCGCGCGCGCGCGC"  # Tm > 53°C
        good_tm = self.generate_valid_primer(17)  # 适中的Tm
        
        self.assertFalse(is_valid_primer(low_tm))
        self.assertFalse(is_valid_primer(high_tm))
        length, gc, tm = self.check_primer_properties(good_tm)
        self.assertTrue(is_valid_primer(good_tm))

    def test_poly_n_repeats(self):
        # 测试连续碱基重复限制
        poly_a = "ATGCAAAAATGCTAGCT"  # 5个连续A
        poly_t = "ATGCTTTTTGCTAGCT"  # 5个连续T
        poly_g = "ATGCGGGGGCTAGCT"  # 5个连续G
        poly_c = "ATGCCCCCCTAGCT"  # 5个连续C
        no_poly = self.generate_valid_primer(17)  # 无连续重复
        
        self.assertFalse(is_valid_primer(poly_a))
        self.assertFalse(is_valid_primer(poly_t))
        self.assertFalse(is_valid_primer(poly_g))
        self.assertFalse(is_valid_primer(poly_c))
        length, gc, tm = self.check_primer_properties(no_poly)
        self.assertTrue(is_valid_primer(no_poly))

    def test_custom_parameters(self):
        # 测试自定义参数
        primer = "ATGCTAGCTAGCTGATCG"
        # 使用更严格的参数
        self.assertFalse(is_valid_primer(
            primer,
            min_length=20,  # 更长的最小长度
            gc_min=50.0,    # 更高的最小GC含量
            tm_min=60.0     # 更高的最小熔点
        ))
        
        # 使用更宽松的参数
        self.assertTrue(is_valid_primer(
            primer,
            min_length=15,  # 更短的最小长度
            gc_min=30.0,    # 更低的最小GC含量
            gc_max=70.0,    # 更高的最大GC含量
            tm_min=45.0,    # 更低的最小熔点
            tm_max=65.0     # 更高的最大熔点
        ))

if __name__ == '__main__':
    unittest.main()
