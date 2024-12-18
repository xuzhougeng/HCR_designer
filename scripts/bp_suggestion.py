from typing import List, Set, Dict, Tuple

def reverse_complement(sequence: str) -> str:
    """获取序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

class BridgeProbeSystem:
    def __init__(self, file_path: str, k: int = 3):
        """
        初始化系统
        :param file_path: probe_table.txt的路径
        :param k: kmer的长度
        """
        self.k = k
        self.probes = self._load_probes(file_path)
        self.kmer_sets = {}
        self._calculate_all_kmer_sets()
    
    def _load_probes(self, file_path: str) -> Dict[str, str]:
        """从文本文件加载探针数据"""
        probes = {}
        with open(file_path, 'r') as f:
            for line in f:
                bp_id, sequence = line.strip().split()
                probes[bp_id] = sequence
        return probes
    
    def _get_kmers(self, sequence: str) -> Set[str]:
        """获取序列的所有k-mers集合，包括反向互补序列的k-mers"""
        # 获取原序列的k-mers
        kmers = set(sequence[i:i+self.k] for i in range(len(sequence)-self.k+1))
        # 获取反向互补序列的k-mers
        rc_sequence = reverse_complement(sequence)
        rc_kmers = set(rc_sequence[i:i+self.k] for i in range(len(rc_sequence)-self.k+1))
        # 合并两个集合
        return kmers.union(rc_kmers)
    
    def _calculate_all_kmer_sets(self):
        """计算所有序列的kmer集合"""
        for bp_id, sequence in self.probes.items():
            self.kmer_sets[bp_id] = self._get_kmers(sequence)
    
    def find_dissimilar_probes(self, input_bp_ids: List[str], n: int) -> List[Tuple[str, str]]:
        """
        找到与输入BP_ID没有kmer重叠的探针
        :param input_bp_ids: 输入的BP_ID列表
        :param n: 需要返回的BP_ID数量
        :return: 不相似的BP_ID和序列的元组列表
        """
        # 获取输入序列的所有kmer
        input_kmers = set()
        for bp_id in input_bp_ids:
            if bp_id not in self.kmer_sets:
                raise ValueError(f"BP_ID {bp_id} not found in database")
            input_kmers.update(self.kmer_sets[bp_id])
        
        # 找到没有重叠的探针
        dissimilar_probes = []
        for bp_id, kmers in self.kmer_sets.items():
            if bp_id not in input_bp_ids:
                # 如果没有重叠，intersection为空集
                if not input_kmers.intersection(kmers):
                    dissimilar_probes.append((bp_id, self.probes[bp_id]))
                    if len(dissimilar_probes) == n:
                        break
        
        return dissimilar_probes
    
    def find_conflicting_probes(self, bp_id: str) -> List[Tuple[str, str, int, List[str]]]:
        """
        查找与指定BP_ID有k-mer重叠的其他探针
        
        Args:
            bp_id: 要查询的BP_ID
            
        Returns:
            有k-mer重叠的BP_ID和序列的元组列表
        """
        if bp_id not in self.kmer_sets:
            raise ValueError(f"BP_ID {bp_id} not found in database")
        
        # 获取输入序列的kmers
        input_kmers = self.kmer_sets[bp_id]
        
        # 查找有重叠的探针
        conflicting_probes = []
        for other_bp_id, kmers in self.kmer_sets.items():
            if other_bp_id != bp_id:
                # 如果有重叠，intersection不为空集
                overlap = input_kmers.intersection(kmers)
                if overlap:
                    # 检查重叠的k-mers中是否包含反向互补序列
                    rc_overlap = set()
                    for kmer in overlap:
                        rc_kmer = reverse_complement(kmer)
                        if rc_kmer in overlap:
                            rc_overlap.add(kmer)
                            rc_overlap.add(rc_kmer)
                    
                    conflicting_probes.append((
                        other_bp_id, 
                        self.probes[other_bp_id],
                        len(overlap),  # 添加重叠的k-mer数量
                        list(rc_overlap)[:3] if rc_overlap else list(overlap)[:3]  # 优先显示反向互补的k-mer对
                    ))
        
        # 按重叠k-mer数量降序排序
        conflicting_probes.sort(key=lambda x: x[2], reverse=True)
        return conflicting_probes
    
    def analyze_sequences_conflicts(self, sequences: List[str]) -> List[Tuple[str, List[Tuple[str, str, int, List[str]]]]]:
        """
        分析输入序列与已有探针的冲突情况
        
        Args:
            sequences: 要分析的序列列表
        
        Returns:
            每个输入序列的冲突情况列表，每项包含：
            - 输入序列
            - 与该序列有冲突的探针列表(BP_ID, 序列, 重叠数量, 重叠k-mer示例)
        """
        results = []
        
        for seq in sequences:
            # 获取输入序列的kmers
            seq_kmers = self._get_kmers(seq)
            
            # 查找有重叠的探针
            conflicts = []
            for bp_id, kmers in self.kmer_sets.items():
                overlap = seq_kmers.intersection(kmers)
                if overlap:
                    # 检查重叠的k-mers中是否包含反向互补序列
                    rc_overlap = set()
                    for kmer in overlap:
                        rc_kmer = reverse_complement(kmer)
                        if rc_kmer in overlap:
                            rc_overlap.add(kmer)
                            rc_overlap.add(rc_kmer)
                    
                    conflicts.append((
                        bp_id,
                        self.probes[bp_id],
                        len(overlap),
                        list(rc_overlap)[:3] if rc_overlap else list(overlap)[:3]  # 优先显示反向互补的k-mer对
                    ))
            
            # 按重叠k-mer数量降序排序
            conflicts.sort(key=lambda x: x[2], reverse=True)
            results.append((seq, conflicts))
        
        return results


def design_multiple_probes(probe_table_file: str, 
                           input_bp_ids: List[str]=[], 
                           n: int=1,
                           k: int=9) -> List[Tuple[str, str]]:
    """
    基于已有探针设计N个新的探针
    
    Args:
        input_bp_ids: 已有的BP_ID列表
        n: 需要设计的新探针数量
    
    Returns:
        新设计的BP_ID和序列的元组列表
    """
    # 用于存储结果
    result_probes = []
    current_bp_ids = input_bp_ids.copy()
    
    system = BridgeProbeSystem(probe_table_file, k=k)

    # 迭代设计n个探针
    while len(result_probes) < n:
        # 每次找一个不相似的探针
        dissimilar = system.find_dissimilar_probes(current_bp_ids, 1)
        
        # 如果找不到更多不相似的探针,退出循环
        if not dissimilar:
            break
            
        # 将找到的探针添加到结果中
        new_probe = dissimilar[0]
        result_probes.append(new_probe)
        
        # 将新探针的BP_ID添加到current_bp_ids中,用下一轮迭代
        current_bp_ids.append(new_probe[0])
        
    return result_probes


def query_conflicting_probes(probe_table_file: str, bp_id: str, k: int = 9) -> List[Tuple[str, str, int, List[str]]]:
    """
    查询与指定BP_ID有冲突的探针
    
    Args:
        probe_table_file: probe_table.txt的路径
        bp_id: 要查询的BP_ID
        k: kmer的长度
    
    Returns:
        有冲突的BP_ID、序列、重叠k-mer数量和示例k-mer的元组列表
    """
    system = BridgeProbeSystem(probe_table_file, k=k)
    return system.find_conflicting_probes(bp_id)


# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser(description="Bridge probe suggestion system")
#     # resources/probe_table.txt
#     parser.add_argument("-p", "--probe_table", type=str, required=True, help="Probe table file")
#     parser.add_argument("-k", "--kmer", type=int, required=True, help="Kmer length")
#     parser.add_argument("-i", "--input", type=str, required=True, help="Input BP_ID list")
#     parser.add_argument("-n", "--number", type=int, required=True, help="Number of dissimilar probes to return")
#     args = parser.parse_args()
    
#     # 使用示例
#     system = BridgeProbeSystem(args.probe_table, k=args.kmer)
#     input_bp_ids = args.input.split(",")
#     n = args.number
    
#     results = system.find_dissimilar_probes(input_bp_ids, n)
    
#     print(f"Found {len(results)} dissimilar probes:")
#     for bp_id, sequence in results:
#         print(f"BP_ID: {bp_id}, Sequence: {sequence}")

if __name__ == "__main__":
    import sys
    probe_table_file = sys.argv[1]
    result_probes = design_multiple_probes(probe_table_file,n=50,k=9)
    for bp_id, sequence in result_probes:
        print(f"BP_ID: {bp_id}, Sequence: {sequence}")
