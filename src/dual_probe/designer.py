"""双探针设计核心逻辑模块"""
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
import logging
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from collections import Counter

from .config import DualProbeConfig
from ..common.sequence_utils import (
    has_hairpin,
    calculate_gc_content, 
    check_poly_n,
    has_dimer_issues,
    calculate_tm,
    is_valid_probe
)

logger = logging.getLogger(__name__)

def greedy(position_list):
    """Using greedy algorithm to select the maximum number of non-overlapping probes
    
    Args:
        position_list: a list of probe positions, like [[1, 12], [5, 17], [20, 32]]
        
    Returns:
        list: Selected non-overlapping probe positions
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



@dataclass
class ProbeSet:
    """探针组合结果类"""
    left_probe: Tuple[str, int, int, float, float]  # 序列,起始位置,结束位置,Tm值,GC含量
    right_probe: Tuple[str, int, int, float, float]

class DualProbeDesigner:
    """双探针设计器"""
    def __init__(self, config: DualProbeConfig):
        self.config = config

    def _generate_all_probes(self, seq: str, probe_size: int) -> Dict[int, str]:
        """生成所有可能的探针"""
        probes = {}
        for i in range(len(seq) - probe_size + 1):
            probes[i] = seq[i: i + probe_size]
        return probes

    def _filter_by_polyN(self, probes: Dict[int, str], polyN: int = 5) -> Dict[int, str]:
        """过滤含有多聚体N的探针"""
        filtered_probes = {}
        for pos, probe in probes.items():
            if check_poly_n(probe, polyN):
                filtered_probes[pos] = probe
                
        return filtered_probes

    def _filter_by_gc(self, probes: Dict[int, str]) -> Dict[int, str]:
        """过滤GC含量不合适的探针"""
        filtered_probes = {}
        
        for pos, probe in probes.items():
            left_seq = probe[:self.config.left_probe_length]
            right_seq = probe[-self.config.right_probe_length:]
            
            left_gc = calculate_gc_content(left_seq)
            right_gc = calculate_gc_content(right_seq)
            
            if (self.config.min_gc <= left_gc <= self.config.max_gc and
                self.config.min_gc <= right_gc <= self.config.max_gc):
                filtered_probes[pos] = probe
                
        return filtered_probes

    def _filter_by_tm(self, probes: Dict[int, str]) -> Dict[int, str]:
        """过滤Tm值不合适的探针"""
        filtered_probes = {}
        
        for pos, probe in probes.items():
            left_seq = probe[:self.config.left_probe_length]
            right_seq = probe[-self.config.right_probe_length:]
            
            left_tm = calculate_tm(left_seq)
            right_tm = calculate_tm(right_seq)
            
            if (self.config.min_tm <= left_tm <= self.config.max_tm and
                self.config.min_tm <= right_tm <= self.config.max_tm):
                filtered_probes[pos] = probe
                
        return filtered_probes

    def _filter_by_low_complexity(self, 
                                   probes: Dict[int, str], 
                                   seq: str) -> Dict[int, str]:
        """过滤具有重复k-mer的探针"""
        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return ''.join([complement[base] for base in reversed(seq)])
        
        k = self.config.kmer_size
        seq_rc = reverse_complement(seq)
        
        # 生成k-mer计数
        seq_kmers = Counter([seq[i:i+k] for i in range(len(seq) - k + 1)])
        rc_kmers = Counter([seq_rc[i:i+k] for i in range(len(seq_rc) - k + 1)])
        
        # 合并重复k-mer
        repeat_kmers = set()
        for kmer, count in seq_kmers.items():
            if count > 1:
                repeat_kmers.add(kmer)
        for kmer, count in rc_kmers.items():
            if count > 1:
                repeat_kmers.add(kmer)
        
        # 过滤探针
        filtered_probes = {}
        for pos, probe in probes.items():
            probe_kmers = {probe[i:i+k] for i in range(len(probe) - k + 1)}
            if not probe_kmers.intersection(repeat_kmers):
                filtered_probes[pos] = probe
                
        return filtered_probes

    def _select_non_overlapping_probes(self, 
                                       probes: Dict[int, str],
                                       method: str = 'dp') -> Dict[int, str]:
        """选择互不重叠的探针集合"""
        if not probes:
            return {}
            
        # 构建位置区间列表
        intervals = []
        for pos, probe in probes.items():
            end = pos + len(probe) + self.config.min_probe_spacing
            intervals.append([pos, end])
            
        # 按结束位置排序
        intervals.sort(key=lambda x: x[1])
        
        if method == 'greedy':
            selected = greedy(intervals)
        elif method == 'dp' :
            selected = dp(intervals)
        else:
            raise ValueError(f"不支持的方法: {method}")
        
        # 根据选中的区间过滤探针
        filtered_probes = {}
        selected_starts = {interval[0] for interval in selected}
        for pos, probe in probes.items():
            if pos in selected_starts:
                filtered_probes[pos] = probe
                
        return filtered_probes

    def design_probes(self, sequence: str) -> List[ProbeSet]:
        """设计双探针组合"""
        # 清理输入序列
        sequence = self._get_sequence(sequence)
        
        # 计算总探针长度
        total_length = (self.config.left_probe_length + 
                       self.config.gap_size + 
                       self.config.right_probe_length)
        
        # 生成并过滤探针
        probes = self._generate_all_probes(sequence, total_length)
        logger.info(f"初始探针数量: {len(probes)}")
        
        probes = self._filter_by_polyN(probes, self.config.max_poly_n)
        logger.info(f"多聚体N过滤后: {len(probes)}")
        
        probes = self._filter_by_gc(probes)
        logger.info(f"GC含量过滤后: {len(probes)}")
        
        probes = self._filter_by_tm(probes)
        logger.info(f"Tm值过滤后: {len(probes)}")
        
        probes = self._filter_by_kmer(probes, sequence)
        logger.info(f"k-mer过滤后: {len(probes)}")
        
        probes = self._select_non_overlapping_probes(probes)
        logger.info(f"重叠过滤后: {len(probes)}")
        
        # 转换为ProbeSet列表
        probe_sets = []
        for start_pos, probe_seq in probes.items():
            left_seq = probe_seq[:self.config.left_probe_length]
            right_seq = probe_seq[-self.config.right_probe_length:]
            
            probe_set = ProbeSet(
                left_probe=(
                    left_seq,
                    start_pos,
                    start_pos + self.config.left_probe_length,
                    calculate_tm(left_seq),
                    calculate_gc_content(left_seq)
                ),
                right_probe=(
                    right_seq,
                    start_pos + self.config.left_probe_length + self.config.gap_size,
                    start_pos + total_length,
                    calculate_tm(right_seq),
                    calculate_gc_content(right_seq)
                )
            )
            probe_sets.append(probe_set)
            
        logger.info(f"成功设计 {len(probe_sets)} 个双探针组合")
        return probe_sets


