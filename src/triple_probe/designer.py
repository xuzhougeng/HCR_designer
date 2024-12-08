"""三探针设计核心逻辑模块"""
from dataclasses import dataclass
from typing import List, Tuple, Optional
import logging
from Bio.Seq import Seq

from .config import TripleProbeConfig
from ..common.sequence_utils import (
    has_hairpin,
    calculate_gc_content,
    check_poly_n,
    has_dimer_issues,
    calculate_tm,
    is_valid_probe
)

logger = logging.getLogger(__name__)

@dataclass
class ProbeSet:
    """探针组合结果类"""
    left_probe: Tuple[str, int, int, float, float]  # 序列,起始位置,结束位置,Tm值,GC含量
    middle_probe: Tuple[str, int, int, float, float]
    right_probe: Tuple[str, int, int, float, float]

class TripleProbeDesigner:
    """三探针设计器"""
    def __init__(self, config: TripleProbeConfig):
        self.config = config
        
    def design_probes(self, sequence: str) -> List[ProbeSet]:
        """设计三探针组合
        
        Args:
            sequence: 目标序列
            
        Returns:
            List[ProbeSet]: 探针组合列表
        """
        sequence_length = len(sequence)
        probe_sets = []
        used_positions = set()
        
        pos = 0
        while pos < sequence_length - 3 * self.config.min_length:
            logger.debug(f"Searching for probe set starting at position {pos}")
            
            # 寻找左引物
            left_result = self._find_next_valid_probe(
                sequence, pos, sequence_length, used_positions
            )
            if not left_result:
                pos += 1
                continue
                
            # 寻找中间引物（必须紧接着左引物）
            middle_result = self._find_next_valid_probe(
                sequence, left_result[2], sequence_length, used_positions,
                must_start_at=left_result[2]  # 确保从左引物结束位置开始
            )
            if not middle_result:
                pos += 1
                continue
                
            # 检查左引物和中间引物的兼容性
            if not self._check_primer_compatibility(left_result, middle_result):
                pos += 1
                continue
                
            # 寻找右引物（与中间引物之间保持固定的gap）
            right_start = middle_result[2] + self.config.r_m_gap
            if right_start >= sequence_length - self.config.min_length:
                pos += 1
                continue
                
            right_result = self._find_next_valid_probe(
                sequence, right_start, sequence_length, used_positions
            )
            if not right_result:
                pos += 1
                continue
                
            # 检查三个引物的兼容性
            if self._check_primer_compatibility(left_result, middle_result, right_result):
                logger.info(f"Found valid probe set:")
                logger.info(f"Left probe: {left_result[0]}, Tm: {left_result[3]:.1f}")
                logger.info(f"Middle probe: {middle_result[0]}, Tm: {middle_result[3]:.1f}")
                logger.info(f"Right probe: {right_result[0]}, Tm: {right_result[3]:.1f}")
                
                probe_set = ProbeSet(left_result, middle_result, right_result)
                probe_sets.append(probe_set)
                
                # 更新已使用的位置
                for result in [left_result, middle_result, right_result]:
                    for i in range(result[1], result[2]):
                        used_positions.add(i)
                
                # 从左引物后继续搜索
                pos = left_result[2]
            else:
                pos += 1
                
        return probe_sets
    
    def _find_next_valid_probe(self, sequence: str, start_pos: int, 
                              end_pos: int, used_positions: set,
                              must_start_at: int = None) -> Optional[Tuple]:
        """在指定范围内找到第一个有效的引物
        
        Parameters:
        -----------
        start_pos : int
            起始位置
        end_pos : int
            结束位置
        used_pos : set
            已使用的位置集合
        must_start_at : int, optional
            必须从这个位置开始（用于确保Middle紧接Left）
        """
        # 如果指定了必须的起始位置，就只检查从这个位置开始的引物
        if must_start_at is not None:
            start_pos = must_start_at
            
        for pos in range(start_pos, end_pos - self.config.min_length + 1):
            # 如果指定了必须的起始位置，但当前位置不是该位置，则跳过
            if must_start_at is not None and pos != must_start_at:
                break
                
            # 检查位置是否已被使用
            if any(p in used_pos for p in range(pos, pos + self.config.max_length)):
                continue
                
            # 尝试不同长度的引物
            for length in range(self.config.min_length, 
                                min(self.config.max_length + 1, end_pos - pos + 1)):
                probe = sequence[pos:pos + length]
                try:
                    # 检查引物是否符合基本要求
                    if not is_valid_probe(probe, 
                                        min_length=self.config.min_length,
                                        max_length=self.config.max_length,
                                        gc_min=self.config.gc_min,
                                        gc_max=self.config.gc_max,
                                        tm_min=self.config.tm_min,
                                        tm_max=self.config.tm_max,
                                        poly_n=self.config.poly_n):
                        continue
                    # 在返回之前进行反向互补
                    rev_comp_probe = str(Seq(probe).reverse_complement())
                    logger.debug(f"Found valid probe: {rev_comp_probe}, Tm: {tm:.1f}, GC: {gc:.1f}")
                    return (rev_comp_probe, pos, pos + length, tm, gc)
                except Exception as e:
                    logger.error(f"计算引物参数时出错: {e}")
                    continue
            return None

    def _check_probe_compatibility(self, probe1, probe2, probe3=None):
        """检查引物之间是否存在互补性问题"""
        probes = [p[0] for p in [probe1, probe2, probe3] if p is not None]
        return not has_dimer_issues(probes, self.config.min_complementary_length)
