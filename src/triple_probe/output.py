"""三探针输出处理模块"""
from typing import List, Tuple
from pathlib import Path

import csv
import logging
from Bio.SeqUtils import MeltingTemp as mt
from ..common.blast_utils import analyze_blast_results

from .designer import ProbeSet

logger = logging.getLogger(__name__)

class ProbeOutputHandler:
    """探针输出处理器"""
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def save_probe_sets(self, probe_sets: List[ProbeSet], 
                       task_name: str, output_prefix: str,
                       blast_db: str = None):
        """保存探针组合结果
        
        Args:
            probe_sets: 探针组合列表
            task_name: 任务名称
            output_prefix: 输出文件前缀
            blast_db: BLAST数据库路径（可选）
        """
        # 保存详细结果
        logger.info(f"保存详细结果")
        self._save_detailed_results(probe_sets, output_prefix, blast_db)
        
        # 保存BED格式
        self._save_bed_format(probe_sets, task_name, output_prefix)
        
        # 保存CSV格式
        self._save_csv_format(probe_sets, task_name, output_prefix)

    def _save_detailed_results(self, probe_sets: List[ProbeSet], 
                             output_prefix: str,
                             blast_db: str = None):
        """保存详细结果
        
        Args:
            probe_sets: 探针组合列表
            output_prefix: 输出文件前缀
            blast_db: BLAST数据库路径（可选）
        """
        output_file = self.output_dir / f"{output_prefix}_detailed.txt"
        
        with open(output_file, 'w') as f:
            f.write("三探针设计结果报告\n")
            f.write("=" * 80 + "\n\n")
            
            for i, probe_set in enumerate(probe_sets, 1):
                f.write(f"探针组合 {i}:\n")
                f.write("-" * 80 + "\n")
                
                # 写入每个探针的详细信息
                self._write_probe_details(f, "Left", probe_set.left_probe, blast_db)
                self._write_probe_details(f, "Middle", probe_set.middle_probe, blast_db)
                self._write_probe_details(f, "Right", probe_set.right_probe, blast_db)
                
                # 写入探针组合的评估信息
                self._write_probe_set_evaluation(f, probe_set)
                f.write("\n" + "=" * 80 + "\n\n")
    
    def _write_probe_details(self, file, probe_type: str, 
                            probe_info: Tuple,
                            blast_db: str = None):
        """写入单个探针详细信息
        
        Args:
            file: 输出文件对象
            probe_type: 探针类型（Left/Middle/Right）
            probe_info: 探针信息元组
            blast_db: BLAST数据库路径（可选）
        """
        sequence, start, end, tm, gc = probe_info
        
        file.write(f"\n{probe_type} Primer:\n")
        file.write(f"序列: {sequence}\n")
        file.write(f"位置: {start + 1}-{end}\n")  # 转换为1-based坐标
        file.write(f"长度: {len(sequence)}bp\n")
        file.write(f"GC含量: {gc:.1f}%\n")
        file.write(f"Tm值: {tm:.1f}°C\n")
        
        # 如果提供了BLAST数据库，添加BLAST分析结果
        if blast_db:
            blast_results = self._run_blast_analysis(sequence, blast_db)
            self._write_blast_results(file, blast_results)
            
        file.write("\n")
    
    def _write_probe_set_evaluation(self, file, probe_set: ProbeSet):
        """写入探针组合的评估信息
        
        Args:
            file: 输出文件对象
            probe_set: 探针组合
        """
        file.write("\n探针组合评估:\n")
        
        # 计算间隔
        left_to_middle = probe_set.middle_probe[1] - probe_set.left_probe[2]
        middle_to_right = probe_set.right_probe[1] - probe_set.middle_probe[2]
        
        file.write(f"Left-Middle间隔: {left_to_middle}bp\n")
        file.write(f"Middle-Right间隔: {middle_to_right}bp\n")
        
        # 计算总长度
        total_length = (probe_set.right_probe[2] - probe_set.left_probe[1])
        file.write(f"总跨度: {total_length}bp\n")
        
        # 计算平均Tm和GC
        avg_tm = (probe_set.left_probe[3] + 
                 probe_set.middle_probe[3] + 
                 probe_set.right_probe[3]) / 3
        avg_gc = (probe_set.left_probe[4] + 
                 probe_set.middle_probe[4] + 
                 probe_set.right_probe[4]) / 3
        
        file.write(f"平均Tm值: {avg_tm:.1f}°C\n")
        file.write(f"平均GC含量: {avg_gc:.1f}%\n")

    def _save_bed_format(self, probe_sets: List[ProbeSet], 
                        task_name: str, output_prefix: str):
        """保存BED格式结果
        
        Args:
            probe_sets: 探针组合列表
            task_name: 任务名称
            output_prefix: 输出文件前缀
        """
        bed_file = self.output_dir / f"{output_prefix}.bed"
        
        with open(bed_file, 'w') as f:
            for i, probe_set in enumerate(probe_sets, 1):
                # 写入左探针
                f.write(f"{task_name}\t{probe_set.left_probe[1]}\t"
                       f"{probe_set.left_probe[2]}\tL-{i}\t0\t+\n")
                
                # 写入中间探针
                f.write(f"{task_name}\t{probe_set.middle_probe[1]}\t"
                       f"{probe_set.middle_probe[2]}\tM-{i}\t0\t+\n")
                
                # 写入右探针
                f.write(f"{task_name}\t{probe_set.right_probe[1]}\t"
                       f"{probe_set.right_probe[2]}\tR-{i}\t0\t+\n")

    def _save_csv_format(self, probe_sets: List[ProbeSet], 
                        task_name: str, output_prefix: str):
        """保存CSV格式结果
        
        Args:
            probe_sets: 探针组合列表
            task_name: 任务名称
            output_prefix: 输出文件前缀
        """
        csv_file = self.output_dir / f"{output_prefix}.csv"
        
        headers = ['Set_ID', 'Primer_Type', 'Sequence', 'Start', 'End', 
                  'Length', 'Tm', 'GC_Content']
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            
            for i, probe_set in enumerate(probe_sets, 1):
                # 写入每个探针的信息
                for probe_type, probe_info in [
                    ('Left', probe_set.left_probe),
                    ('Middle', probe_set.middle_probe),
                    ('Right', probe_set.right_probe)
                ]:
                    sequence, start, end, tm, gc = probe_info
                    writer.writerow([
                        f"{task_name}_{i}",
                        probe_type,
                        sequence,
                        start + 1,  # 转换为1-based坐标
                        end,
                        len(sequence),
                        f"{tm:.1f}",
                        f"{gc:.1f}"
                    ])

    def _run_blast_analysis(self, sequence: str, blast_db: str) -> dict:
        """运行BLAST分析
        
        Args:
            sequence: 探针序列
            blast_db: BLAST数据库路径
            
        Returns:
            dict: BLAST分析结果
        """
        # 这里需要实现BLAST分析逻辑
        # 可以调用之前实现的blast相关函数
        try:
            # 获取错配统计和详细匹配信息
            mismatch_stats, detailed_matches = analyze_blast_results(sequence, blast_db)
            
            return {
                'mismatch_stats': mismatch_stats,
                'detailed_matches': detailed_matches
            }
            
        except Exception as e:
            logger.error(f"BLAST分析出错: {e}")
            return {
                'mismatch_stats': {},
                'detailed_matches': [],
                'alignment_details': []
            }

    def _write_blast_results(self, file, blast_results: dict):
        """写入BLAST分析结果
        
        Args:
            file: 输出文件对象
            blast_results: BLAST分析结果
        """
        mismatch_stats = blast_results['mismatch_stats']
        detailed_matches = blast_results['detailed_matches']
        
        file.write("\nBLAST分析结果:\n")
        
        # 写入错配统计
        file.write("\nBLAST Analysis:\n")
        file.write("Mismatches\tCount\tDescription\n")
        file.write("-" * 50 + "\n")
        for mismatches, count in sorted(mismatch_stats.items()):
            description = "Perfect match" if mismatches == 0 else f"{mismatches} mismatch(es)"
            file.write(f"{mismatches}\t\t{count}\t\t{description}\n")

        # 写入详细匹配信息
        file.write("\nDetailed Matches:\n")
        file.write("Subject ID\tMismatches\tGaps\tTotal Mismatches\n")
        file.write("-" * 70 + "\n")
        for subject_id, mismatches, gaps in detailed_matches:
            total = mismatches + gaps
            file.write(f"{subject_id}\t{mismatches}\t{gaps}\t{total}\n")