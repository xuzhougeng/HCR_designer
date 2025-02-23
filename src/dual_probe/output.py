"""双探针输出处理模块"""
from typing import List, Tuple
from pathlib import Path
import csv
import logging
from ..common.blast_utils import analyze_blast_results
import json
from ..common.probeset_select import select_probe_sets

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
        """保存探针组合结果，并根据质量得分排序
        
        Args:
            probe_sets: 探针组合列表
            task_name: 任务名称
            output_prefix: 输出文件前缀
            blast_db: BLAST数据库路径（可选）
            
        Returns:
            List[ProbeSet]: 按质量得分排序后的探针组合列表
        """
        # 如果有BLAST数据库，先计算所有探针的BLAST结果并排序
        if blast_db:
            # 存储所有探针的BLAST结果，避免重复计算
            blast_cache = {}
            scored_probes = []
            
            for probe_set in probe_sets:
                # 获取或计算BLAST结果
                left_seq = probe_set.left_probe[0]
                right_seq = probe_set.right_probe[0]
                
                if left_seq not in blast_cache:
                    blast_cache[left_seq] = self._run_blast_analysis(left_seq, blast_db)
                if right_seq not in blast_cache:
                    blast_cache[right_seq] = self._run_blast_analysis(right_seq, blast_db)
                
                # 计算探针组合得分
                blast_results_dict = {
                    'left': blast_cache[left_seq],
                    'right': blast_cache[right_seq]
                }
                score = self._calculate_probe_set_score(blast_results_dict)
                scored_probes.append((probe_set, score, blast_results_dict))
            
            # 按得分排序（得分越低越好）
            scored_probes.sort(key=lambda x: x[1])
            
            # 保存详细结果
            self._save_detailed_results_with_cache(scored_probes, output_prefix)
            
            # 更新排序后的probe_sets
            probe_sets = [p[0] for p in scored_probes]
        else:
            # 如果没有BLAST数据库，直接保存结果
            self._save_detailed_results(probe_sets, output_prefix, None)
        
        # 保存其他格式
        self._save_bed_format(probe_sets, task_name, output_prefix)
        self._save_csv_format(probe_sets, task_name, output_prefix)
        
        return probe_sets

    def _save_detailed_results(self, probe_sets: List[ProbeSet], 
                             output_prefix: str,
                             blast_db: str = None):
        """保存详细结果"""
        output_file = self.output_dir / f"{output_prefix}_detailed.txt"
        
        with open(output_file, 'w') as f:
            f.write("双探针设计结果报告\n")
            f.write("=" * 80 + "\n\n")
            
            for i, probe_set in enumerate(probe_sets, 1):
                f.write(f"探针组合 {i}:\n")
                f.write("-" * 80 + "\n")
                
                # 收集每个探针的BLAST结果
                blast_results_dict = {
                    'left': self._write_primer_details(f, "Left", probe_set.left_probe, blast_db),
                    'right': self._write_primer_details(f, "Right", probe_set.right_probe, blast_db)
                }
                
                # 计算探针组合得分
                if blast_db:
                    probe_set_score = self._calculate_probe_set_score(blast_results_dict)
                else:
                    probe_set_score = None
                
                # 写入探针组合的评估信息
                self._write_probe_set_evaluation(f, probe_set, probe_set_score)
                f.write("\n" + "=" * 80 + "\n\n")
    
    def _write_primer_details(self, file, primer_type: str, 
                            primer_info: Tuple,
                            blast_db: str = None) -> dict:
        """写入单个探针详细信息并返回BLAST结果
        
        Args:
            file: 输出文件对象
            primer_type: 探针类型（Left/Right）
            primer_info: 探针信息元组
            blast_db: BLAST数据库路径（可选）
            
        Returns:
            dict: BLAST分析结果
        """
        sequence, start, end, tm, gc = primer_info
        blast_results = {}
        
        file.write(f"\n{primer_type} Primer:\n")
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
        return blast_results
    
    def _calculate_probe_set_score(self, blast_results_dict: dict) -> float:
        """计算探针组合的总体得分
        
        评分规则:
        1. 每个探针的错配数量（越少越好）
        2. 探针的特异性（唯一匹配更好）
        3. 连续错配的惩罚
        
        Args:
            blast_results_dict: 包含每个探针BLAST结果的字典
            
        Returns:
            float: 探针组合的得分（分数越低越好）
        """
        total_score = 0
        
        for probe_type, blast_results in blast_results_dict.items():
            mismatch_stats = blast_results.get('mismatch_stats', {})
            detailed_matches = blast_results.get('detailed_matches', [])
            
            # 计算错配分数
            if mismatch_stats:
                min_mismatches = min(mismatch_stats.keys())
                mismatch_score = min_mismatches * 10  # 每个错配加10分
            else:
                mismatch_score = 50  # 如果没有匹配，给予较高的惩罚分数
                
            # 计算多重匹配惩罚
            multiple_match_penalty = max(0, len(detailed_matches) - 1) * 5
            
            # 连续错配惩罚
            consecutive_penalty = 0
            for match in detailed_matches:
                if "consecutive_mismatches" in match:
                    consecutive_penalty += 20
                    
            probe_score = mismatch_score + multiple_match_penalty + consecutive_penalty
            total_score += probe_score
        
        return total_score
    
    def _write_probe_set_evaluation(self, file, probe_set: ProbeSet, probe_set_score: float = None):
        """写入探针组合的评估信息"""
        file.write("\n探针组合评估:\n")
        
        # 计算间隔
        gap = probe_set.right_probe[1] - probe_set.left_probe[2]
        file.write(f"探针间隔: {gap}bp\n")
        
        # 计算总长度
        total_length = probe_set.right_probe[2] - probe_set.left_probe[1]
        file.write(f"总跨度: {total_length}bp\n")
        
        # 计算平均Tm和GC
        avg_tm = (probe_set.left_probe[3] + probe_set.right_probe[3]) / 2
        avg_gc = (probe_set.left_probe[4] + probe_set.right_probe[4]) / 2
        
        file.write(f"平均Tm值: {avg_tm:.1f}°C\n")
        file.write(f"平均GC含量: {avg_gc:.1f}%\n")
        
        # 如果有探针组合得分，写入得分信息
        if probe_set_score is not None:
            file.write(f"探针组合质量得分: {probe_set_score:.1f} (分数越低越好)\n")

    def _save_bed_format(self, probe_sets: List[ProbeSet], 
                        task_name: str, output_prefix: str):
        """保存BED格式结果"""
        bed_file = self.output_dir / f"{output_prefix}.bed"
        
        with open(bed_file, 'w') as f:
            for i, probe_set in enumerate(probe_sets, 1):
                # 写入左探针
                f.write(f"{task_name}\t{probe_set.left_probe[1]}\t"
                       f"{probe_set.left_probe[2]}\tL-{i}\t0\t+\n")
                
                # 写入右探针
                f.write(f"{task_name}\t{probe_set.right_probe[1]}\t"
                       f"{probe_set.right_probe[2]}\tR-{i}\t0\t+\n")

    def _save_csv_format(self, probe_sets: List[ProbeSet], 
                        task_name: str, output_prefix: str):
        """保存CSV格式结果"""
        csv_file = self.output_dir / f"{output_prefix}.csv"
        
        headers = ['Set_ID', 'Primer_Type', 'Sequence', 'Start', 'End', 
                  'Length', 'Tm', 'GC_Content']
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            
            for i, probe_set in enumerate(probe_sets, 1):
                # 写入每个探针的信息
                for primer_type, primer_info in [
                    ('Left', probe_set.left_probe),
                    ('Right', probe_set.right_probe)
                ]:
                    sequence, start, end, tm, gc = primer_info
                    writer.writerow([
                        f"{task_name}_{i}",
                        primer_type,
                        sequence,
                        start + 1,  # 转换为1-based坐标
                        end,
                        len(sequence),
                        f"{tm:.1f}",
                        f"{gc:.1f}"
                    ])

    def _run_blast_analysis(self, sequence: str, blast_db: str) -> dict:
        """运行BLAST分析"""
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

    def _save_detailed_results_with_cache(self, scored_probes: List[Tuple], output_prefix: str):
        """使用缓存的BLAST结果保存详细信息
        
        Args:
            scored_probes: 包含(probe_set, score, blast_results_dict)的列表
            output_prefix: 输出文件前缀
        """
        # 保存文本格式
        output_file = self.output_dir / f"{output_prefix}_detailed.txt"
        
        # 创建JSON数据结构
        json_data = {
            "report_title": "双探针设计结果报告",
            "probe_sets": []
        }
        
        with open(output_file, 'w') as f:
            f.write("双探针设计结果报告\n")
            f.write("=" * 80 + "\n\n")
            
            for i, (probe_set, score, blast_results_dict) in enumerate(scored_probes, 1):
                # 文本格式输出
                f.write(f"探针组合 {i}:\n")
                f.write("-" * 80 + "\n")
                
                # 写入每个探针的详细信息
                self._write_primer_details_with_cache(f, "Left", probe_set.left_probe, 
                                                    blast_results_dict['left'])
                self._write_primer_details_with_cache(f, "Right", probe_set.right_probe, 
                                                    blast_results_dict['right'])
                
                # 写入探针组合的评估信息
                self._write_probe_set_evaluation(f, probe_set, score)
                f.write("\n" + "=" * 80 + "\n\n")
                
                # 收集JSON数据
                probe_set_data = {
                    "id": i,
                    "left_probe": self._get_primer_json_data(probe_set.left_probe, blast_results_dict['left']),
                    "right_probe": self._get_primer_json_data(probe_set.right_probe, blast_results_dict['right']),
                    "evaluation": self._get_evaluation_json_data(probe_set, score)
                }
                json_data["probe_sets"].append(probe_set_data)
        
        # 保存JSON格式
        json_output_file = self.output_dir / f"{output_prefix}_detailed.json"
        with open(json_output_file, 'w') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
    
    def _get_primer_json_data(self, primer_info: Tuple, blast_results: dict) -> dict:
        """生成探针的JSON数据
        
        Args:
            primer_info: 探针信息元组
            blast_results: BLAST分析结果
            
        Returns:
            dict: 探针的JSON格式数据
        """
        sequence, start, end, tm, gc = primer_info
        data = {
            "sequence": sequence,
            "position": {
                "start": start + 1,  # 转换为1-based坐标
                "end": end
            },
            "length": len(sequence),
            "gc_content": f"{gc:.1f}",
            "tm": f"{tm:.1f}"
        }
        
        if blast_results:
            data["blast_results"] = {
                "mismatch_stats": blast_results['mismatch_stats'],
                "detailed_matches": [
                    {
                        "subject_id": subject_id,
                        "mismatches": mismatches,
                        "gaps": gaps,
                        "total_mismatches": mismatches + gaps
                    }
                    for subject_id, mismatches, gaps in blast_results['detailed_matches']
                ]
            }
        
        return data
    
    def _get_evaluation_json_data(self, probe_set: ProbeSet, probe_set_score: float = None) -> dict:
        """生成评估信息的JSON数据
        
        Args:
            probe_set: 探针组合
            probe_set_score: 探针组合的质量得分（可选）
            
        Returns:
            dict: 评估信息的JSON格式数据
        """
        # 计算间隔
        gap = probe_set.right_probe[1] - probe_set.left_probe[2]
        
        # 计算总长度
        total_length = probe_set.right_probe[2] - probe_set.left_probe[1]
        
        # 计算平均Tm和GC
        avg_tm = (probe_set.left_probe[3] + probe_set.right_probe[3]) / 2
        avg_gc = (probe_set.left_probe[4] + probe_set.right_probe[4]) / 2
        
        data = {
            "gap": gap,
            "total_span": total_length,
            "average_tm": f"{avg_tm:.1f}",
            "average_gc_content": f"{avg_gc:.1f}"
        }
        
        if probe_set_score is not None:
            data["quality_score"] = f"{probe_set_score:.1f}"
        
        return data

    def _write_primer_details_with_cache(self, file, primer_type: str, 
                                       primer_info: Tuple,
                                       blast_results: dict):
        """使用缓存的BLAST结果写入探针详细信息
        
        Args:
            file: 输出文件对象
            primer_type: 探针类型
            primer_info: 探针信息元组
            blast_results: 缓存的BLAST结果
        """
        sequence, start, end, tm, gc = primer_info
        
        file.write(f"\n{primer_type} Primer:\n")
        file.write(f"序列: {sequence}\n")
        file.write(f"位置: {start + 1}-{end}\n")
        file.write(f"长度: {len(sequence)}bp\n")
        file.write(f"GC含量: {gc:.1f}%\n")
        file.write(f"Tm值: {tm:.1f}°C\n")
        
        if blast_results:
            self._write_blast_results(file, blast_results)
            
        file.write("\n")

    def save_filtered_probe_sets(self, json_data: dict, output_prefix: str):
        """
        基于select_probe_sets对探针组合进行筛选并输出报告
        
        Args:
            json_data: 包含探针组合数据的字典
            output_prefix: 输出文件前缀
        """
        # 使用select_probe_sets进行筛选
        filtered_data = select_probe_sets(json_data)
        
        # 保存筛选后的完整JSON数据
        json_output_file = self.output_dir / f"{output_prefix}_filtered.json"
        with open(json_output_file, 'w', encoding='utf-8') as f:
            json.dump(filtered_data, f, indent=2, ensure_ascii=False)
        
        # 生成筛选报告
        report_file = self.output_dir / f"{output_prefix}_filtered_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("探针组合筛选报告\n")
            f.write("=" * 80 + "\n\n")
            
            # 统计信息
            total_sets = len(filtered_data['probe_sets'])
            selected_sets = sum(1 for probe_set in filtered_data['probe_sets'] if probe_set.get('keep', False))
            
            f.write(f"总探针组合数: {total_sets}\n")
            f.write(f"筛选后保留数: {selected_sets}\n")
            f.write("\n" + "-" * 80 + "\n\n")
            
            # 输出保留的探针组合详细信息
            f.write("保留的探针组合:\n")
            f.write("-" * 80 + "\n")
            for probe_set in filtered_data['probe_sets']:
                if probe_set.get('keep', False):
                    f.write(f"\n探针组合 ID: {probe_set['id']}\n")
                    
                    # 左探针信息
                    left_probe = probe_set['left_probe']
                    f.write("\n左探针:\n")
                    f.write(f"序列: {left_probe['sequence']}\n")
                    f.write(f"位置: {left_probe['position']['start']}-{left_probe['position']['end']}\n")
                    f.write(f"GC含量: {left_probe['gc_content']}%\n")
                    f.write(f"Tm值: {left_probe['tm']}°C\n")
                    
                    # 右探针信息
                    right_probe = probe_set['right_probe']
                    f.write("\n右探针:\n")
                    f.write(f"序列: {right_probe['sequence']}\n")
                    f.write(f"位置: {right_probe['position']['start']}-{right_probe['position']['end']}\n")
                    f.write(f"GC含量: {right_probe['gc_content']}%\n")
                    f.write(f"Tm值: {right_probe['tm']}°C\n")
                    
                    # 评估信息
                    eval_info = probe_set['evaluation']
                    f.write("\n评估信息:\n")
                    f.write(f"质量得分: {eval_info['quality_score']}\n")
                    f.write(f"探针间隔: {eval_info['gap']}bp\n")
                    f.write(f"总跨度: {eval_info['total_span']}bp\n")
                    f.write(f"平均Tm值: {eval_info['average_tm']}°C\n")
                    f.write(f"平均GC含量: {eval_info['average_gc_content']}%\n")
                    f.write(f"筛选原因: {probe_set['reason']}\n")
                    f.write("\n" + "-" * 80 + "\n")
            
            # 输出未保留的探针组合简要信息
            f.write("\n\n未保留的探针组合:\n")
            f.write("-" * 80 + "\n")
            for probe_set in filtered_data['probe_sets']:
                if not probe_set.get('keep', False):
                    f.write(f"\n探针组合 ID: {probe_set['id']}\n")
                    f.write(f"质量得分: {probe_set['evaluation']['quality_score']}\n")
                    f.write(f"未保留原因: {probe_set['reason']}\n")
                    f.write("-" * 40 + "\n")
        
        # 保存筛选后的BED格式
        bed_file = self.output_dir / f"{output_prefix}_filtered.bed"
        task_name = output_prefix.split('_filtered')[0]  # 获取原始task_name
        with open(bed_file, 'w') as f:
            for probe_set in filtered_data['probe_sets']:
                if probe_set.get('keep', False):
                    set_id = probe_set['id']
                    # 写入左探针
                    f.write(f"{task_name}\t{probe_set['left_probe']['position']['start']}\t"
                           f"{probe_set['left_probe']['position']['end']}\tL-{set_id}\t0\t+\n")
                    
                    # 写入右探针
                    f.write(f"{task_name}\t{probe_set['right_probe']['position']['start']}\t"
                           f"{probe_set['right_probe']['position']['end']}\tR-{set_id}\t0\t+\n")
        
        # 保存筛选后的CSV格式
        csv_file = self.output_dir / f"{output_prefix}_filtered.csv"
        headers = ['Set_ID', 'Primer_Type', 'Sequence', 'Start', 'End', 
                  'Length', 'Tm', 'GC_Content']
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            
            for probe_set in filtered_data['probe_sets']:
                if probe_set.get('keep', False):
                    set_id = probe_set['id']
                    # 写入每个探针的信息
                    for probe_type, probe_info in [
                        ('Left', probe_set['left_probe']),
                        ('Right', probe_set['right_probe'])
                    ]:
                        writer.writerow([
                            f"{task_name}_{set_id}",
                            probe_type,
                            probe_info['sequence'],
                            probe_info['position']['start'],
                            probe_info['position']['end'],
                            len(probe_info['sequence']),
                            probe_info['tm'],
                            probe_info['gc_content']
                        ])