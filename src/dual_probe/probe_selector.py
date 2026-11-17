"""ProbeSelector：按 TPM 从双探针候选中挑选探针对，修正 L 探针 5' 端碱基，输出 Excel。

移植自 ProbeSelector_Share/probe_filter.py。get_target_pairs() 至 select_probes() 的算法与原脚本
逐行一致（仅修正 dual_probe.csv 中任务名含下划线时 Set_ID 解析失败）；原 main() 拆成可传路径、
可返回内存流的函数，基因名与文件夹名的匹配规则不变。

输入（每个基因一个文件夹，文件夹名 = 基因表 Name 列）：
    dual_probe.bed            探针位置（制表符分隔，6 列以上）
    dual_probe_dual.csv       探针序列（列：Probe ID, Sequence）
    dual_probe.csv            可选，含 Set_ID / Primer_Type / Sequence / Tm / Length
    dual_probe_detailed.txt   可选，探针打分
基因信息表 CSV 必需列：Name, GeneID, BP_ID, TPM。
"""
import os
import csv
import io
import re
from collections import Counter
from typing import Dict, List, Tuple, Any, Optional, Union

def get_target_pairs(tpm):
    if tpm >= 100:
        return 4
    elif tpm >= 50:
        return 5
    elif tpm >= 30:
        return 6
    elif tpm >= 10:
        return 7
    elif tpm >= 5:
        return 8
    elif tpm >= 1:
        return 9
    else:
        return 10

def get_nucleotide_priority(l_seq, r_seq):
    if len(l_seq) == 0:
        return 8
    
    l_first = l_seq[0]
    
    if l_first in ['A', 'T']:
        return 5
    
    if l_first in ['C', 'G']:
        if len(l_seq) >= 2 and l_seq[1] in ['A', 'T']:
            return 6
        
        gc_count = 1
        for i in range(1, len(l_seq)):
            if l_seq[i] in ['C', 'G']:
                gc_count += 1
            else:
                break
        
        r_last = r_seq[-1] if len(r_seq) > 0 else ''
        
        if gc_count >= 2 and r_last in ['A', 'T']:
            return 6
        elif gc_count == 2 and r_last in ['C', 'G']:
            return 7
        elif gc_count >= 3 and r_last in ['C', 'G']:
            return 8
    
    return 8

def parse_dual_probe_bed(filepath: str) -> List[Dict[str, Union[int, str]]]:
    probes: Dict[int, Dict[str, Optional[Dict[str, int]]]] = {}
    gene_name = ''
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                gene_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                probe_name = parts[3]
                strand = parts[5]
                
                probe_type = probe_name.split('-')[0]
                probe_num = int(probe_name.split('-')[1])
                
                if probe_num not in probes:
                    probes[probe_num] = {'L': None, 'R': None}
                
                if probe_type == 'L':
                    probes[probe_num]['L'] = {'start': start, 'end': end}
                elif probe_type == 'R':
                    probes[probe_num]['R'] = {'start': start, 'end': end}
    
    probe_pairs: List[Dict[str, Union[int, str]]] = []
    for num in sorted(probes.keys()):
        pair = probes[num]
        if pair['L'] and pair['R']:
            probe_pairs.append({
                'num': num,
                'gene': gene_name,
                'L_start': pair['L']['start'],
                'L_end': pair['L']['end'],
                'R_start': pair['R']['start'],
                'R_end': pair['R']['end']
            })
    
    return probe_pairs

def parse_dual_probe_dual_csv(filepath: str) -> Dict[str, str]:
    probes: Dict[str, str] = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            probe_id = row['Probe ID']
            sequence = row['Sequence']
            probes[probe_id] = sequence
    return probes

def parse_dual_probe_csv(filepath: str) -> Dict[int, Dict[str, Optional[Dict[str, Union[str, float, int]]]]]:
    """解析dual_probe.csv，提取原始靶标序列和TM值"""
    probe_info: Dict[int, Dict[str, Optional[Dict[str, Union[str, float, int]]]]] = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            set_id = row['Set_ID']
            primer_type = row['Primer_Type']
            sequence = row['Sequence']
            tm = float(row['Tm']) if row['Tm'] else None
            length = int(row['Length']) if row['Length'] else len(sequence)
            
            # 任务名可含下划线（如 ABCG13_demo），编号取最后一段
            probe_num = int(set_id.rsplit('_', 1)[1])
            
            if probe_num not in probe_info:
                probe_info[probe_num] = {'L': None, 'R': None}
            
            if primer_type == 'Left':
                probe_info[probe_num]['L'] = {'sequence': sequence, 'tm': tm, 'length': length}
            elif primer_type == 'Right':
                probe_info[probe_num]['R'] = {'sequence': sequence, 'tm': tm, 'length': length}
    
    return probe_info

def calculate_tm(sequence: str) -> float:
    """计算DNA序列的熔解温度（Tm值），使用Nearest-Neighbor近似方法"""
    sequence = sequence.upper()
    length = len(sequence)
    
    if length == 0:
        return 0.0
    
    # 统计GC含量
    gc_count = sequence.count('G') + sequence.count('C')
    
    # 使用简化的Tm计算公式（适用于短序列）
    if length >= 14:
        tm = 64.9 + 41 * (gc_count - 16.4) / length
    else:
        # 短序列使用不同的公式
        tm = 2 * (sequence.count('A') + sequence.count('T')) + 4 * (sequence.count('G') + sequence.count('C'))
    
    return tm

def parse_dual_probe_detailed(filepath: str) -> Dict[int, float]:
    probe_scores: Dict[int, float] = {}
    current_probe: Optional[int] = None
    score = 0.0
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            match = re.search(r'探针组合 (\d+):', line)
            if match:
                if current_probe is not None:
                    probe_scores[current_probe] = score
                current_probe = int(match.group(1))
                score = 0.0
                continue
            
            match = re.search(r'质量得分: ([\d.]+)', line)
            if match:
                score = float(match.group(1))
        
        if current_probe is not None:
            probe_scores[current_probe] = score
    
    return probe_scores

def check_gap_requirement(candidate: Dict[str, Union[int, str]], selected: List[Dict[str, Union[int, str]]]) -> bool:
    for sp in selected:
        gap_before = candidate['L_start'] - sp['R_end']
        gap_after = sp['L_start'] - candidate['R_end']
        if 0 <= gap_before < 18 or 0 <= gap_after < 18:
            return False
    return True

def correct_probe_sequence(l_seq: str, r_seq: str, l_tm: Optional[float] = None, r_tm: Optional[float] = None) -> Tuple[str, str, str]:
    """
    修正探针序列，使L探针5'端以A或T开头，提高探针连接效率
    
    规则：
    1. L探针5'端单个C/G开头：将C/G转移至R探针3'端，或从R探针3'端转移A/T至L探针5'端
    2. L探针5'端两个C/G开头：将两个C/G转移至R探针3'端（需第3个碱基为A/T），或从R探针3'端转移A/T
    3. L探针5'端3个或更多C/G：放弃调整
    4. 调整后L/R探针TM值差异不能超过3℃
    
    参数：
    - l_seq, r_seq: L/R探针序列
    - l_tm, r_tm: 原始TM值（可选，如果提供则用于验证调整后的TM差异）
    
    返回：(修正后的L序列, 修正后的R序列, 修正说明)
    """
    if len(l_seq) == 0 or len(r_seq) == 0:
        return l_seq, r_seq, ''
    
    # 计算L探针5'端连续的C/G数量
    gc_count = 0
    for i in range(len(l_seq)):
        if l_seq[i] in ['C', 'G']:
            gc_count += 1
        else:
            break
    
    if gc_count == 0:
        # L探针已以A/T开头，无需调整
        return l_seq, r_seq, ''
    
    if gc_count >= 3:
        # L探针5'端有3个或更多C/G
        # 如果R探针3'端以A/T结尾，可以将其转移至L探针5'端
        if r_seq[-1] in ['A', 'T']:
            new_l = r_seq[-1] + l_seq
            new_r = r_seq[:-1]
            new_l_tm = calculate_tm(new_l)
            new_r_tm = calculate_tm(new_r)
            tm_diff = abs(new_l_tm - new_r_tm)
            
            if tm_diff <= 3.0:
                return new_l, new_r, f"R->L {r_seq[-1]}"
            else:
                return l_seq, r_seq, '跳过(TM>3)'
        else:
            # R探针3'端以C/G结尾，放弃调整
            return l_seq, r_seq, '跳过(L探针5端3个以上C/G且R端CG)'
    
    # 计算原始TM差异
    if l_tm is not None and r_tm is not None:
        original_tm_diff = abs(l_tm - r_tm)
    else:
        original_tm_diff = 0
    
    # 最大允许TM差异
    max_tm_diff = 3.0
    
    if gc_count == 1:
        # L探针5'端单个C/G开头
        l_first = l_seq[0]
        
        correction_options = []
        
        # 方案1：将R探针3'端的A或T转移至L探针5'端最前面（优先方案）
        if r_seq[-1] in ['A', 'T']:
            new_l = r_seq[-1] + l_seq
            new_r = r_seq[:-1]
            new_l_tm = calculate_tm(new_l)
            new_r_tm = calculate_tm(new_r)
            tm_diff = abs(new_l_tm - new_r_tm)
            correction_options.append({
                'new_l': new_l,
                'new_r': new_r,
                'note': f"R->L {r_seq[-1]}",
                'tm_diff': tm_diff,
                'priority': 1
            })
        
        # 方案2：将L探针5'端的C/G转移至R探针3'端
        new_l = l_seq[1:]
        new_r = r_seq + l_first
        new_l_tm = calculate_tm(new_l)
        new_r_tm = calculate_tm(new_r)
        tm_diff = abs(new_l_tm - new_r_tm)
        correction_options.append({
            'new_l': new_l,
            'new_r': new_r,
            'note': f"L->R {l_first}",
            'tm_diff': tm_diff,
            'priority': 2
        })
        
        # 选择TM差异最小且不超过3℃的方案
        valid_options = [opt for opt in correction_options if opt['tm_diff'] <= max_tm_diff]
        
        if valid_options:
            # 优先选择方案1（从R转移A/T到L）
            for opt in valid_options:
                if opt.get('priority', 2) == 1:
                    return opt['new_l'], opt['new_r'], opt['note']
            # 如果方案1无效，选择TM差异最小的方案
            best_opt = min(valid_options, key=lambda x: x['tm_diff'])
            return best_opt['new_l'], best_opt['new_r'], best_opt['note']
        else:
            # 所有方案TM差异都超过3℃，放弃调整
            return l_seq, r_seq, '跳过(TM>3)'
    
    elif gc_count == 2:
        # L探针5'端两个C/G开头
        # 检查第3个碱基是否为A/T
        if len(l_seq) >= 3 and l_seq[2] in ['A', 'T']:
            l_first_two = l_seq[:2]
            
            correction_options = []
            
            # 方案1：将R探针3'端的A或T转移至L探针5'端最前面（优先方案）
            if r_seq[-1] in ['A', 'T']:
                new_l = r_seq[-1] + l_seq
                new_r = r_seq[:-1]
                new_l_tm = calculate_tm(new_l)
                new_r_tm = calculate_tm(new_r)
                tm_diff = abs(new_l_tm - new_r_tm)
                correction_options.append({
                    'new_l': new_l,
                    'new_r': new_r,
                    'note': f"R->L {r_seq[-1]}",
                    'tm_diff': tm_diff,
                    'priority': 1
                })
            
            # 方案2：将L探针5'端两个C/G转移至R探针3'端
            new_l = l_seq[2:]
            new_r = r_seq + l_first_two
            new_l_tm = calculate_tm(new_l)
            new_r_tm = calculate_tm(new_r)
            tm_diff = abs(new_l_tm - new_r_tm)
            correction_options.append({
                'new_l': new_l,
                'new_r': new_r,
                'note': f"L->R {l_first_two}",
                'tm_diff': tm_diff,
                'priority': 2
            })
            
            # 选择TM差异最小且不超过3℃的方案
            valid_options = [opt for opt in correction_options if opt['tm_diff'] <= max_tm_diff]
            
            if valid_options:
                # 优先选择方案1（从R转移A/T到L）
                for opt in valid_options:
                    if opt.get('priority', 2) == 1:
                        return opt['new_l'], opt['new_r'], opt['note']
                # 如果方案1无效，选择TM差异最小的方案
                best_opt = min(valid_options, key=lambda x: x['tm_diff'])
                return best_opt['new_l'], best_opt['new_r'], best_opt['note']
            else:
                # 所有方案TM差异都超过3℃，放弃调整
                return l_seq, r_seq, '跳过(TM>3)'
        else:
            # L探针第3个碱基不是A/T，无法调整
            return l_seq, r_seq, '跳过(非A/T)'
    
    return l_seq, r_seq, ''

def select_probes(probe_pairs: List[Dict[str, Union[int, str]]], 
                  dual_probe_data: Dict[str, str], 
                  probe_scores: Dict[int, float], 
                  gene_name: str, 
                  bp_id: str, 
                  target_count: int) -> List[Dict[str, Union[int, str, float]]]:
    if len(probe_pairs) == 0:
        return []
    
    if len(probe_pairs) <= target_count:
        for p in probe_pairs:
            l_probe_id = f"{gene_name}-{p['num']}-L"
            r_probe_id = f"{gene_name}-{p['num']}-R-{bp_id}"
            l_seq = dual_probe_data.get(l_probe_id, '')
            r_seq = dual_probe_data.get(r_probe_id, '')
            p['priority'] = get_nucleotide_priority(l_seq, r_seq)
            p['score'] = probe_scores.get(p['num'], float('inf'))
        return sorted(probe_pairs, key=lambda x: (x['L_start'], x['priority']))
    
    for p in probe_pairs:
        l_probe_id = f"{gene_name}-{p['num']}-L"
        r_probe_id = f"{gene_name}-{p['num']}-R-{bp_id}"
        l_seq = dual_probe_data.get(l_probe_id, '')
        r_seq = dual_probe_data.get(r_probe_id, '')
        p['priority'] = get_nucleotide_priority(l_seq, r_seq)
        p['score'] = probe_scores.get(p['num'], float('inf'))
    
    probes_by_pos = sorted(probe_pairs, key=lambda x: x['L_start'])
    first_start = probes_by_pos[0]['L_start']
    last_end = probes_by_pos[-1]['R_end']
    gene_length = last_end - first_start
    
    interval = gene_length / target_count
    
    bins = []
    for i in range(target_count):
        bin_start = first_start + interval * i
        bin_end = first_start + interval * (i + 1)
        bins.append((bin_start, bin_end))
    
    probe_pairs_by_num = sorted(probe_pairs, key=lambda x: x['num'])
    
    selected = []
    used_nums = set()
    
    for bin_idx in range(target_count):
        if len(selected) >= target_count:
            break
        
        bin_start, bin_end = bins[bin_idx]
        
        candidates_in_bin = []
        for p in probe_pairs_by_num:
            if p['num'] in used_nums:
                continue
            
            probe_center = (p['L_start'] + p['R_end']) / 2
            if bin_start <= probe_center <= bin_end:
                candidates_in_bin.append(p)
        
        if not candidates_in_bin:
            for p in probe_pairs_by_num:
                if p['num'] in used_nums:
                    continue
                candidates_in_bin.append(p)
        
        candidates_in_bin.sort(key=lambda x: (x['score'], x['priority'], x['num']))
        
        for p in candidates_in_bin:
            if check_gap_requirement(p, selected):
                selected.append(p)
                used_nums.add(p['num'])
                break
    
    if len(selected) >= target_count:
        selected_sorted = sorted(selected, key=lambda x: (x['L_start'], x['priority']))
        return selected_sorted
    
    remaining = target_count - len(selected)
    remaining_probes = []
    for p in probe_pairs:
        if p['num'] in used_nums:
            continue
        remaining_probes.append(p)
    
    remaining_probes.sort(key=lambda x: (x['score'], x['priority'], x['num']))
    
    for p in remaining_probes:
        if len(selected) >= target_count:
            break
        if check_gap_requirement(p, selected):
            selected.append(p)
            used_nums.add(p['num'])
    
    if len(selected) < target_count:
        remaining = target_count - len(selected)
        remaining_probes = []
        for p in probe_pairs:
            if p['num'] in used_nums:
                continue
            remaining_probes.append(p)
        
        remaining_probes.sort(key=lambda x: (x['num']))
        
        for p in remaining_probes:
            selected.append(p)
            used_nums.add(p['num'])
            remaining -= 1
            if remaining == 0:
                break
    
    selected_sorted = sorted(selected, key=lambda x: (x['L_start'], x['priority']))
    
    return selected_sorted


# ---------------------------------------------------------------------------
# 以下由原 main() 拆分而来：文件名改为参数，结果可写文件或返回内存流。
# 原脚本遇到问题时静默跳过；这里照旧跳过，同时把原因记入 warnings 供页面/命令行显示。
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = ('Name', 'GeneID', 'BP_ID', 'TPM')
DEFAULT_OUTPUT_NAME = 'selected_probes_v3.xlsx'

# 提示模板（中文为源文本，网页按界面语言翻译；占位符与 locales/messages.en.json 一致）
WARNINGS = {
    'tpm_not_number': '{name}：TPM “{value}” 不是数字，按 0 处理（{pairs} 对）。',
    'duplicate_name': '基因表中 {name} 出现 {count} 次，使用最后一行。',
    'duplicate_gene_id': 'GeneID {gene_id} 在基因表中出现 {count} 次：{names}。',
    'folder_missing': '{name}：没有同名的结果文件夹，已跳过。',
    'files_missing': '{name}：缺少 dual_probe.bed 或 dual_probe_dual.csv，已跳过。',
    'no_probe_pairs': '{name}：dual_probe.bed 中没有完整的 L/R 探针对，已跳过。',
    'fewer_pairs': '{name}：目标 {target} 对，只有 {selected} 对可选。',
    'missing_sequences': '{name}：{count} 条探针在 dual_probe_dual.csv 中找不到（请核对 BP_ID {bp_id}）。',
    'folder_not_in_table': '结果文件夹 {name} 不在基因表中，已跳过。',
}


# 页面展示的 TPM 分档：每档取下限代入 get_target_pairs()，展示与算法始终一致
TPM_TIERS = (('< 1', 0.0), ('1 – 5', 1.0), ('5 – 10', 5.0), ('10 – 30', 10.0),
             ('30 – 50', 30.0), ('50 – 100', 50.0), ('≥ 100', 100.0))


def tpm_rule_rows():
    return [(label, get_target_pairs(tpm)) for label, tpm in TPM_TIERS]


class ProbeSelectorError(ValueError):
    """基因表或输入文件无法使用。"""


def warning(code, **values):
    return {'code': code, 'values': values}


def format_warning(item, translate=None):
    template = WARNINGS[item['code']]
    if translate is not None:
        return translate(template, **item['values'])
    return template.format(**item['values'])


def parse_tpm_table(text: str) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, Any]]]:
    """读取基因信息表（原 main() 开头部分）。返回 (tpm_data, warnings)。"""
    reader = csv.DictReader(io.StringIO(text.lstrip('﻿')))
    columns = [(c or '').strip() for c in (reader.fieldnames or [])]
    missing = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing:
        raise ProbeSelectorError('基因表缺少必需列：' + ', '.join(missing) + '（需要 Name, GeneID, BP_ID, TPM）。')
    reader.fieldnames = columns

    tpm_data = {}
    warnings = []
    name_counts = Counter()
    gene_names = {}
    for row in reader:
        name = (row['Name'] or '').strip()
        if not name:
            continue
        gene_id = (row['GeneID'] or '').strip()
        bp_id = (row['BP_ID'] or '').strip()
        tpm_val = (row['TPM'] or '').strip()
        try:
            tpm_val = float(tpm_val)
        except ValueError:
            warnings.append(warning('tpm_not_number', name=name, value=tpm_val, pairs=get_target_pairs(0.0)))
            tpm_val = 0.0
        tpm_data[name] = {'TPM': tpm_val, 'GeneID': gene_id, 'BP_ID': bp_id}
        name_counts[name] += 1
        if gene_id:
            gene_names.setdefault(gene_id, []).append(name)

    for name, count in name_counts.items():
        if count > 1:
            warnings.append(warning('duplicate_name', name=name, count=count))
    for gene_id, names in gene_names.items():  # analyze_duplicates.py 的检查
        if len(names) > 1:
            warnings.append(warning('duplicate_gene_id', gene_id=gene_id, count=len(names), names=', '.join(names)))
    if not tpm_data:
        raise ProbeSelectorError('基因表中没有基因（Name 列为空）。')
    return tpm_data, warnings


def load_tpm_table(filepath: str):
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        return parse_tpm_table(f.read())


def discover_gene_dirs(base_dir: str) -> Dict[str, str]:
    """原脚本的规则：base_dir 下每个子文件夹是一个基因，文件夹名即基因名。"""
    return {name: os.path.join(base_dir, name) for name in sorted(os.listdir(base_dir))
            if os.path.isdir(os.path.join(base_dir, name))}


def collect_gene_results(gene_dirs: Dict[str, str], tpm_data: Dict[str, Dict[str, Any]]):
    """对每个基因挑选并修正探针（原 main() 主循环）。返回 (gene_results, warnings)。"""
    gene_results = []
    warnings = []

    for folder_name in sorted(set(gene_dirs) - set(tpm_data)):
        warnings.append(warning('folder_not_in_table', name=folder_name))
    for folder_name in tpm_data:
        if folder_name not in gene_dirs:
            warnings.append(warning('folder_missing', name=folder_name))

    for folder_name, folder_path in gene_dirs.items():
        if folder_name not in tpm_data:
            continue

        dual_probe_bed = os.path.join(folder_path, 'dual_probe.bed')
        dual_probe_dual_csv = os.path.join(folder_path, 'dual_probe_dual.csv')
        dual_probe_csv = os.path.join(folder_path, 'dual_probe.csv')
        dual_probe_detailed = os.path.join(folder_path, 'dual_probe_detailed.txt')

        if not os.path.exists(dual_probe_bed) or not os.path.exists(dual_probe_dual_csv):
            warnings.append(warning('files_missing', name=folder_name))
            continue

        probe_pairs = parse_dual_probe_bed(dual_probe_bed)
        dual_probe_data = parse_dual_probe_dual_csv(dual_probe_dual_csv)

        # 解析dual_probe.csv获取原始靶标序列和TM值
        probe_tm_data = {}
        if os.path.exists(dual_probe_csv):
            probe_tm_data = parse_dual_probe_csv(dual_probe_csv)

        probe_scores = {}
        if os.path.exists(dual_probe_detailed):
            probe_scores = parse_dual_probe_detailed(dual_probe_detailed)

        if not probe_pairs:
            warnings.append(warning('no_probe_pairs', name=folder_name))
            continue

        tpm_val = tpm_data[folder_name]['TPM']
        gene_id = tpm_data[folder_name]['GeneID']
        bp_id = tpm_data[folder_name]['BP_ID']
        target_count = get_target_pairs(tpm_val)

        selected = select_probes(probe_pairs, dual_probe_data, probe_scores, folder_name, bp_id, target_count)

        probes_for_gene = []
        missing_sequences = 0
        for probe in selected:
            probe_num = probe['num']
            probe_score = probe_scores.get(probe_num, 'N/A')

            l_probe_id = f"{folder_name}-{probe_num}-L"
            r_probe_id = f"{folder_name}-{probe_num}-R-{bp_id}"

            l_seq = dual_probe_data.get(l_probe_id, '')
            r_seq = dual_probe_data.get(r_probe_id, '')

            # 获取原始TM值
            l_tm = None
            r_tm = None
            if probe_num in probe_tm_data:
                if probe_tm_data[probe_num]['L']:
                    l_tm = probe_tm_data[probe_num]['L']['tm']
                if probe_tm_data[probe_num]['R']:
                    r_tm = probe_tm_data[probe_num]['R']['tm']

            # 修正探针序列（考虑TM值）
            corrected_l, corrected_r, correction_note = correct_probe_sequence(l_seq, r_seq, l_tm, r_tm)

            if l_probe_id in dual_probe_data:
                probes_for_gene.append({
                    'ProbeID': l_probe_id,
                    'Sequence': corrected_l,
                    'Score': probe_score,
                    'CorrectionNote': correction_note
                })
            else:
                missing_sequences += 1

            if r_probe_id in dual_probe_data:
                probes_for_gene.append({
                    'ProbeID': r_probe_id,
                    'Sequence': corrected_r,
                    'Score': probe_score,
                    'CorrectionNote': ''
                })
            else:
                missing_sequences += 1

        if len(selected) < target_count:
            warnings.append(warning('fewer_pairs', name=folder_name, target=target_count, selected=len(selected)))
        if missing_sequences:
            warnings.append(warning('missing_sequences', name=folder_name, count=missing_sequences, bp_id=bp_id or '-'))

        gene_results.append({
            'GeneName': folder_name,
            'GeneID': gene_id,
            'TPM': tpm_val,
            'BP_ID': bp_id,
            'TargetCount': target_count,
            'SelectedPairs': len(selected),
            'Corrected': sum(1 for p in probes_for_gene if p['CorrectionNote'] and not p['CorrectionNote'].startswith('跳过')),
            'Probes': probes_for_gene
        })

    return gene_results, warnings


def build_workbook(gene_results) -> 'Workbook':
    """生成与原脚本相同版式的 Selected Probes 工作表。"""
    # 延迟导入：缺少 openpyxl 时只影响本功能，不影响整个网站启动
    from openpyxl import Workbook
    from openpyxl.styles import Alignment

    wb = Workbook()
    ws = wb.active
    ws.title = "Selected Probes"

    headers = ['GeneInfo', 'ProbeName', 'ProbeID', 'Sequence', 'Length', 'Score', 'Correction']
    ws.append(headers)

    current_row = 2

    for gene in gene_results:
        gene_name = gene['GeneName']
        gene_id = gene['GeneID']
        tpm_val = gene['TPM']
        bp_id = gene['BP_ID']
        target_count = gene['TargetCount']
        probes = gene['Probes']

        num_rows = len(probes)

        try:
            bp_num = int(bp_id.replace('BP_', ''))
        except ValueError:  # 原脚本在此处报错退出；非 BP_数字 格式时原样保留
            bp_num = bp_id.replace('BP_', '')
        probe_name = f"{gene_name}_{target_count}X_BP{bp_num}"

        middle_row = current_row + num_rows // 2 - 1

        for i, probe in enumerate(probes):
            row = current_row + i

            ws.cell(row=row, column=1, value='')
            ws.cell(row=row, column=2, value='')

            if row == middle_row - 1:
                ws.cell(row=row, column=1, value=gene_name)
            elif row == middle_row:
                ws.cell(row=row, column=1, value=gene_id)
                ws.cell(row=row, column=2, value=probe_name)
            elif row == middle_row + 1:
                ws.cell(row=row, column=1, value=tpm_val)

            ws.cell(row=row, column=3, value=probe['ProbeID'])
            ws.cell(row=row, column=4, value=probe['Sequence'])
            ws.cell(row=row, column=5, value=len(probe['Sequence']) if probe['Sequence'] else 0)
            ws.cell(row=row, column=6, value=probe['Score'])
            ws.cell(row=row, column=7, value=probe.get('CorrectionNote', ''))

        current_row += num_rows

    for col in ws.columns:
        max_length = 0
        column = col[0].column_letter
        for cell in col:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(str(cell.value))
            except:
                pass
        adjusted_width = max_length + 2
        ws.column_dimensions[column].width = adjusted_width

    for row in ws.iter_rows():
        for cell in row:
            cell.alignment = Alignment(horizontal='left', vertical='center')

    return wb


def workbook_bytes(wb) -> bytes:
    buffer = io.BytesIO()
    wb.save(buffer)
    return buffer.getvalue()


def run(tpm_file: str, base_dir: str, output_file: str):
    """命令行流程：等同原 main()，路径由参数给出。"""
    tpm_data, warnings = load_tpm_table(tpm_file)
    gene_results, more = collect_gene_results(discover_gene_dirs(base_dir), tpm_data)
    build_workbook(gene_results).save(output_file)
    return gene_results, warnings + more
