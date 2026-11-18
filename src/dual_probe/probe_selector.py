"""ProbeSelector：按 TPM 从双探针候选中挑选探针对，修正 L 探针 5' 端碱基，输出 Excel。

移植自 ProbeSelector_Share/probe_filter.py，并按生成它的提示词修正两处：
- 修正后 L/R 的 Tm 差（<=3 ℃）只按靶标序列计算（提示词 10.7），用与设计工具相同的最近邻 Tm；
  原脚本用含接头的完整探针和简化公式，几乎不会拦下修正。
- L 5' 端 >=3 个 C/G 且 R 3' 端为 C/G 的探针对只在其他探针对不够时使用（提示词第 9 步）；
  原脚本按打分优先，会在有更好选择时也选中它们。
其余（TPM 分档、区间等分、18 bp 间距、碱基优先级、修正方案）与原脚本一致；另修正任务名含下划线时
dual_probe.csv 的 Set_ID 解析失败。原 main() 拆成可传路径、可返回内存流的函数，基因名与文件夹名的
匹配规则不变。

输入（每个基因一个文件夹，文件夹名 = 基因表 Name 列）：
    dual_probe.bed            探针位置（制表符分隔，6 列以上）
    dual_probe_dual.csv       探针序列（列：Probe ID, Sequence）
    dual_probe.csv            可选，含 Set_ID / Primer_Type / Sequence / Tm / Length
    dual_probe_detailed.txt   可选，探针打分
    target.fasta              可选，基因序列；按 GeneID 匹配 TPM 表时用其标题
基因信息表 CSV 必需列：Name, GeneID, BP_ID, TPM；也可用探针数量推荐页导出的 gene_id, tpm 表。
"""
import os
import csv
import io
import re
from collections import Counter
from typing import Dict, List, Tuple, Any, Optional, Union

from src.common.sequence_utils import calculate_tm as target_tm  # 最近邻 Tm，与设计工具一致

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

# get_nucleotide_priority() 的最低档（8）：L 5' 端 >=3 个 C/G 且 R 3' 端为 C/G，提示词第 9 步只作兜底
LAST_RESORT_PRIORITY = 8


def check_gap_requirement(candidate: Dict[str, Union[int, str]], selected: List[Dict[str, Union[int, str]]]) -> bool:
    for sp in selected:
        gap_before = candidate['L_start'] - sp['R_end']
        gap_after = sp['L_start'] - candidate['R_end']
        if 0 <= gap_before < 18 or 0 <= gap_after < 18:
            return False
    return True

def correct_probe_sequence(l_seq: str, r_seq: str, l_tm: Optional[float] = None, r_tm: Optional[float] = None,
                           l_flank: int = 0, r_flank: int = 0) -> Tuple[str, str, str]:
    """
    修正探针序列，使L探针5'端以A或T开头，提高探针连接效率
    
    规则：
    1. L探针5'端单个C/G开头：将C/G转移至R探针3'端，或从R探针3'端转移A/T至L探针5'端
    2. L探针5'端两个C/G开头：将两个C/G转移至R探针3'端（需第3个碱基为A/T），或从R探针3'端转移A/T
    3. L探针5'端3个或更多C/G：放弃调整
    4. 调整后L/R靶标序列的TM值差异不能超过3℃（提示词 10.7：只算靶标序列，不含接头/桥探针；
       与设计工具相同的最近邻 Tm）
    
    参数：
    - l_seq, r_seq: L/R探针序列（完整合成序列；L 的 5' 端、R 的 3' 端为靶标序列，两者相接）
    - l_tm, r_tm: 原始TM值（可选，如果提供则用于验证调整后的TM差异）
    - l_flank, r_flank: L 探针 3' 端、R 探针 5' 端非靶标部分的长度，计算 Tm 前去掉
    
    返回：(修正后的L序列, 修正后的R序列, 修正说明)
    """
    def target_tm_diff(new_l, new_r):
        l_target = new_l[:len(new_l) - l_flank] if l_flank else new_l
        r_target = new_r[r_flank:]
        try:
            return abs(target_tm(l_target) - target_tm(r_target))
        except ValueError:  # 无法计算（如含非 ACGT 碱基）时不调整
            return float('inf')

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
            tm_diff = target_tm_diff(new_l, new_r)
            
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
            tm_diff = target_tm_diff(new_l, new_r)
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
        tm_diff = target_tm_diff(new_l, new_r)
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
                tm_diff = target_tm_diff(new_l, new_r)
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
            tm_diff = target_tm_diff(new_l, new_r)
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
            if p['num'] in used_nums or p['priority'] >= LAST_RESORT_PRIORITY:
                continue
            
            probe_center = (p['L_start'] + p['R_end']) / 2
            if bin_start <= probe_center <= bin_end:
                candidates_in_bin.append(p)
        
        if not candidates_in_bin:
            for p in probe_pairs_by_num:
                if p['num'] in used_nums or p['priority'] >= LAST_RESORT_PRIORITY:
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
        if p['num'] in used_nums or p['priority'] >= LAST_RESORT_PRIORITY:
            continue
        remaining_probes.append(p)
    
    remaining_probes.sort(key=lambda x: (x['score'], x['priority'], x['num']))
    
    for p in remaining_probes:
        if len(selected) >= target_count:
            break
        if check_gap_requirement(p, selected):
            selected.append(p)
            used_nums.add(p['num'])
    
    # 提示词第 9 步：前面仍选不够时，才使用 L 5' 端 >=3 个 C/G 且 R 3' 端为 C/G 的探针对（同样遵守间距）
    last_resort = sorted((p for p in probe_pairs if p['num'] not in used_nums and p['priority'] >= LAST_RESORT_PRIORITY),
                         key=lambda x: (x['score'], x['num']))
    for p in last_resort:
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
    'gene_id_unknown': '{name}：结果文件夹中没有 target.fasta，无法按 GeneID 匹配 TPM，已跳过。',
    'no_tpm_for_gene': '{name}：TPM 表中没有 {gene_id}，已跳过。',
    'duplicate_tpm_gene': 'TPM 表中 {gene_id} 出现 {count} 次，使用最后一行。',
}

# 第二种基因表：探针数量推荐页导出的 CSV（gene_id,tpm,probe_number），按 GeneID 匹配结果文件夹；
# probe_number 不使用，探针对数仍由 get_target_pairs() 决定。
RECOMMEND_COLUMNS = ('gene_id', 'tpm')

# 本流程合成序列的非靶标部分长度（scripts/create_split_probe.py）：
# L = 靶标 + TCTCGTTTTCTGAACTTAGC + CGGTATCAAG；R = CTGTTTAAGA + 随机碱基 + 桥探针改写（19 nt）+ 靶标
DEFAULT_L_FLANK = 30
DEFAULT_R_FLANK = 30


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
        raise ProbeSelectorError('基因表缺少必需列：' + ', '.join(missing)
                                 + '（需要 Name, GeneID, BP_ID, TPM；也可直接用探针数量推荐页导出的 gene_id, tpm 表）。')
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


def is_recommend_table(text: str) -> bool:
    """表头含 gene_id 与 tpm（不区分大小写）且没有 Name 列时，视为探针数量推荐页导出的表。"""
    header = next(csv.reader(io.StringIO(text.lstrip('\ufeff'))), [])
    columns = {(c or '').strip().lower() for c in header}
    return set(RECOMMEND_COLUMNS) <= columns and 'name' not in columns


def describe_gene_dir(folder_path: str, name: str) -> Tuple[Optional[str], Optional[str]]:
    """结果文件夹的 (GeneID, BP_ID)：GeneID 取 target.fasta 标题，BP_ID 取 dual_probe_dual.csv 中 R 探针编号后缀。"""
    gene_id = bp_id = None
    fasta = os.path.join(folder_path, 'target.fasta')
    if os.path.exists(fasta):
        with open(fasta, 'r', encoding='utf-8') as f:
            first = f.readline().strip()
        if first.startswith('>') and first[1:].split():
            gene_id = first[1:].split()[0]
    dual = os.path.join(folder_path, 'dual_probe_dual.csv')
    if os.path.exists(dual):
        pattern = re.compile(re.escape(name) + r'-\d+-R-(.+)')
        for probe_id in parse_dual_probe_dual_csv(dual):
            match = pattern.fullmatch(probe_id)
            if match:
                bp_id = match.group(1)
                break
    return gene_id, bp_id


def tpm_by_gene_id(text: str, gene_dirs: Dict[str, str]):
    """把推荐页的 TPM 表（每个基因一行，可能覆盖整个转录组）对应到结果文件夹。返回 (tpm_data, warnings)。"""
    reader = csv.DictReader(io.StringIO(text.lstrip('\ufeff')))
    columns = {(c or '').strip().lower(): c for c in (reader.fieldnames or [])}
    rows, counts = {}, Counter()
    for row in reader:
        gene_id = (row[columns['gene_id']] or '').strip()
        if gene_id:
            rows[gene_id.casefold()] = (row[columns['tpm']] or '').strip()
            counts[gene_id.casefold()] += 1
    if not rows:
        raise ProbeSelectorError('TPM 表中没有基因（gene_id 列为空）。')

    tpm_data, warnings = {}, []
    for name, folder_path in gene_dirs.items():
        gene_id, bp_id = describe_gene_dir(folder_path, name)
        if not gene_id:
            warnings.append(warning('gene_id_unknown', name=name))
            continue
        key = gene_id.casefold()
        if key not in rows:
            warnings.append(warning('no_tpm_for_gene', name=name, gene_id=gene_id))
            continue
        if counts[key] > 1:  # 只提示用到的基因，整张转录组表的其他行不逐条报告
            warnings.append(warning('duplicate_tpm_gene', gene_id=gene_id, count=counts[key]))
        try:
            tpm_val = float(rows[key])
        except ValueError:
            warnings.append(warning('tpm_not_number', name=name, value=rows[key], pairs=get_target_pairs(0.0)))
            tpm_val = 0.0
        tpm_data[name] = {'TPM': tpm_val, 'GeneID': gene_id, 'BP_ID': bp_id or ''}
    return tpm_data, warnings


def select_from_gene_dirs(table_text: str, gene_dirs: Dict[str, str]):
    """读取任一格式的基因表并对 gene_dirs 中的基因挑选。返回 (gene_results, warnings)。"""
    if is_recommend_table(table_text):
        tpm_data, warnings = tpm_by_gene_id(table_text, gene_dirs)
        gene_dirs = {name: path for name, path in gene_dirs.items() if name in tpm_data}
    else:
        tpm_data, warnings = parse_tpm_table(table_text)
    gene_results, more = collect_gene_results(gene_dirs, tpm_data)
    return gene_results, warnings + more


def target_flanks(l_seq: str, r_seq: str, arms) -> Tuple[int, int]:
    """L 探针 3' 端、R 探针 5' 端的非靶标长度：优先按 dual_probe.csv 中的原始靶标序列推出。"""
    l_flank = DEFAULT_L_FLANK if len(l_seq) > DEFAULT_L_FLANK else 0
    r_flank = DEFAULT_R_FLANK if len(r_seq) > DEFAULT_R_FLANK else 0
    if arms:
        if arms['L'] and arms['L']['sequence'] and l_seq.startswith(arms['L']['sequence']):
            l_flank = len(l_seq) - len(arms['L']['sequence'])
        if arms['R'] and arms['R']['sequence'] and r_seq.endswith(arms['R']['sequence']):
            r_flank = len(r_seq) - len(arms['R']['sequence'])
    return l_flank, r_flank


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

            # 修正探针序列（考虑TM值；只比较靶标序列的 Tm）
            l_flank, r_flank = target_flanks(l_seq, r_seq, probe_tm_data.get(probe_num))
            corrected_l, corrected_r, correction_note = correct_probe_sequence(l_seq, r_seq, l_tm, r_tm,
                                                                               l_flank, r_flank)

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
    """命令行流程：等同原 main()，路径由参数给出；基因表可为四列表或推荐页导出的 TPM 表。"""
    with open(tpm_file, 'r', encoding='utf-8-sig') as f:
        table_text = f.read()
    gene_results, warnings = select_from_gene_dirs(table_text, discover_gene_dirs(base_dir))
    build_workbook(gene_results).save(output_file)
    return gene_results, warnings
