import json
import copy
import argparse
from collections import defaultdict
if __package__ in (None, ''):
    # Preserve the existing direct-script CLI as well as package imports.
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common.paired_specificity import paired_specificity_score, specificity_rejection

REJECTION_REASONS = {
    'specificity_not_evaluated': 'Not selected: paired specificity evidence is missing or incomplete.',
    'no_exact_on_target_pair': 'Not selected: no exact pair was found on the intended target.',
    'cross_species_exact_pair': 'Excluded: an adjacent full-length exact pair exists in another species.',
    'cross_species_near_pair': 'Excluded: a near-exact pair in another species is within the policy threshold.',
    'unclassified_exact_pair': 'Excluded: an exact pair exists on a non-target transcript with unclassified species.',
    'unclassified_near_pair': 'Excluded: a near-exact pair on an unclassified transcript is within the policy threshold.',
    'same_species_exact_pair': 'Excluded: an exact pair exists on a non-allowed transcript of the target species.',
    'same_species_near_pair': 'Excluded: a near-exact pair on a non-allowed same-species transcript is within the policy threshold.',
}

def get_probe_intervals(probe_set):
    """
    获取探针组合中所有探针的区间（start, end）
    """
    intervals = []
    for probe_key in ['left_probe', 'middle_probe', 'right_probe']:
        probe = probe_set.get(probe_key)
        if probe and 'position' in probe:
            start = probe['position']['start']
            end = probe['position']['end']
            intervals.append((start, end))
    return intervals

def intervals_too_close(interval1, interval2, min_distance=15):
    """
    检查两个区间是否距离小于 min_distance
    两个区间 [s1, e1] 和 [s2, e2] 距离 >= d 当且仅当：
    e1 < s2 - d 或 e2 < s1 - d
    """
    s1, e1 = interval1
    s2, e2 = interval2
    # 如果两个区间重叠或距离 < min_distance，则认为太近
    return not (e1 < s2 - min_distance or e2 < s1 - min_distance)

def is_set_too_close(new_intervals, selected_probe_intervals, min_distance=15):
    """
    检查新探针组合的探针是否与已选择的探针组合的探针距离太近
    """
    for new_interval in new_intervals:
        for selected_interval in selected_probe_intervals:
            if intervals_too_close(new_interval, selected_interval, min_distance):
                return True
    return False

def find_closest_selected_set(probe_set, selected_sets, probe_intervals_dict, min_distance=15):
    """
    找到与给定探针组合距离太近的已选择探针组合
    返回第一个发现的太近的探针组合的 id
    """
    new_intervals = get_probe_intervals(probe_set)
    for selected_set in selected_sets:
        selected_id = selected_set['id']
        selected_intervals = probe_intervals_dict[selected_id]
        if is_set_too_close(new_intervals, selected_intervals, min_distance):
            return selected_id
    return None

def calculate_position_range(probe_set):
    """
    计算探针组合的起始和结束位置范围
    """
    intervals = get_probe_intervals(probe_set)
    if not intervals:
        return None, None
    
    start_positions = [interval[0] for interval in intervals]
    end_positions = [interval[1] for interval in intervals]
    return min(start_positions), max(end_positions)

def calculate_coverage_score(selected_probe_sets):
    """
    计算探针组的覆盖度得分
    考虑探针组在序列全长的分布情况
    """
    if not selected_probe_sets:
        return 0
    
    # 获取所有探针组的位置范围
    positions = []
    for probe_set in selected_probe_sets:
        start, end = calculate_position_range(probe_set)
        if start is not None and end is not None:
            positions.append((start, end))
    
    if not positions:
        return 0
    
    # 按起始位置排序
    positions.sort()
    
    # 计算覆盖的范围
    min_pos = min([pos[0] for pos in positions])
    max_pos = max([pos[1] for pos in positions])
    total_range = max_pos - min_pos
    
    # 计算探针组的均匀分布得分
    if len(positions) <= 1:
        return 1
    
    # 计算相邻探针组之间的距离
    gaps = []
    for i in range(1, len(positions)):
        gap = positions[i][0] - positions[i-1][1]
        gaps.append(gap)
    
    # 理想情况下，所有间隔应该相等
    ideal_gap = total_range / (len(positions) - 1)
    gap_variance = sum((gap - ideal_gap) ** 2 for gap in gaps) / len(gaps)
    
    # 归一化得分 (越小越好)
    normalized_score = 1 / (1 + gap_variance / 10000)
    
    return normalized_score

def evaluate_specificity(probe_set):
    """
    评估探针组的特异性
    优先使用同一转录本上的双臂配对结果，旧格式使用单臂统计。
    返回用于排序的得分（越小越好）；缺少评估时返回 None。
    """
    if 'specificity' in probe_set:
        return paired_specificity_score(probe_set['specificity'])
    specificity_score = 0
    
    for probe_key in ['left_probe', 'middle_probe', 'right_probe']:
        probe = probe_set.get(probe_key)
        if not probe:
            continue
        if 'blast_results' not in probe or not probe['blast_results'].get('mismatch_stats'):
            return None
        
        mismatch_stats = probe['blast_results']['mismatch_stats']
        
        # 对完全匹配的脱靶位置进行更严格的惩罚
        exact_matches = int(mismatch_stats.get("0", mismatch_stats.get(0, 0)))
        one_mismatch = int(mismatch_stats.get("1", mismatch_stats.get(1, 0)))
        
        # 完全匹配的脱靶权重要高
        specificity_score += exact_matches * 10 + one_mismatch * 2
    
    return specificity_score

def select_probe_sets(data, max_selected=5, gc_junction_penalty=20.0, min_distance=15):
    """
    主处理函数：筛选探针组合并添加 keep 和 reason 字段

    策略：
    1. 要求特异性检查时，排除证据缺失、无目标完全配对及跨物种完全配对。
    2. 保留原有连接位点过滤，按实际质量分数与配对特异性分数排序。
    3. 按指定间距贪心选择，不保证最大覆盖度。

    Args:
        data: 输入的探针数据
        max_selected: 要选择的探针组合数量，默认为5
        gc_junction_penalty: Left探针G/C连接处的惩罚分数（硬性过滤时此参数无效）
        min_distance: 两个已选组合的臂之间至少保留的未占用碱基数
    """
    if max_selected < 0 or min_distance < 0:
        raise ValueError('max_selected and min_distance must be non-negative')
    probe_sets = data['probe_sets']
    if len({p['id'] for p in probe_sets}) != len(probe_sets):
        raise ValueError('Probe set IDs must be unique')
    rejections = {}
    data['selection_parameters'] = {'max_selected': max_selected,
                                    'min_distance': min_distance,
                                    'gc_junction_penalty': gc_junction_penalty}

    # 提取每个探针组合的探针区间，存储以便后续使用
    probe_intervals_dict = {}
    for probe_set in probe_sets:
        probe_intervals_dict[probe_set['id']] = get_probe_intervals(probe_set)

    # 计算每个探针组的位置范围和特异性得分
    probe_set_info = []
    for probe_set in probe_sets:
        # 硬性过滤：直接排除has_gc_junction=True的探针（前两个碱基都是G/C）
        specificity_score = evaluate_specificity(probe_set)
        probe_set['evaluation']['specificity_score'] = specificity_score
        if data.get('specificity_required'):
            rejection = specificity_rejection(probe_set.get('specificity'))
            if rejection:
                rejections[probe_set['id']] = rejection
                continue
        if probe_set['evaluation'].get('has_gc_junction', False):
            continue

        start, end = calculate_position_range(probe_set)
        if start is None or end is None:
            continue

        # 同时考虑质量得分和特异性评估
        quality_score = float(probe_set['evaluation'].get('quality_score') or 0)

        # 综合得分 = 质量得分 + 特异性得分（两者都是越小越好）
        combined_score = quality_score + (specificity_score or 0)

        # 软惩罚：如果第一个碱基是G/C（但第二个不是），增加一定惩罚
        # 这样会优先选择A/T开头的探针，但不会完全排除C/G单碱基开头的
        left_seq = probe_set['left_probe']['sequence']
        if len(left_seq) >= 1 and left_seq[0] in 'GC':
            combined_score += gc_junction_penalty  # 使用可配置的惩罚分数

        probe_set_info.append({
            'id': probe_set['id'],
            'start': start,
            'end': end,
            'quality_score': quality_score,
            'specificity_score': specificity_score,
            'combined_score': combined_score,
            'probe_set': probe_set
        })

    # 按综合得分排序（得分越低越优先）
    # 优先选择质量最好的探针，不考虑覆盖度
    probe_set_info.sort(key=lambda x: x['combined_score'])

    # 贪心选择：按得分顺序选择，确保探针间距离足够
    best_selection = []
    selected_intervals = []

    for probe_info in (probe_set_info if max_selected else []):
        current_intervals = probe_intervals_dict[probe_info['id']]

        # 检查是否与已选择的探针组距离太近
        if not is_set_too_close(current_intervals, selected_intervals, min_distance):
            best_selection.append(probe_info)
            selected_intervals.extend(current_intervals)

            # 如果已经选择了足够数量的探针组，就停止
            if len(best_selection) >= max_selected:
                break
    
    # 如果没有找到足够数量的探针组，就使用贪心算法找出最多的探针组
    if not best_selection and max_selected:
        # 按质量分数排序
        sorted_probe_sets = sorted(probe_set_info, key=lambda x: x['combined_score'])
        
        # 初始化已选择的探针组和所有已选择的探针区间
        selected_sets = []
        selected_probe_intervals = []
        
        # 选择满足条件的探针组合
        for probe_info in sorted_probe_sets:
            current_intervals = probe_intervals_dict[probe_info['id']]
            # 检查新探针组合是否与已选择的探针组合距离太近
            if not is_set_too_close(current_intervals, selected_probe_intervals, min_distance):
                selected_sets.append(probe_info)
                selected_probe_intervals.extend(current_intervals)
                
                # 如果已经选择了足够数量的探针组，就停止
                if len(selected_sets) >= max_selected:
                    break
        
        best_selection = selected_sets
    
    # 确保选择了尽可能多的探针组
    if len(best_selection) < max_selected:
        print(f"Warning: Only {len(best_selection)} probe sets were selected, less than required {max_selected}.")

    # 提取选择的探针组ID
    selected_ids = {info['id'] for info in best_selection}

    # 为每个探针组合添加 keep 和 reason 字段
    for probe_set in probe_sets:
        probe_id = probe_set['id']

        # 特异性排除优先于质量和距离，避免把明确的跨物种匹配隐藏为其他原因。
        if probe_id in rejections:
            probe_set['keep'] = False
            probe_set['reason_code'] = rejections[probe_id]
            probe_set['reason'] = REJECTION_REASONS.get(
                rejections[probe_id], f'Excluded by the target policy: {rejections[probe_id]}.')
        elif probe_set['evaluation'].get('has_gc_junction', False):
            probe_set['keep'] = False
            probe_set['reason_code'] = 'gc_junction'
            probe_set['reason'] = "Filtered: Left probe starts with GC/CC/CG/GG (inhibits SplintR ligation)."
        elif probe_id in selected_ids:
            probe_set['keep'] = True
            probe_set['reason_code'] = 'selected'
            probe_set['reason'] = ("Selected by quality and observed paired-match ranking; experimental specificity is unverified."
                                   if probe_set.get('specificity', {}).get('status') == 'completed'
                                   else "Selected by quality; paired specificity was not evaluated.")
        else:
            probe_set['keep'] = False
            # 查找距离太近的已选择探针组合
            closest_selected_id = find_closest_selected_set(
                probe_set,
                [info['probe_set'] for info in best_selection],
                probe_intervals_dict, min_distance
            )
            if closest_selected_id is not None:
                probe_set['reason_code'] = 'too_close'
                probe_set['reason'] = f"Too close to selected set {closest_selected_id}."
            else:
                probe_set['reason_code'] = 'selection_limit_or_rank'
                probe_set['reason'] = "Not selected within the requested count after ranking."
    
    return data

def main():
    parser = argparse.ArgumentParser(description='选择探针组合')
    parser.add_argument('input_file', help='输入的JSON文件路径')
    parser.add_argument('-n', '--num_probes', type=int, default=5, help='要选择的探针组合数量（默认：5）')
    parser.add_argument('--min-distance', type=int, default=15, help='已选组合之间的最少空闲碱基数（默认：15）')
    
    args = parser.parse_args()
    
    # 加载数据
    with open(args.input_file, 'r') as f:
        data = json.load(f)
    
    # 处理数据
    updated_data = select_probe_sets(data, args.num_probes, min_distance=args.min_distance)
    
    # 输出结果
    print(json.dumps(updated_data, indent=2))

if __name__ == "__main__":
    main()
