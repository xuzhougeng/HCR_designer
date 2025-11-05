import json
import copy
import argparse
from collections import defaultdict

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

def find_closest_selected_set(probe_set, selected_sets, probe_intervals_dict):
    """
    找到与给定探针组合距离太近的已选择探针组合
    返回第一个发现的太近的探针组合的 id
    """
    new_intervals = get_probe_intervals(probe_set)
    for selected_set in selected_sets:
        selected_id = selected_set['id']
        selected_intervals = probe_intervals_dict[selected_id]
        if is_set_too_close(new_intervals, selected_intervals):
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
    基于探针的脱靶情况（mismatch_stats）
    返回特异性得分（越小越好）
    """
    specificity_score = 0
    
    for probe_key in ['left_probe', 'middle_probe', 'right_probe']:
        probe = probe_set.get(probe_key)
        if not probe or 'blast_results' not in probe or 'mismatch_stats' not in probe['blast_results']:
            continue
        
        mismatch_stats = probe['blast_results']['mismatch_stats']
        
        # 对完全匹配的脱靶位置进行更严格的惩罚
        exact_matches = int(mismatch_stats.get("0", 0))
        one_mismatch = int(mismatch_stats.get("1", 0))
        
        # 完全匹配的脱靶权重要高
        specificity_score += exact_matches * 10 + one_mismatch * 2
    
    return specificity_score

def select_probe_sets(data, max_selected=5):
    """
    主处理函数：筛选探针组合并添加 keep 和 reason 字段
    优先考虑间隔合适和覆盖度，然后考虑特异性
    
    Args:
        data: 输入的探针数据
        max_selected: 要选择的探针组合数量，默认为5
    """
    probe_sets = data['probe_sets']
    
    # 提取每个探针组合的探针区间，存储以便后续使用
    probe_intervals_dict = {}
    for probe_set in probe_sets:
        probe_intervals_dict[probe_set['id']] = get_probe_intervals(probe_set)
    
    # 计算每个探针组的位置范围和特异性得分
    probe_set_info = []
    for probe_set in probe_sets:
        start, end = calculate_position_range(probe_set)
        if start is None or end is None:
            continue
            
        # 同时考虑质量得分和特异性评估
        quality_score = float(probe_set['evaluation']['quality_score'])
        specificity_score = evaluate_specificity(probe_set)

        # 综合得分 = 质量得分 + 特异性得分（两者都是越小越好）
        combined_score = quality_score + specificity_score

        # 如果Left探针3'端具有强碱基连接(第一个碱基是G/C)，降低综合得分以优先选择
        # strong_junction=True时减5分，使其更容易被选中
        if probe_set['evaluation'].get('strong_junction', False):
            combined_score -= 5.0

        probe_set_info.append({
            'id': probe_set['id'],
            'start': start,
            'end': end,
            'quality_score': quality_score,
            'specificity_score': specificity_score,
            'combined_score': combined_score,
            'probe_set': probe_set
        })
    
    # 按位置排序
    probe_set_info.sort(key=lambda x: x['start'])
    
    # 使用动态规划找出最佳组合
    # 尝试所有可能的起始探针组
    best_selection = []
    best_coverage_score = 0
    best_combined_score = float('inf')
    
    for start_idx in range(len(probe_set_info)):
        # 从每个可能的起点开始选择
        current_selection = []
        selected_intervals = []
        
        # 贪心选择，从当前起点遍历所有探针组
        for i in range(start_idx, len(probe_set_info)):
            current_probe_info = probe_set_info[i]
            current_intervals = probe_intervals_dict[current_probe_info['id']]
            
            # 检查是否与已选择的探针组距离太近
            if not is_set_too_close(current_intervals, selected_intervals):
                current_selection.append(current_probe_info)
                selected_intervals.extend(current_intervals)
                
                # 如果已经选择了足够数量的探针组，就停止
                if len(current_selection) >= max_selected:
                    break
        
        # 评估当前选择
        if len(current_selection) == max_selected:
            coverage_score = calculate_coverage_score([info['probe_set'] for info in current_selection])
            combined_score_sum = sum(info['combined_score'] for info in current_selection)
            
            # 优先考虑覆盖度，其次考虑综合得分
            if coverage_score > best_coverage_score or (coverage_score == best_coverage_score and combined_score_sum < best_combined_score):
                best_selection = current_selection
                best_coverage_score = coverage_score
                best_combined_score = combined_score_sum
    
    # 如果没有找到足够数量的探针组，就使用贪心算法找出最多的探针组
    if not best_selection:
        # 按质量分数排序
        sorted_probe_sets = sorted(probe_set_info, key=lambda x: x['combined_score'])
        
        # 初始化已选择的探针组和所有已选择的探针区间
        selected_sets = []
        selected_probe_intervals = []
        
        # 选择满足条件的探针组合
        for probe_info in sorted_probe_sets:
            current_intervals = probe_intervals_dict[probe_info['id']]
            # 检查新探针组合是否与已选择的探针组合距离太近
            if not is_set_too_close(current_intervals, selected_probe_intervals):
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
        if probe_id in selected_ids:
            probe_set['keep'] = True
            probe_set['reason'] = "Selected for optimal coverage and specificity."
        else:
            probe_set['keep'] = False
            # 查找距离太近的已选择探针组合
            closest_selected_id = find_closest_selected_set(
                probe_set, 
                [info['probe_set'] for info in best_selection], 
                probe_intervals_dict
            )
            if closest_selected_id is not None:
                probe_set['reason'] = f"Too close to selected set {closest_selected_id}."
            else:
                probe_set['reason'] = "Not selected due to lower specificity or suboptimal coverage."
    
    return data

def main():
    parser = argparse.ArgumentParser(description='选择探针组合')
    parser.add_argument('input_file', help='输入的JSON文件路径')
    parser.add_argument('-n', '--num_probes', type=int, default=5, help='要选择的探针组合数量（默认：5）')
    
    args = parser.parse_args()
    
    # 加载数据
    with open(args.input_file, 'r') as f:
        data = json.load(f)
    
    # 处理数据
    updated_data = select_probe_sets(data, args.num_probes)
    
    # 输出结果
    print(json.dumps(updated_data, indent=2))

if __name__ == "__main__":
    main()