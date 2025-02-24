import json
import copy
import argparse

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

def select_probe_sets(data, max_selected=5):
    """
    主处理函数：筛选探针组合并添加 keep 和 reason 字段
    Args:
        data: 输入的探针数据
        max_selected: 要选择的探针组合数量，默认为5
    """
    probe_sets = data['probe_sets']
    
    # 提取每个探针组合的探针区间，存储以便后续使用
    probe_intervals_dict = {}
    for probe_set in probe_sets:
        probe_intervals_dict[probe_set['id']] = get_probe_intervals(probe_set)
    
    # 按 quality_score 排序，quality_score 越低特异性越好
    sorted_probe_sets = sorted(probe_sets, key=lambda x: float(x['evaluation']['quality_score']))
    
    # 初始化已选择的探针组合和所有已选择的探针区间
    selected_sets = []
    selected_probe_intervals = []
    
    # 选择满足条件的探针组合
    for probe_set in sorted_probe_sets:
        if len(selected_sets) >= max_selected:
            break
        new_intervals = probe_intervals_dict[probe_set['id']]
        # 检查新探针组合是否与已选择的探针组合距离太近
        if not is_set_too_close(new_intervals, selected_probe_intervals):
            selected_sets.append(probe_set)
            selected_probe_intervals.extend(new_intervals)
    
    # 确保选择了 5 个探针组合
    if len(selected_sets) < max_selected:
        print(f"Warning: Only {len(selected_sets)} probe sets were selected, less than required {max_selected}.")
    
    # 为每个探针组合添加 keep 和 reason 字段
    selected_ids = {probe_set['id'] for probe_set in selected_sets}
    for probe_set in probe_sets:
        probe_id = probe_set['id']
        if probe_id in selected_ids:
            probe_set['keep'] = True
            probe_set['reason'] = "Selected for optimal specificity and spacing."
        else:
            probe_set['keep'] = False
            # 查找距离太近的已选择探针组合
            closest_selected_id = find_closest_selected_set(probe_set, selected_sets, probe_intervals_dict)
            if closest_selected_id is not None:
                probe_set['reason'] = f"Too close to selected set {closest_selected_id}."
            else:
                probe_set['reason'] = "Not selected due to lower specificity compared to selected sets."
    
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