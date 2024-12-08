
def create_triplet_probe(primer_sets, bridge_probe):
    """
    创建三合一探针, 用于TCR

    Parameters:
    -----------
    primer_sets : list
        探针组合列表，每个元素包含三个探针信息 (L, M, R)
    bridge_probe : str
        桥接探针序列，长度必须为19个碱基

    Returns:
    --------
    list: 转换后的三合一探针列表，每个元素包含：
        L: L + N + brigde_probe[0:16] + (bridge_probe[17:19]的互补序列) + AAGATA
        M: ACATTA + M
        R: R + TAATGTTATCTT
    """
    if not bridge_probe or len(bridge_probe) != 19:
        raise ValueError("桥接探针必须为19个碱基长度")

    # 定义碱基互补对应关系
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    
    # 随机碱基选择
    random_base = random.choice(['A', 'T', 'C', 'G'])
    
    # 获取bridge_probe最后两个碱基的互补序列
    bridge_end_complement = ''.join(complement[base] for base in bridge_probe[17:19])
    
    triplet_probes = []
    for left, middle, right in primer_sets:
        # 使用primer_sets中的序列信息
        L = left[0]  # 获取左引物序列
        M = middle[0]  # 获取中间引物序列
        R = right[0]  # 获取右引物序列
        
        # 构建三个探针
        L_probe = L + random_base + bridge_probe[0:17] + bridge_end_complement + "AAGATA"
        M_probe = "ACATTA" + M
        R_probe = R + "TAATGTTATCTT"
        
        triplet_probes.append((L_probe, M_probe, R_probe))


        # 创建 README 文件
        readme_path = os.path.join(output_dir, "readme.txt")
        with open(readme_path, 'w') as f:
            f.write(f"TCR Probe Design Results\n")
            f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Probe Length: {min_length}-{max_length} bp\n")
            f.write(f"- GC Content: {gc_min}%-{gc_max}%\n")
            f.write(f"- Melting Temperature: {tm_min}°C-{tm_max}°C\n")
            f.write(f"- Minimum Gap: {min_gap} bp\n")
            f.write(f"- Minimum Complementary Length: {min_complementary_length} bp\n")
            if blast_db:
                f.write(f"- Reference Genome: {os.path.basename(blast_db)}\n")
        files.append(readme_path)

    return triplet_probes

