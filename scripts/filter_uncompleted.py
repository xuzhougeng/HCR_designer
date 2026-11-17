"""从全量基因表中剔除已合成基因并去重，输出待做表（ProbeSelector 辅助脚本，原 filter_uncompleted.py）。

用法：
    python -m scripts.filter_uncompleted [--all RLK_660.csv] [--synthesized RLK_101-258_Synthesized.csv]
                                         [--output RLK_Uncompleted.csv] [--duplicates-report RLK_Duplicates_Report.csv]
全量表按 GeneID 列匹配，已合成表按 “Gene ID” 列匹配；默认文件名与原脚本相同。
"""
import argparse
import csv
from collections import Counter


def main(argv=None):
    parser = argparse.ArgumentParser(description='从全量基因表中剔除已合成基因，输出待做表')
    parser.add_argument('--all', dest='rlk660_path', default='RLK_660.csv', help='全量基因表（列 GeneID）')
    parser.add_argument('--synthesized', dest='synthesized_path', default='RLK_101-258_Synthesized.csv',
                        help='已合成基因表（列 Gene ID）')
    parser.add_argument('--output', dest='output_path', default='RLK_Uncompleted.csv', help='待做表输出路径')
    parser.add_argument('--duplicates-report', dest='duplicates_report_path', default='RLK_Duplicates_Report.csv',
                        help='重复 GeneID 报告输出路径')
    args = parser.parse_args(argv)
    rlk660_path = args.rlk660_path
    synthesized_path = args.synthesized_path
    output_path = args.output_path
    duplicates_report_path = args.duplicates_report_path

    all_rows = []
    all_gene_ids = []
    with open(rlk660_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            gene_id = row['GeneID'].strip()
            all_gene_ids.append(gene_id)
            all_rows.append(row)

    gene_counter = Counter(all_gene_ids)
    duplicate_gene_ids = {k for k, v in gene_counter.items() if v > 1}

    print(f"=== {rlk660_path} 分析报告 ===")
    print(f"总行数: {len(all_rows)}")
    print(f"唯一 GeneID 数: {len(gene_counter)}")
    print(f"重复 GeneID 数: {len(duplicate_gene_ids)}")

    if duplicate_gene_ids:
        print("\n重复的基因列表:")
        for gene_id in sorted(duplicate_gene_ids):
            count = gene_counter[gene_id]
            occurrences = [i+2 for i, gid in enumerate(all_gene_ids) if gid == gene_id]
            print(f"  {gene_id}: 出现 {count} 次，行号: {', '.join(map(str, occurrences))}")

        with open(duplicates_report_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['GeneID', '出现次数', '行号'])
            for gene_id in sorted(duplicate_gene_ids):
                count = gene_counter[gene_id]
                occurrences = [i+2 for i, gid in enumerate(all_gene_ids) if gid == gene_id]
                writer.writerow([gene_id, count, ', '.join(map(str, occurrences))])
        print(f"\n重复基因详细报告已保存到: {duplicates_report_path}")

    synthesized_gene_ids = set()
    with open(synthesized_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row['Gene ID'].strip()
            if gene_id:
                synthesized_gene_ids.add(gene_id)

    print(f"\n=== 合成状态分析 ===")
    print(f"已合成基因数: {len(synthesized_gene_ids)}")

    seen_gene_ids = set()
    uncompleted_rows = []
    for row in all_rows:
        gene_id = row['GeneID'].strip()
        if gene_id in synthesized_gene_ids:
            continue
        if gene_id in seen_gene_ids:
            continue
        seen_gene_ids.add(gene_id)
        uncompleted_rows.append(row)

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(uncompleted_rows)

    print(f"未合成基因数: {len(seen_gene_ids)}")
    print(f"输出文件行数: {len(uncompleted_rows)}")
    print(f"\n未合成基因表格已保存到: {output_path}")


if __name__ == '__main__':
    main()
