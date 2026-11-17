"""检查基因表中的重复 GeneID（ProbeSelector 辅助脚本，原 analyze_duplicates.py）。

用法：python -m scripts.analyze_duplicates [基因表.csv]   （默认 RLK_660.csv，需含 GeneID, Name, Class, TPM 列）
"""
import argparse
import csv
from collections import Counter


def main(argv=None):
    parser = argparse.ArgumentParser(description='检查基因表中的重复 GeneID')
    parser.add_argument('table', nargs='?', default='RLK_660.csv', help='基因表 CSV（默认 RLK_660.csv）')
    args = parser.parse_args(argv)

    with open(args.table, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = []
        gene_ids = []
        for idx, row in enumerate(reader, start=2):
            gene_id = row['GeneID'].strip()
            gene_ids.append(gene_id)
            rows.append((idx, gene_id, row['Name'], row['Class'], row['TPM']))

    counter = Counter(gene_ids)
    duplicates = {k: v for k, v in counter.items() if v > 1}

    print(f'{args.table} 总行数（不含表头）:', len(gene_ids))
    print(f'{args.table} 唯一 GeneID 数:', len(counter))
    print(f'{args.table} 重复的 GeneID 数:', len(duplicates))
    print()
    print('--- 重复的 GeneID 列表（共', len(duplicates), '个）---')
    print('{:<15} {:<8} {}'.format('GeneID', '出现次数', '行号和详细信息'))
    print('-' * 80)

    for gene_id in sorted(duplicates.keys()):
        count = duplicates[gene_id]
        occurrences = [(r[0], r[2], r[3], r[4]) for r in rows if r[1] == gene_id]
        info = '; '.join(['行{}({},{},{})'.format(r[0], r[1], r[2], r[3]) for r in occurrences])
        print('{:<15} {:<8} {}'.format(gene_id, count, info))


if __name__ == '__main__':
    main()
