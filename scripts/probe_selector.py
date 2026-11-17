"""命令行版 ProbeSelector：按 TPM 挑选双探针并输出 Excel（原 probe_filter.py 的 main()）。

用法（在项目根目录）：
    python -m scripts.probe_selector --genes RLK359-458.csv --input-dir 结果目录 [-o selected_probes_v3.xlsx]

结果目录下每个基因一个文件夹，文件夹名与基因表 Name 列一致（即双探针批量结果 ZIP 解压后的结构）。
"""
import argparse
import os
import sys

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.dual_probe.probe_selector import DEFAULT_OUTPUT_NAME, format_warning, run


def main(argv=None):
    parser = argparse.ArgumentParser(description='按 TPM 从双探针候选中挑选探针对并输出 Excel')
    parser.add_argument('--genes', required=True, help='基因信息表 CSV（必需列 Name, GeneID, BP_ID, TPM）')
    parser.add_argument('--input-dir', default='.', help='包含各基因结果文件夹的目录（默认当前目录）')
    parser.add_argument('-o', '--output', help=f'输出 xlsx（默认 <input-dir>/{DEFAULT_OUTPUT_NAME}，重复运行会覆盖）')
    args = parser.parse_args(argv)

    output_file = args.output or os.path.join(args.input_dir, DEFAULT_OUTPUT_NAME)
    gene_results, warnings = run(args.genes, args.input_dir, output_file)
    for item in warnings:
        print('提示：' + format_warning(item))
    print(f"共 {len(gene_results)} 个基因")
    print(f"筛选完成，结果已保存到 {output_file}")


if __name__ == '__main__':
    main()
