#!/usr/bin/env python3
import random
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import argparse
from typing import List, Tuple
import re
from functools import lru_cache
import multiprocessing
from multiprocessing import Pool
from functools import partial

# 定义碱基互补对照表
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

@lru_cache(maxsize=1024)
def has_hairpin(sequence: str, hairpin_len: int = 4, min_matches: int = 3) -> bool:
    """
    Check if sequence can form hairpin structure.
    Uses caching to avoid recalculating for the same sequence.
    """
    # Get the first and last N bases
    front = sequence[:hairpin_len]
    back = sequence[-hairpin_len:]
    
    # Reverse complement the back sequence
    back_rc = ''.join(COMPLEMENT[base] for base in reversed(back))
    
    # Count matching bases
    matches = sum(1 for f, b in zip(front, back_rc) if f == b)
    
    return matches >= min_matches

def has_poly_n(sequence: str, max_repeat: int = 4) -> bool:
    """Check if sequence contains poly-N (repeated nucleotides)."""
    # 使用更快的字符串搜索
    return any(base * max_repeat in sequence for base in 'ATCG')

def generate_random_sequences(length: int, count: int = 2000) -> List[str]:
    """
    Generate multiple random DNA sequences of given length.
    
    Args:
        length: Length of sequences to generate
        count: Number of sequences to generate
    
    Returns:
        List of random DNA sequences
    """
    sequences = []
    bases = ['A', 'T', 'G', 'C']
    
    # 直接生成指定数量的随机序列
    for _ in range(count):
        seq = ''.join(random.choices(bases, k=length))
        sequences.append(seq)
    
    return sequences

@lru_cache(maxsize=1024)
def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence. Uses caching to avoid recalculating."""
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)

@lru_cache(maxsize=1024)
def calculate_tm(sequence: str) -> float:
    """Calculate melting temperature. Uses caching to avoid recalculating."""
    return mt.Tm_NN(Seq(sequence))

def process_sequence_chunk(args, min_gc, max_gc, max_poly_n, min_tm, max_tm, hairpin_len, min_hairpin_matches):
    """Process a chunk of sequences in parallel."""
    sequences, chunk_stats = args
    valid_chunk = []
    
    for seq in sequences:
        chunk_stats['total'] += 1
        
        # 快速GC含量预检查
        gc_count = seq.count('G') + seq.count('C')
        gc_content = gc_count / len(seq)
        if not (min_gc <= gc_content <= max_gc):
            chunk_stats['gc'] += 1
            continue
        
        # 检查连续碱基
        if has_poly_n(seq, max_poly_n):
            chunk_stats['poly_n'] += 1
            continue
            
        # 计算精确的Tm值
        temp = calculate_tm(seq)
        if not (min_tm <= temp <= max_tm):
            chunk_stats['tm'] += 1
            continue
            
        # 检查发卡结构
        if has_hairpin(seq, hairpin_len, min_hairpin_matches):
            chunk_stats['hairpin'] += 1
            continue
            
        # 序列通过所有检查
        chunk_stats['valid'] += 1
        valid_chunk.append((seq, temp, gc_content))
    
    return valid_chunk, chunk_stats

def filter_valid_probes(sequences: List[str], target_count: int = 1000, max_poly_n: int = 4,
                     min_tm: float = 53.0, max_tm: float = 56.0,
                     min_gc: float = 0.4, max_gc: float = 0.6,
                     hairpin_len: int = 4, min_hairpin_matches: int = 3,
                     num_cores: int = 8) -> Tuple[List[Tuple[str, float, float]], dict]:
    """Filter sequences based on various criteria using multiprocessing.
    
    Args:
        sequences: List of DNA sequences to filter
        target_count: Number of valid probes to find before stopping (not used anymore)
        max_poly_n: Maximum number of consecutive identical nucleotides allowed
        min_tm: Minimum melting temperature allowed
        max_tm: Maximum melting temperature allowed
        min_gc: Minimum GC content allowed
        max_gc: Maximum GC content allowed
        hairpin_len: Length of sequence to check for hairpin structures
        min_hairpin_matches: Minimum number of matches to consider as hairpin
        num_cores: Number of CPU cores to use for parallel processing (default: 8)
    
    Returns:
        Tuple containing list of valid probes and statistics dictionary
    """
    # 初始化统计信息
    stats = {'total': 0, 'poly_n': 0, 'gc': 0, 'tm': 0, 'hairpin': 0, 'valid': 0}
    
    # 使用指定的核心数，但不超过系统可用核心数
    available_cores = multiprocessing.cpu_count()
    num_cores = min(num_cores, available_cores)
    chunk_size = max(1, len(sequences) // (num_cores * 4))  # 每个进程处理多个chunk以平衡负载
    
    # 将序列分成chunks
    sequence_chunks = []
    chunk_stats_list = []
    for i in range(0, len(sequences), chunk_size):
        chunk = sequences[i:min(i + chunk_size, len(sequences))]
        chunk_stats = {'total': 0, 'poly_n': 0, 'gc': 0, 'tm': 0, 'hairpin': 0, 'valid': 0}
        sequence_chunks.append((chunk, chunk_stats))
        chunk_stats_list.append(chunk_stats)

    # 创建部分函数，固定其他参数
    process_func = partial(process_sequence_chunk,
                         min_gc=min_gc, max_gc=max_gc,
                         max_poly_n=max_poly_n,
                         min_tm=min_tm, max_tm=max_tm,
                         hairpin_len=hairpin_len,
                         min_hairpin_matches=min_hairpin_matches)

    # 使用进程池并行处理
    valid_probes = []
    with Pool(processes=num_cores) as pool:
        results = pool.map(process_func, sequence_chunks)
        
        # 合并结果，保留所有有效探针
        for chunk_valid_probes, chunk_stats in results:
            valid_probes.extend(chunk_valid_probes)
            for key in stats:
                stats[key] += chunk_stats[key]
    
    return valid_probes, stats

def load_existing_probes(input_file: str) -> List[str]:
    """Load existing probe sequences from a file."""
    probes = []
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # 移除空白字符和注释
                line = line.strip()
                if line and not line.startswith('#'):
                    # 如果行包含多列（如温度、GC含量等），只取第一列作为序列
                    probe = line.split()[0].strip()
                    probes.append(probe)
        print(f"Loaded {len(probes)} existing probes from {input_file}")
    except FileNotFoundError:
        print(f"Warning: Input file {input_file} not found")
    return probes

def main():
    parser = argparse.ArgumentParser(description='Generate DNA probes with specific criteria')
    parser.add_argument('-n', '--num_probes', type=int, default=1000,
                       help='Number of probes to generate')
    parser.add_argument('-l', '--length', type=int, default=25,
                       help='Length of each probe (final length will be this value)')
    parser.add_argument('-o', '--output', type=str, default='probes.txt',
                       help='Output file name (default: probes.txt)')
    parser.add_argument('-i', '--input', type=str,
                       help='Input file containing existing probes (optional)')
    parser.add_argument('-p', '--max_poly_n', type=int, default=4,
                       help='Maximum allowed consecutive identical nucleotides (default: 4)')
    parser.add_argument('--min_tm', type=float, default=53.0,
                       help='Minimum melting temperature in Celsius (default: 52.0)')
    parser.add_argument('--max_tm', type=float, default=56.0,
                       help='Maximum melting temperature in Celsius (default: 57.0)')
    parser.add_argument('--min_gc', type=float, default=0.4,
                       help='Minimum GC content (default: 0.4)')
    parser.add_argument('--max_gc', type=float, default=0.6,
                       help='Maximum GC content (default: 0.6)')
    parser.add_argument('--hairpin_len', type=int, default=4,
                       help='Length of sequence ends to check for hairpin structure (default: 4)')
    parser.add_argument('--min_hairpin_matches', type=int, default=3,
                       help='Minimum number of matching bases to consider as hairpin (default: 3)')
    parser.add_argument('--num_cores', type=int, default=8,
                       help='Number of CPU cores to use for parallel processing (default: 8)')
    
    args = parser.parse_args()
    
    print(f"Generating initial pool of sequences with length {args.length}...")
    print(f"Temperature range: {args.min_tm}-{args.max_tm}°C")
    print(f"GC content range: {args.min_gc*100:.0f}%-{args.max_gc*100:.0f}%")
    print(f"Checking hairpin with {args.hairpin_len} bases and minimum {args.min_hairpin_matches} matches")
    
    # 加载现有探针（如果提供了输入文件）
    existing_probes = []
    if args.input:
        existing_probes = load_existing_probes(args.input)
        # 验证现有探针的长度
        valid_existing = [p for p in existing_probes if len(p) == args.length]
        if len(valid_existing) != len(existing_probes):
            print(f"Warning: {len(existing_probes) - len(valid_existing)} probes were skipped due to incorrect length")
        existing_probes = valid_existing
    
    # 计算需要额外生成的探针数量
    remaining_probes = max(0, args.num_probes - len(existing_probes))
    
    valid_probes = []
    stats = {'total': 0, 'poly_n': 0, 'gc': 0, 'tm': 0, 'hairpin': 0, 'valid': 0}
    
    # 如果有现有探针，先验证它们
    if existing_probes:
        print(f"Validating {len(existing_probes)} existing probes...")
        existing_valid_probes, existing_stats = filter_valid_probes(
            existing_probes,
            target_count=len(existing_probes),
            max_poly_n=args.max_poly_n,
            min_tm=args.min_tm,
            max_tm=args.max_tm,
            min_gc=args.min_gc,
            max_gc=args.max_gc,
            hairpin_len=args.hairpin_len,
            min_hairpin_matches=args.min_hairpin_matches,
            num_cores=args.num_cores
        )
        valid_probes.extend(existing_valid_probes)
        for key in stats:
            stats[key] += existing_stats[key]
        
        print(f"Found {len(existing_valid_probes)} valid probes from existing sequences")
        remaining_probes = max(0, args.num_probes - len(valid_probes))
    
    # 如果还需要更多探针，生成新的
    if remaining_probes > 0:
        print(f"Generating {remaining_probes} additional probes...")
        # 直接生成随机序列，不考虑GC含量
        sequences = generate_random_sequences(args.length, count=remaining_probes * 1000)
        print(f"Generated {len(sequences)} sequences, filtering for valid probes...")
        
        new_valid_probes, new_stats = filter_valid_probes(
            sequences,
            max_poly_n=args.max_poly_n,
            min_tm=args.min_tm,
            max_tm=args.max_tm,
            min_gc=args.min_gc,
            max_gc=args.max_gc,
            hairpin_len=args.hairpin_len,
            min_hairpin_matches=args.min_hairpin_matches,
            num_cores=args.num_cores
        )
        valid_probes.extend(new_valid_probes)
        for key in stats:
            stats[key] += new_stats[key]
    
    if len(valid_probes) < args.num_probes:
        print(f"\nWarning: Only found {len(valid_probes)} valid probes out of {args.num_probes} requested")
    
    # 输出结果
    print(f"\nStatistics:")
    print(f"Total sequences processed: {stats['total']}")
    print(f"Failed poly-N filter: {stats['poly_n']}")
    print(f"Failed GC content filter: {stats['gc']}")
    print(f"Failed Tm filter: {stats['tm']}")
    print(f"Failed hairpin filter: {stats['hairpin']}")
    print(f"Valid probes: {stats['valid']}")
    
    print(f"\nWriting {len(valid_probes)} probes to {args.output}")
    with open(args.output, 'w') as f:
        for probe, tm, gc in valid_probes:
            f.write(f"{probe}\t{tm:.2f}\t{gc:.2f}\n")
    
    print(f"\nGeneration complete!")

if __name__ == "__main__":
    main()
