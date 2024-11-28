#!/usr/bin/env python3
"""
Filter probes based on k-mer matching with a genome sequence.
Usage: python filter_probe.py genome.fasta probe.txt [k_size]
"""

import sys
import os
import hashlib
import pickle
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO

def generate_kmers(sequence, k):
    """Generate k-mers from a sequence."""
    return {sequence[i:i+k] for i in range(len(sequence)-k+1)}

def get_cache_path(genome_file, k):
    """Generate cache file path based on genome file hash and k value."""
    # Calculate MD5 hash of the genome file
    md5_hash = hashlib.md5()
    with open(genome_file, 'rb') as f:
        # Read in chunks to handle large files
        for chunk in iter(lambda: f.read(4096), b''):
            md5_hash.update(chunk)
    
    # Create cache directory if it doesn't exist
    cache_dir = Path(os.path.dirname(os.path.abspath(__file__))) / 'cache'
    cache_dir.mkdir(exist_ok=True)
    
    # Create cache filename using hash and k value
    cache_file = cache_dir / f"kmer_cache_{md5_hash.hexdigest()}_k{k}.pkl"
    return cache_file

def read_genome_kmers(genome_file, k=9):
    """Read genome file and count k-mer frequencies.
    
    Returns a Counter object with k-mer frequencies.
    Uses caching to avoid recomputing k-mers for the same file.
    """
    cache_file = get_cache_path(genome_file, k)
    
    # Try to load from cache first
    if cache_file.exists():
        print(f"Loading k-mers from cache: {cache_file}", file=sys.stderr)
        try:
            with open(cache_file, 'rb') as f:
                kmer_counts = pickle.load(f)
            return kmer_counts
        except Exception as e:
            print(f"Error loading cache: {e}", file=sys.stderr)
    
    # Use Counter to count k-mer frequencies
    kmer_counts = Counter()
    
    print(f"Computing k-mers for {genome_file}", file=sys.stderr)
    # Count k-mers in both forward and reverse complement sequences
    with open(genome_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq).upper()
            rev_seq = str(record.seq.reverse_complement()).upper()
            
            # Count k-mers in forward sequence
            for i in range(len(seq)-k+1):
                kmer = seq[i:i+k]
                kmer_counts[kmer] += 1
                
            # Count k-mers in reverse complement sequence
            for i in range(len(rev_seq)-k+1):
                kmer = rev_seq[i:i+k]
                kmer_counts[kmer] += 1
    
    print(f"Total unique k-mers: {len(kmer_counts)}", file=sys.stderr)
    print(f"Max k-mer frequency: {max(kmer_counts.values())}", file=sys.stderr)
    
    # Save to cache
    try:
        print(f"Saving k-mers to cache: {cache_file}", file=sys.stderr)
        with open(cache_file, 'wb') as f:
            pickle.dump(kmer_counts, f)
    except Exception as e:
        print(f"Error saving cache: {e}", file=sys.stderr)
    
    return kmer_counts

def get_kmer_info(sequence, kmer_counts, k):
    """Get the highest frequency k-mer and its count from a sequence."""
    seq_kmers = generate_kmers(sequence.upper(), k)
    if not seq_kmers:
        return "NA", 0
    
    # Find the k-mer with highest frequency
    max_kmer = max(seq_kmers, key=lambda x: kmer_counts.get(x, 0))
    return max_kmer, kmer_counts.get(max_kmer, 0)

def process_probes(probe_file, kmer_counts, k):
    """Process probes and add k-mer information.
    
    Returns a list of probe lines with added k-mer information.
    """
    processed_probes = []
    with open(probe_file) as f:
        for line in f:
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            parts = line.split()
            if not parts:
                continue
                
            sequence = parts[0]
            max_kmer, count = get_kmer_info(sequence, kmer_counts, k)
            
            # Add k-mer and count information
            processed_probes.append(f"{line}\t{max_kmer}\t{count}")
    
    return processed_probes

def main():
    if len(sys.argv) < 3:
        print("Usage: python filter_probe.py genome.fasta probe.txt [k_size]")
        sys.exit(1)

    genome_file = sys.argv[1]
    probe_file = sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 9

    # Get k-mer frequencies
    kmer_counts = read_genome_kmers(genome_file, k)
    
    # Process probes and add k-mer information
    processed_probes = process_probes(probe_file, kmer_counts, k)
    
    # Print header
    print(f"#Original_Content\tMax_Kmer\tKmer_Count")
    
    # Output processed probes
    for probe in processed_probes:
        print(probe)

if __name__ == "__main__":
    main()
