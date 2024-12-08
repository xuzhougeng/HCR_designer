#!/usr/bin/env python3

import argparse
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Tuple
import sys

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    return (gc_count / total_length) * 100 if total_length > 0 else 0

def calculate_tm(sequence: str) -> float:
    """Calculate melting temperature of a sequence using nearest neighbor method."""
    seq = Seq(sequence)
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN4)  # Using nearest neighbor method

def process_fasta(input_file: str, output_file: str):
    """Process FASTA file and add GC content and Tm to headers."""
    modified_records = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        sequence = str(record.seq)
        gc_content = calculate_gc_content(sequence)
        tm = calculate_tm(sequence)
        
        # Modify the description to include GC content and Tm
        new_description = f"{record.description} GC={gc_content:.1f}% Tm={tm:.1f}°C"
        
        # Create new record with modified description
        new_record = SeqRecord(
            record.seq,
            id=record.id,
            description=new_description.replace(record.id + " ", "", 1)
        )
        modified_records.append(new_record)
    
    # Write modified records to output file
    SeqIO.write(modified_records, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description='Calculate GC content and Tm for sequences in FASTA file')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', type=str, default='output.fasta', 
                        help='Output FASTA file name (default: output.fasta)')
    
    args = parser.parse_args()
    
    # Process FASTA file
    process_fasta(args.input, args.output)
    print(f"\nModified FASTA file has been written to {args.output}")

if __name__ == "__main__":
    main()
