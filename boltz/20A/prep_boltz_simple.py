#!/usr/bin/env python3
"""
Simple Boltz YAML preparation script for 5DK3 ProteinMPNN outputs.
No external dependencies required (no pyyaml needed).

Usage:
    python prep_boltz_simple.py <distance> <aa_length>
    
Example:
    python prep_boltz_simple.py 10A 10AA
    
Output directory: ./boltz_inputs/{aa_length}/
"""

import os
import sys
import glob

# === FIXED CHAIN SEQUENCES (from 5dk3_cleaned.pdb) ===
CHAIN_A = "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
CHAIN_F = "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"

# === CONFIGURATION ===
MPNN_BASE = "../../../ProteinMPNN/5dk3_scripts"
OUTPUT_BASE = "./boltz_inputs"


def parse_fasta(filepath):
    """Parse FASTA file, return list of (header, sequence)."""
    entries = []
    with open(filepath, 'r') as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            entries.append((header, "".join(seq)))
    return entries


def write_yaml_manual(chain_a, chain_b, chain_f, chain_g, outpath):
    """
    Write Boltz YAML manually without pyyaml dependency.
    """
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {chain_a}
      msa: empty
  - protein:
      id: B
      sequence: {chain_b}
      msa: empty
  - protein:
      id: F
      sequence: {chain_f}
      msa: empty
  - protein:
      id: G
      sequence: {chain_g}
      msa: empty
"""
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        f.write(yaml_content)


def get_design_name(header, idx):
    """Extract sample number from header or generate name."""
    if "sample=" in header:
        for part in header.split(','):
            if "sample=" in part:
                return f"sample_{part.split('=')[1].strip()}"
    return "original" if idx == 0 else f"seq_{idx}"


def main():
    if len(sys.argv) < 3:
        print("Usage: python prep_boltz_simple.py <distance> <aa_length>")
        print("Example: python prep_boltz_simple.py 10A 10AA")
        print("\nAvailable distances: 10A, 15A, 20A, 30A, 40A, 50A")
        print("Available AA lengths: 10AA, 20AA, 30AA, 40AA, 50AA")
        print("\nOutput: ./boltz_inputs/{aa_length}/*.yaml")
        sys.exit(1)
    
    distance = sys.argv[1]  # e.g., "10A"
    aa_length = sys.argv[2]  # e.g., "10AA"
    
    # Build input path
    input_dir = os.path.join(MPNN_BASE, distance, f"{distance}_{aa_length}_outputs", "seqs")
    
    # Output directory includes aa_length prefix
    output_dir = os.path.join(OUTPUT_BASE, aa_length)
    
    if not os.path.exists(input_dir):
        print(f"Error: Input directory not found: {input_dir}")
        sys.exit(1)
    
    fasta_files = sorted(glob.glob(os.path.join(input_dir, "run_*.fa")))
    print(f"Found {len(fasta_files)} FASTA files in {input_dir}")
    print(f"Output directory: {output_dir}")
    
    count = 0
    for fasta_path in fasta_files:
        run_name = os.path.basename(fasta_path).replace(".fa", "")
        entries = parse_fasta(fasta_path)
        
        for i, (header, sequence) in enumerate(entries):
            # ProteinMPNN outputs: "designed_B/designed_G"
            parts = sequence.split('/')
            if len(parts) != 2:
                print(f"  Warning: Unexpected format in {run_name}, entry {i}")
                continue
            
            chain_b, chain_g = parts[0], parts[1]
            design_name = get_design_name(header, i)
            filename = f"{run_name}_{design_name}.yaml"
            outpath = os.path.join(output_dir, filename)
            
            write_yaml_manual(CHAIN_A, chain_b, CHAIN_F, chain_g, outpath)
            count += 1
    
    print(f"Created {count} YAML files in {output_dir}/")


if __name__ == "__main__":
    main()
