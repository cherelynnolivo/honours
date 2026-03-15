#!/usr/bin/env python3
import sys
from Bio.PDB import MMCIFParser

if len(sys.argv) < 2:
    print("Usage: python count_residues.py <file.cif>")
    sys.exit(1)

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("protein", sys.argv[1])

total = 0
for model in structure:
    for chain in model:
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        count = len(residues)
        total += count
        print(f"Chain {chain.id}: {count} residues")

print(f"\nTotal: {total} residues")
