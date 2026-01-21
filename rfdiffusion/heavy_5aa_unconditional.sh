#!/bin/bash

# Set your RFDiffusion directory path here
RFDIFF_DIR="/srv/scratch/z5358491/RFdiffusion/5DK3_scripts"
PDB_FILE="5DK3.pdb"
OUTPUT_BASE="5dk3_outputs"

# Download PDB if needed
if [ ! -f "$PDB_FILE" ]; then
    echo "Downloading 5DK3.pdb..."
    wget https://files.rcsb.org/download/5DK3.pdb -O $PDB_FILE
fi

# Create output directory
mkdir -p $OUTPUT_BASE

cd $RFDIFF_DIR

echo "============================================================================"
echo "5DK3 ANTIBODY Single chain, fixed-length linker"
echo "============================================================================"
echo ""
echo "What this does:"
echo "- Takes only the heavy chain (B)"
echo "- Keeps residues 1-122 (variable region)"
echo "- Adds exactly 5 new amino acids (unconditional)"
echo "- Keeps residues 123-444 (constant region)"
echo "- This creates a simple linear split with new linker"
echo ""

pwd
python ../scripts/run_inference.py inference.input_pdb=$PDB_FILE inference.output_prefix=$OUTPUT_BASE/heavy_5aa_unconditional/run inference.num_designs=3 contigmap.contigs='[B1-122/5/B123-226/B233-444]'
# Residue 230 has no coordinate data (hinge region), and stops with the AssertionError. Need to skip over the unmodeled region. Sometimes in highly flexible regions like the hinge of an antibody, could not be reliably modeled from the experimental data.

echo "Output: Check $OUTPUT_BASE/heavy_5aa_unconditional*.pdb"
echo "You should see 3 designs with a 5aa linker between VH and CH1"
