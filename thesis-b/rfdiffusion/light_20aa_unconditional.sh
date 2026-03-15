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
echo "5DK3 ANTIBODY single chain, fixed-length linker"
echo "============================================================================"
echo ""
echo "What this does:"
echo "- Takes the light chain (A)"
echo "- Keeps residues within the variable region"
echo "- Adds exactly 20 new amino acids (unconditional)"
echo "- Keeps residues within the constant region"
echo "- This creates a split in the light chain and adds a new linker"
echo ""

pwd
python ../scripts/run_inference.py inference.input_pdb=$PDB_FILE inference.output_prefix=$OUTPUT_BASE/light_20aa_unconditional/run inference.num_designs=3 contigmap.contigs='[A1-114/20/115-218]'

echo "Output: Check $OUTPUT_BASE/heavy_light_20aa_unconditional*.pdb"
echo "You should see 3 designs with a 20aa linker between VH and CH1 of the light chain"
