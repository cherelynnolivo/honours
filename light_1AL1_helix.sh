#!/bin/bash

# Set your RFDiffusion directory path here
RFDIFF_DIR="/srv/scratch/z5358491/RFdiffusion/5DK3_scripts"
PDB_FILE="5DK3.pdb"
HELIX_FILE="1AL1.pdb"
OUTPUT_BASE="5dk3_outputs"

# Download PDB if needed
if [ ! -f "$PDB_FILE" ]; then
    echo "Downloading 5DK3.pdb..."
    wget https://files.rcsb.org/download/5DK3.pdb -O $PDB_FILE
fi

if [ ! -f "$HELIX_FILE" ]; then
    echo "Downloading 1AL1.pdb (helix template)..."
    wget https://files.rcsb.org/download/1AL1.pdb -O $HELIX_FILE
fi

# Create output directory
mkdir -p $OUTPUT_BASE

cd $RFDIFF_DIR

echo "============================================================================"
echo "5DK3 ANTIBODY single light chain with 1AL1 helix linker"
echo "============================================================================"
echo ""
echo "What this does:"
echo "- Takes the light chain (A)"
echo "- Keeps residues within the variable region (1-114)"
echo "- Inserts the 1AL1 helix structure as a linker"
echo "- Keeps residues within the constant region (115-218)"
echo "- This creates a split in the light chain and adds a helix linker"
echo ""

pwd
python ../scripts/run_inference.py inference.input_pdb=$PDB_FILE inference.output_prefix=$OUTPUT_BASE/light_1AL1_helix/run inference.num_designs=3 contigmap.contigs='[A1-114/0 B1-13/0 A115-218]' pdb_contig=$HELIX_FILE

echo "Output: Check $OUTPUT_BASE/light_1AL1_helix*.pdb"
echo "You should see 3 designs with a simple helix between VL and CL of the light chain"
