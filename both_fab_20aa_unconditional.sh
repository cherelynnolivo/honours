#!/bin/bash

RFDIFF_DIR="/srv/scratch/z5358491/RFdiffusion/5DK3_scripts"
PDB_FILE="5DK3.pdb"
OUTPUT_BASE="5dk3_outputs"

if [ ! -f "$PDB_FILE" ]; then
    echo "Downloading 5DK3.pdb..."
    wget https://files.rcsb.org/download/5DK3.pdb -O $PDB_FILE
fi

mkdir -p $OUTPUT_BASE
cd $RFDIFF_DIR

echo "============================================================================"
echo "5DK3 ANTIBODY - Both Fab (Heavy + Light chains together)"
echo "============================================================================"
echo ""
echo "What this does:"
echo "- Heavy chain (B): VH + linker + CH1"
echo "- Light chain (A): VL + linker + CL"
echo "- Both chains modified simultaneously"
echo ""

pwd
python ../scripts/run_inference.py \
    inference.input_pdb=$PDB_FILE \
    inference.output_prefix=$OUTPUT_BASE/both_fab_5aa_unconditional/run \
    inference.num_designs=4 \
    contigmap.contigs='[A1-114/5/A115-218/0 B1-122/5/B123-226/B233-444/0 F1-114/5/F115-218/0 G1-122/5/G123-226/G236-444]'

echo "Output: Check $OUTPUT_BASE/both_fab_5aa_unconditional*.pdb"