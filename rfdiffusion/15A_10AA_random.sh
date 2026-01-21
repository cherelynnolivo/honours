#!/bin/bash

# Set your RFDiffusion directory path here
RFDIFF_DIR="/srv/scratch/z5358491/RFdiffusion/5DK3_scripts"
PDB_FILE="./5dk3_structures/5dk3_15A.pdb"
OUTPUT_BASE="5dk3_full_15A_10AA_outputs"

# Create output directory
mkdir -p $OUTPUT_BASE

cd $RFDIFF_DIR

pwd
python ../scripts/run_inference.py inference.input_pdb=$PDB_FILE \
  inference.output_prefix=$OUTPUT_BASE/run \
  inference.num_designs=10 \
  contigmap.contigs='[A1-218/0 B1-226/B233-238/10/B239-444/0 F1-218/0 G1-229/G236-238/10/G239-444/0]'
