#!/bin/bash

folder_with_pdbs="/srv/scratch/z5358491/RFdiffusion/5DK3_scripts/thesis_c/insertion/10A/10A_40aa_outputs"
output_dir="/srv/scratch/z5358491/ProteinMPNN/thesis_c/insertion/10A/10A_40AA_outputs"
mkdir -p $output_dir
mkdir -p $output_dir/seqs

# Parse the PDB into JSONL format
path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
python /srv/scratch/z5358491/ProteinMPNN//helper_scripts/parse_multiple_chains.py \
  --input_path=$folder_with_pdbs \
  --output_path=$path_for_parsed_chains

# Specify which chains to design
chains_to_design="B G"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
python /srv/scratch/z5358491/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
  --input_path=$path_for_parsed_chains \
  --output_path=$path_for_assigned_chains \
  --chain_list "$chains_to_design"

path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
python /srv/scratch/z5358491/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py \
  --input_path=$path_for_parsed_chains \
  --output_path=$path_for_fixed_positions \
  --chain_list "$chains_to_design" \
  --specify_non_fixed \
  --position_list "233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272, 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269"

# Run ProteinMPNN
python /srv/scratch/z5358491/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path $path_for_parsed_chains \
  --chain_id_jsonl $path_for_assigned_chains \
  --fixed_positions_jsonl $path_for_fixed_positions \
  --out_folder $output_dir \
  --num_seq_per_target 8 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
