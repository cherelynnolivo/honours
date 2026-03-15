#!/bin/bash
folder_with_pdbs="../../../RFdiffusion/5DK3_scripts/5dk3_full_10A_20AA_outputs/"
output_dir="./10A_20AA_outputs"
mkdir -p $output_dir
# Parse the PDB into JSONL format
path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
python ../../helper_scripts/parse_multiple_chains.py \
  --input_path=$folder_with_pdbs \
  --output_path=$path_for_parsed_chains
# Specify which chains to design
chains_to_design="B G"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
python ../../helper_scripts/assign_fixed_chains.py \
  --input_path=$path_for_parsed_chains \
  --output_path=$path_for_assigned_chains \
  --chain_list "$chains_to_design"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
python ../../helper_scripts/make_fixed_positions_dict.py \
  --input_path=$path_for_parsed_chains \
  --output_path=$path_for_fixed_positions \
  --chain_list "$chains_to_design" \
  --specify_non_fixed \
  --position_list "233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252, 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252"
# Run ProteinMPNN
python ../../protein_mpnn_run.py \
  --jsonl_path $path_for_parsed_chains \
  --chain_id_jsonl $path_for_assigned_chains \
  --fixed_positions_jsonl $path_for_fixed_positions \
  --out_folder $output_dir \
  --num_seq_per_target 8 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
