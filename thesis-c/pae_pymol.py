# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
# ]
# ///
"""
Generate a PyMOL colouring script from Boltz-2 PAE data.

Computes mean PAE between each residue and a target chain,
then outputs a .pml script that colours the structure by that value.

Usage:
    uv run pae_pymol.py <cif_path> --target <chain> [--output <script.pml>]

Examples:
    # Colour all residues by their mean PAE to chain G
    uv run pae_pymol.py ./predictions/5dk3_base/5dk3_base_model_0.cif --target G

    # Colour by mean PAE to chain A
    uv run pae_pymol.py ./predictions/5dk3_base/5dk3_base_model_0.cif --target A

Then in PyMOL:
    @pae_to_G.pml
"""

import sys
import os
import argparse
import numpy as np
from collections import OrderedDict


def parse_chains_from_cif(cif_path):
    """Parse chain IDs and residue counts from a Boltz-2 output CIF."""
    chains = OrderedDict()
    seen_residues = {}

    in_atom_site = False
    field_indices = {}
    field_count = 0

    with open(cif_path, 'r') as f:
        for line in f:
            line_stripped = line.strip()

            if line_stripped == 'loop_':
                in_atom_site = False
                field_indices = {}
                field_count = 0
                continue

            if line_stripped.startswith('_atom_site.'):
                in_atom_site = True
                field_name = line_stripped.split('.')[1]
                field_indices[field_name] = field_count
                field_count += 1
                continue

            if in_atom_site and (line_stripped.startswith('_') or line_stripped.startswith('#') or line_stripped == ''):
                if line_stripped.startswith('_') and not line_stripped.startswith('_atom_site.'):
                    in_atom_site = False
                if line_stripped.startswith('#') or line_stripped == '':
                    in_atom_site = False
                continue

            if not in_atom_site:
                continue

            tokens = line_stripped.split()
            if len(tokens) < field_count:
                continue

            group_pdb_idx = field_indices.get('group_PDB')
            if group_pdb_idx is not None:
                if tokens[group_pdb_idx] != 'ATOM':
                    continue

            chain_idx = field_indices.get('label_asym_id')
            seq_idx = field_indices.get('label_seq_id')

            if chain_idx is None or seq_idx is None:
                continue

            chain_id = tokens[chain_idx]
            seq_id = tokens[seq_idx]

            key = (chain_id, seq_id)
            if key not in seen_residues:
                seen_residues[key] = True
                if chain_id not in chains:
                    chains[chain_id] = []
                chains[chain_id].append(int(seq_id))

    # Convert to OrderedDict of {chain_id: sorted list of seq_ids}
    for chain_id in chains:
        chains[chain_id] = sorted(set(chains[chain_id]))

    return chains


def find_pae_npz(cif_path):
    """Find the PAE NPZ file corresponding to a CIF."""
    import glob

    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    cif_basename = os.path.basename(cif_path).replace('.cif', '')

    pae_npz = os.path.join(cif_dir, f'pae_{cif_basename}.npz')
    if os.path.exists(pae_npz):
        return pae_npz

    if '_model_' in cif_basename:
        model_num = cif_basename.split('_model_')[-1]
    else:
        model_num = '0'

    patterns = [
        os.path.join(cif_dir, f'*pae*model_{model_num}.npz'),
    ]
    for pattern in patterns:
        matches = glob.glob(pattern)
        matches = [m for m in matches if not os.path.basename(m).startswith('._')]
        if matches:
            return matches[0]

    return None


def main():
    parser = argparse.ArgumentParser(description='Generate PyMOL PAE colouring script')
    parser.add_argument('cif_path', help='Path to Boltz-2 output CIF')
    parser.add_argument('--target', required=True, help='Target chain ID to compute mean PAE against')
    parser.add_argument('--output', default=None, help='Output .pml script path (default: pae_to_<chain>.pml)')
    parser.add_argument('--vmax', type=float, default=30.0, help='Max PAE for colour scale (default: 30)')
    args = parser.parse_args()

    if not os.path.exists(args.cif_path):
        print(f"ERROR: CIF not found: {args.cif_path}")
        sys.exit(1)

    # Parse chains
    chains = parse_chains_from_cif(args.cif_path)
    chain_ids = list(chains.keys())
    print(f"Chains in CIF: {[(c, len(chains[c])) for c in chain_ids]}")

    if args.target not in chains:
        print(f"ERROR: Target chain '{args.target}' not found. Available: {chain_ids}")
        sys.exit(1)

    # Find and load PAE
    npz_path = find_pae_npz(args.cif_path)
    if npz_path is None:
        print("ERROR: Could not find PAE NPZ file.")
        sys.exit(1)
    print(f"PAE file: {npz_path}")

    data = np.load(npz_path)
    pae = data['pae']
    if pae.ndim == 3:
        pae = pae[0]
    print(f"PAE shape: {pae.shape}")

    # Build matrix index ranges for each chain (in CIF file order)
    chain_ranges = OrderedDict()
    pos = 0
    for chain_id in chain_ids:
        n = len(chains[chain_id])
        chain_ranges[chain_id] = (pos, pos + n)
        pos += n

    target_start, target_end = chain_ranges[args.target]
    target_cols = np.arange(target_start, target_end)

    # For each residue in the structure, compute mean PAE to the target chain
    # Use the symmetric measure: mean of pae[i, target] and pae[target, i]
    residue_pae = {}
    for chain_id in chain_ids:
        c_start, c_end = chain_ranges[chain_id]
        seq_ids = chains[chain_id]

        for local_idx, seq_id in enumerate(seq_ids):
            matrix_idx = c_start + local_idx

            # Mean PAE: this residue scored against target chain residues
            pae_forward = pae[matrix_idx, target_cols].mean()
            # Mean PAE: target chain residues scored against this residue
            pae_reverse = pae[target_cols, matrix_idx].mean()
            # Symmetric average
            mean_pae = (pae_forward + pae_reverse) / 2.0

            residue_pae[(chain_id, seq_id)] = mean_pae

    # Generate PyMOL script
    obj_name = os.path.basename(args.cif_path).replace('.cif', '')
    output_path = args.output or f"pae_to_{args.target}.pml"

    with open(output_path, 'w') as f:
        f.write(f"# PAE colouring script\n")
        f.write(f"# Residues coloured by mean PAE to chain {args.target}\n")
        f.write(f"# Blue = low PAE (confident), Red = high PAE (uncertain)\n")
        f.write(f"# Generated from: {os.path.basename(args.cif_path)}\n\n")

        # Set all B-factors to 0 first
        f.write(f"alter {obj_name}, b=0\n\n")

        # Set B-factor to mean PAE value for each residue
        for (chain_id, seq_id), mean_pae in residue_pae.items():
            f.write(f"alter {obj_name} and chain {chain_id} and resi {seq_id}, b={mean_pae:.2f}\n")

        f.write(f"\nrebuild\n")
        f.write(f"spectrum b, blue_white_red, {obj_name}, minimum=0, maximum={args.vmax}\n")
        f.write(f"\n# Optional: show chain boundaries\n")
        for chain_id in chain_ids:
            f.write(f"# select chain_{chain_id}, {obj_name} and chain {chain_id}\n")

        f.write(f"\nprint('Coloured by mean PAE to chain {args.target}')\n")
        f.write(f"print('Blue = confident relative position, Red = uncertain')\n")

    print(f"\nSaved PyMOL script: {output_path}")
    print(f"\nIn PyMOL, load your structure then run:")
    print(f"  @{output_path}")

    # Print summary stats
    print(f"\nMean PAE to chain {args.target} by chain:")
    for chain_id in chain_ids:
        vals = [residue_pae[(chain_id, s)] for s in chains[chain_id]]
        print(f"  {chain_id}: mean={np.mean(vals):.1f}  min={np.min(vals):.1f}  max={np.max(vals):.1f}")


if __name__ == '__main__':
    main()
