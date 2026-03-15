#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["numpy", "matplotlib"]
# ///
"""
Generate a combined PAE heatmap + pLDDT plot for a single Boltz-2 model.

Usage:
    python3 plot_model.py <predictions_dir> <model_number> [output.png]

Example:
    python3 plot_model.py 2a3d/replacement/5dk3_2a3d_replacement_run_3_outputs/boltz_results_5dk3_2a3d_replacement/predictions/5dk3_2a3d_replacement 0
    
    # Or with uv:
    uv run plot_model.py <predictions_dir> <model_number>
"""

import sys
import os
import json
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PALETTE = ["#06b6d4", "#22c55e", "#eab308", "#ec4899", "#a855f7", "#f97316"]


def find_file(directory, pattern):
    """Find a file matching pattern, excluding macOS resource forks."""
    matches = [f for f in glob.glob(os.path.join(directory, pattern))
               if not os.path.basename(f).startswith('._')]
    return matches[0] if matches else None


def parse_chains_from_cif(cif_path):
    """Parse chain boundaries from a Boltz-2 CIF file."""
    chains = []
    if not cif_path or not os.path.exists(cif_path):
        return chains

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    # Find _atom_site column headers
    col_names = []
    in_atom_site = False
    data_start = 0

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('_atom_site.'):
            in_atom_site = True
            col_names.append(stripped.split('.')[1])
        elif in_atom_site and not stripped.startswith('_atom_site.') and not stripped.startswith('#'):
            if stripped and not stripped.startswith('loop_'):
                data_start = i
                break

    if not col_names:
        return chains

    # Find column indices
    try:
        chain_col = col_names.index('auth_asym_id')
    except ValueError:
        try:
            chain_col = col_names.index('label_asym_id')
        except ValueError:
            return chains

    try:
        atom_col = col_names.index('label_atom_id')
    except ValueError:
        try:
            atom_col = col_names.index('auth_atom_id')
        except ValueError:
            return chains

    # Count CA atoms per chain (preserving order of appearance)
    chain_order = []
    chain_counts = {}

    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) <= max(chain_col, atom_col):
            continue
        if parts[0] in ('loop_', '#') or parts[0].startswith('_'):
            break

        atom_name = parts[atom_col]
        chain_id = parts[chain_col]

        if atom_name == 'CA':
            if chain_id not in chain_counts:
                chain_order.append(chain_id)
                chain_counts[chain_id] = 0
            chain_counts[chain_id] += 1

    for cid in chain_order:
        chains.append((cid, chain_counts[cid]))

    return chains


def plot_model(predictions_dir, model_num, output_path=None):
    """Generate combined PAE + pLDDT plot."""
    
    # Find files
    pae_path = find_file(predictions_dir, f"pae_*_model_{model_num}.npz")
    cif_path = find_file(predictions_dir, f"*_model_{model_num}.cif")
    plddt_path = find_file(predictions_dir, f"plddt_*_model_{model_num}.npz")
    conf_path = find_file(predictions_dir, f"confidence_*_model_{model_num}.json")

    print(f"PAE:   {pae_path}")
    print(f"CIF:   {cif_path}")
    print(f"pLDDT: {plddt_path}")
    print(f"Conf:  {conf_path}")

    if not pae_path:
        print("ERROR: No PAE file found!")
        sys.exit(1)

    # Load PAE
    pae_data = np.load(pae_path)
    pae = pae_data['pae']
    if pae.ndim == 3:
        pae = pae[0]
    total = pae.shape[0]

    # Load confidence
    conf = None
    if conf_path:
        try:
            with open(conf_path) as f:
                conf = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load confidence: {e}")

    # Load pLDDT
    plddt_arr = None
    if plddt_path:
        try:
            d = np.load(plddt_path)
            plddt_arr = d[d.files[0]]
            if plddt_arr.ndim == 2:
                plddt_arr = plddt_arr[0]
        except Exception as e:
            print(f"Warning: Could not load pLDDT: {e}")

    # Parse chains
    chains = parse_chains_from_cif(cif_path)
    print(f"Chains: {chains}")

    # Build chain boundaries
    boundaries = []
    offset = 0
    for i, (name, length) in enumerate(chains):
        boundaries.append({
            "name": name,
            "start": offset,
            "end": offset + length,
            "length": length,
            "color": PALETTE[i % len(PALETTE)]
        })
        offset += length

    # If no chains parsed, treat as single block
    if not boundaries:
        boundaries = [{"name": "?", "start": 0, "end": total, "length": total, "color": PALETTE[0]}]

    # --- Create figure ---
    fig = plt.figure(figsize=(18, 8))
    gs = GridSpec(1, 2, width_ratios=[3, 2], wspace=0.3)
    ax_pae = fig.add_subplot(gs[0])
    ax_plddt = fig.add_subplot(gs[1])

    # PAE heatmap
    im = ax_pae.imshow(pae, cmap='viridis', vmin=0, vmax=31.75, aspect='equal',
                        interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax_pae, fraction=0.046, pad=0.04)
    cbar.set_label("Predicted Aligned Error (Å)", fontsize=10)

    ax_pae.set_xlabel("Scored residue index", fontsize=11, fontweight="bold")
    ax_pae.set_ylabel("Aligned residue index", fontsize=11, fontweight="bold")
    ax_pae.set_title("PAE", fontsize=12, fontweight="bold")

    # Draw chain boundary lines
    for b in boundaries:
        for ax_line_pos in [b["start"]]:
            if ax_line_pos > 0:
                ax_pae.axhline(y=ax_line_pos - 0.5, color="white", linewidth=0.8,
                               linestyle="--", alpha=0.7)
                ax_pae.axvline(x=ax_line_pos - 0.5, color="white", linewidth=0.8,
                               linestyle="--", alpha=0.7)

    # Tick positions at chain boundaries
    if boundaries:
        tick_positions = [b["start"] for b in boundaries] + [total]
        ax_pae.set_xticks(tick_positions)
        ax_pae.set_yticks(tick_positions)
        ax_pae.set_xticklabels([str(t) for t in tick_positions], fontsize=9)
        ax_pae.set_yticklabels([str(t) for t in tick_positions], fontsize=9)

    ax_pae.set_xlim(-0.5, total - 0.5)
    ax_pae.set_ylim(total - 0.5, -0.5)

    # pLDDT plot
    if plddt_arr is not None:
        x = np.arange(len(plddt_arr))

        for b in boundaries:
            mask = (x >= b["start"]) & (x < b["end"])
            ax_plddt.fill_between(x[mask], plddt_arr[mask], alpha=0.3, color=b["color"])
            ax_plddt.plot(x[mask], plddt_arr[mask], color=b["color"], linewidth=0.8,
                          label=f"Chain {b['name']} ({b['length']})")
            if b["start"] > 0:
                ax_plddt.axvline(x=b["start"], color="gray", linewidth=0.5,
                                 linestyle="--", alpha=0.5)

        # Background shading for chain regions
        for b in boundaries:
            ax_plddt.axvspan(b["start"], b["end"], alpha=0.05, color=b["color"])

        ax_plddt.set_ylim(0, 1.05)
        ax_plddt.set_xlim(0, len(plddt_arr))
        ax_plddt.set_xlabel("Residue index", fontsize=11, fontweight="bold")
        ax_plddt.set_ylabel("pLDDT", fontsize=11, fontweight="bold")
        ax_plddt.set_title("pLDDT", fontsize=12, fontweight="bold")
        ax_plddt.axhline(y=0.7, color="red", linewidth=0.8, linestyle="--",
                          alpha=0.5, label="0.7 threshold")
        ax_plddt.legend(fontsize=7, loc="lower left")
        ax_plddt.grid(alpha=0.2)

    # Suptitle with stats
    name = os.path.basename(pae_path).replace("pae_", "").replace(".npz", "")
    # Extract variant info from directory name
    dir_basename = os.path.basename(predictions_dir.rstrip('/'))
    
    if conf:
        stats = (f"Conf: {conf.get('confidence_score', 0):.4f}  |  "
                 f"pTM: {conf.get('ptm', 0):.4f}  |  "
                 f"ipTM: {conf.get('iptm', 0):.4f}  |  "
                 f"pLDDT: {conf.get('complex_plddt', 0):.4f}")
        fig.suptitle(f"Boltz-2 — {dir_basename}\n{stats}",
                     fontsize=14, fontweight="bold", y=1.02)
    else:
        fig.suptitle(f"Boltz-2 — {name}", fontsize=14, fontweight="bold", y=1.02)

    plt.tight_layout()

    # Output path
    if not output_path:
        output_path = f"{name}_plot.png"
    
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 plot_model.py <predictions_dir> <model_number> [output.png]")
        print("\nExample:")
        print("  python3 plot_model.py 2a3d/replacement/5dk3_2a3d_replacement_run_3_outputs/"
              "boltz_results_5dk3_2a3d_replacement/predictions/5dk3_2a3d_replacement 0")
        sys.exit(1)

    pred_dir = sys.argv[1]
    model_n = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None

    plot_model(pred_dir, model_n, out)
