# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "matplotlib",
# ]
# ///
"""
PAE Comparison Script for Boltz-2 outputs.

Loads all models in a predictions folder, ranks by confidence score,
and produces two output files:
  1. <prefix>_all_pae_grid.png    — grid of PAE heatmaps for all models,
                                    each with residue-index axes, chain labels,
                                    and all confidence metrics in the title
  2. <prefix>_metrics_summary.png — bar plots of Conf, pTM, ipTM, avg pLDDT

Usage:
    uv run pae_compare.py <predictions_dir> [--top 10] [--vmax 32]

    --top N   Highlights top N in the metrics summary (default: 10)
    --vmax N  PAE colour scale maximum in Å (default: 32)

Example:
    uv run pae_compare.py base/5dk3_base_outputs/boltz_results_5dk3_base/predictions/5dk3_base/
    uv run pae_compare.py base/.../5dk3_base/ --top 5 --vmax 30
"""

import sys
import os
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from collections import OrderedDict


CANONICAL_ORDER = ["A", "B", "F", "G"]


# ─────────────────────────────────────────────────────
# Data loading helpers
# ─────────────────────────────────────────────────────

def parse_chains_from_cif(cif_path):
    chains = OrderedDict()
    seen_residues = {}
    in_atom_site = False
    field_indices = {}
    field_count = 0

    with open(cif_path, "r") as f:
        for line in f:
            ls = line.strip()
            if ls == "loop_":
                in_atom_site = False
                field_indices = {}
                field_count = 0
                continue
            if ls.startswith("_atom_site."):
                in_atom_site = True
                field_indices[ls.split(".")[1]] = field_count
                field_count += 1
                continue
            if in_atom_site and (ls.startswith("_") or ls.startswith("#") or ls == ""):
                if not ls.startswith("_atom_site."):
                    in_atom_site = False
                continue
            if not in_atom_site:
                continue
            tokens = ls.split()
            if len(tokens) < field_count:
                continue
            gpdb = field_indices.get("group_PDB")
            if gpdb is not None and tokens[gpdb] != "ATOM":
                continue
            ci = field_indices.get("label_asym_id")
            si = field_indices.get("label_seq_id")
            if ci is None or si is None:
                continue
            key = (tokens[ci], tokens[si])
            if key not in seen_residues:
                seen_residues[key] = True
                chains.setdefault(tokens[ci], 0)
                chains[tokens[ci]] += 1
    return chains


def find_npz(cif_path, prefix):
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    basename = os.path.basename(cif_path).replace(".cif", "")
    path = os.path.join(cif_dir, f"{prefix}_{basename}.npz")
    if os.path.exists(path):
        return path
    if "_model_" in basename:
        model_num = basename.split("_model_")[-1]
        matches = glob.glob(os.path.join(cif_dir, f"*{prefix}*model_{model_num}.npz"))
        matches = [m for m in matches if not os.path.basename(m).startswith("._")]
        if matches:
            return matches[0]
    return None


def load_pae(cif_path):
    npz_path = find_npz(cif_path, "pae")
    if npz_path is None:
        return None
    data = np.load(npz_path)
    for key in ("pae", "predicted_aligned_error"):
        if key in data:
            arr = data[key]
            return arr[0] if arr.ndim == 3 else arr
    return None


def reorder_to_canonical(pae, chains):
    """Reorder PAE matrix to canonical chain order A, B, F, G."""
    canonical = [c for c in CANONICAL_ORDER if c in chains]
    cif_order = list(chains.keys())

    if cif_order == canonical:
        return pae, OrderedDict((c, chains[c]) for c in canonical)

    starts = {}
    pos = 0
    for chain_id in cif_order:
        starts[chain_id] = pos
        pos += chains[chain_id]

    new_indices = []
    for chain_id in canonical:
        s = starts[chain_id]
        new_indices.extend(range(s, s + chains[chain_id]))
    new_indices = np.array(new_indices)

    return (
        pae[np.ix_(new_indices, new_indices)],
        OrderedDict((c, chains[c]) for c in canonical),
    )


def load_confidence_metrics(cif_path):
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    basename = os.path.basename(cif_path).replace(".cif", "")
    json_path = os.path.join(cif_dir, f"confidence_{basename}.json")
    if not os.path.exists(json_path):
        matches = glob.glob(os.path.join(cif_dir, f"confidence*{basename}*.json"))
        matches = [m for m in matches if not os.path.basename(m).startswith("._")]
        json_path = matches[0] if matches else None
    if json_path is None or not os.path.exists(json_path):
        return {}
    try:
        with open(json_path) as f:
            d = json.load(f)
        # Debug: print keys from first JSON so you can verify key names
        if not hasattr(load_confidence_metrics, "_printed_keys"):
            print(f"  [DEBUG] confidence JSON keys: {list(d.keys())}")
            load_confidence_metrics._printed_keys = True
        # Boltz-2 may store pLDDT under different key names
        plddt = (d.get("complex_plddt")
                 or d.get("plddt")
                 or d.get("avg_plddt")
                 or d.get("mean_plddt"))
        return {
            "confidence_score": d.get("confidence_score"),
            "ptm":              d.get("ptm"),
            "iptm":             d.get("iptm"),
            "avg_plddt":        plddt,
        }
    except Exception:
        return {}


# ─────────────────────────────────────────────────────
# Per-cell drawing
# ─────────────────────────────────────────────────────

def draw_pae_on_ax(ax, pae, chains, metrics, rank, model_num, vmax):
    """
    Draw a fully-labelled PAE heatmap onto ax, matching the single-plot style:
      - residue index ticks at chain boundaries (grey)
      - chain letter labels at midpoints (white, bold)
      - two-line title: rank/model on line 1, all four metrics on line 2
    """
    im = ax.imshow(pae, cmap="viridis", vmin=0, vmax=vmax,
                   aspect="equal", interpolation="nearest")

    # Build boundary positions
    boundaries, pos = [], 0
    for count in chains.values():
        boundaries.append(pos)
        pos += count
    boundaries.append(pos)  # final = total residues

    # (chain boundaries indicated by tick labels only — no overlay lines)

    # Primary ticks: chain letters (A, B, F, G) at midpoints
    midpoints  = [boundaries[i] + c / 2 for i, c in enumerate(chains.values())]
    mid_labels = list(chains.keys())
    ax.set_xticks(midpoints)
    ax.set_xticklabels(mid_labels, fontsize=7, fontweight="bold")
    ax.set_yticks(midpoints)
    ax.set_yticklabels(mid_labels, fontsize=7, fontweight="bold")
    ax.tick_params(axis="both", which="major", length=0, pad=2)

    # Minor ticks: residue indices at chain boundaries (grey, small)
    boundary_pos    = [b - 0.5 for b in boundaries]
    boundary_labels = [str(b) for b in boundaries]
    ax.set_xticks(boundary_pos, minor=True)
    ax.set_xticklabels(boundary_labels, minor=True, fontsize=4, color="grey")
    ax.set_yticks(boundary_pos, minor=True)
    ax.set_yticklabels(boundary_labels, minor=True, fontsize=4, color="grey")
    ax.tick_params(axis="both", which="minor", length=2, pad=1)

    # Title: line 1 = rank + model; line 2 = all four metrics
    conf  = metrics.get("confidence_score")
    ptm   = metrics.get("ptm")
    iptm  = metrics.get("iptm")
    plddt = metrics.get("avg_plddt")

    line1 = f"model_{model_num}"
    parts = []
    if conf  is not None: parts.append(f"Conf={conf:.4f}")
    if ptm   is not None: parts.append(f"pTM={ptm:.4f}")
    if iptm  is not None: parts.append(f"ipTM={iptm:.4f}")
    if plddt is not None: parts.append(f"pLDDT={plddt:.4f}")
    line2 = "  |  ".join(parts)

    subtitle = f"{line1}\n{line2}" if line2 else line1
    ax.set_title(subtitle, fontsize=6, pad=3)

    return im


# ─────────────────────────────────────────────────────
# Grid output  (FIXED colorbar)
# ─────────────────────────────────────────────────────

def make_grid(models, vmax, output_path, prefix):
    n = len(models)
    ncols = min(10, n)
    nrows = (n + ncols - 1) // ncols

    # Cell sizing — extra height for two-line title
    cell_w, cell_h = 3.8, 4.2
    cbar_w = 0.35                       # width reserved for colorbar
    fig_w = ncols * cell_w + cbar_w + 1.0
    fig_h = nrows * cell_h + 0.5        # minimal extra for suptitle

    # Use GridSpec: ncols for heatmaps + 1 narrow column for colorbar
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        nrows, ncols + 1, figure=fig,
        width_ratios=[1] * ncols + [0.05],
        wspace=0.35, hspace=0.55,
        top=0.95,                        # leave just enough for title
    )

    last_im = None
    for idx, m in enumerate(models):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        last_im = draw_pae_on_ax(
            ax, m["pae"], m["chains"], m["metrics"],
            rank=idx + 1, model_num=m["model_num"], vmax=vmax
        )

    # Hide unused grid cells (heatmap columns only)
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        ax.set_visible(False)

    # Colorbar in the dedicated rightmost column
    if last_im is not None:
        cax = fig.add_subplot(gs[:, ncols])
        fig.colorbar(last_im, cax=cax, label="PAE (Å)")

    fig.suptitle(
        f"{prefix}_pae",
        fontsize=13, fontweight="bold",
    )

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved grid: {output_path}")


# ─────────────────────────────────────────────────────
# Metrics summary output
# ─────────────────────────────────────────────────────

def make_metrics_summary(models, output_path, prefix, top_n):
    """2×2 bar chart panel of all model metrics, top N highlighted in orange."""
    x      = list(range(len(models)))
    labels = [str(m["model_num"]) for m in models]
    conf   = [m["metrics"].get("confidence_score") for m in models]
    ptm    = [m["metrics"].get("ptm")              for m in models]
    iptm   = [m["metrics"].get("iptm")             for m in models]
    plddt  = [m["metrics"].get("avg_plddt")        for m in models]

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    fig.suptitle(
        f"{prefix}  —  Confidence metrics across all {len(models)} models\n"
        f"(orange = top {top_n} by confidence score)",
        fontsize=12, fontweight="bold"
    )

    metric_data = [
        (axes[0, 0], conf,  "Confidence score", "steelblue"),
        (axes[0, 1], ptm,   "pTM",              "mediumseagreen"),
        (axes[1, 0], iptm,  "ipTM",             "mediumpurple"),
        (axes[1, 1], plddt, "avg pLDDT",        "tomato"),
    ]
    top_indices = set(range(top_n))

    for ax, vals, ylabel, color in metric_data:
        if all(v is None for v in vals):
            ax.set_visible(False)
            continue
        colors = ["darkorange" if i in top_indices else color for i in range(len(vals))]
        ax.bar(x, [v if v is not None else 0 for v in vals],
               color=colors, edgecolor="none", width=0.8)
        valid = [v for v in vals if v is not None]
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_ylim(0, max(valid) * 1.12)
        ax.axvline(x=top_n - 0.5, color="darkorange",
                   linewidth=1.2, linestyle="--", alpha=0.7)
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for ax in axes[1]:
        ax.set_xlabel("Model (ranked by confidence)", fontsize=10)
        step = max(1, len(models) // 10)
        ax.set_xticks(x[::step])
        ax.set_xticklabels(labels[::step], fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved metrics summary: {output_path}")


# ─────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: uv run pae_compare.py <predictions_dir> [--top 10] [--vmax 32]")
        sys.exit(1)

    pred_dir = sys.argv[1].rstrip("/")
    top_n = 10
    vmax  = 32

    if "--top" in sys.argv:
        top_n = int(sys.argv[sys.argv.index("--top") + 1])
    if "--vmax" in sys.argv:
        vmax = float(sys.argv[sys.argv.index("--vmax") + 1])

    if not os.path.isdir(pred_dir):
        print(f"ERROR: Not a directory: {pred_dir}")
        sys.exit(1)

    cif_files = sorted(
        glob.glob(os.path.join(pred_dir, "*_model_*.cif")),
        key=lambda p: int(os.path.basename(p).split("_model_")[-1].replace(".cif", ""))
    )
    cif_files = [c for c in cif_files if not os.path.basename(c).startswith("._")]

    if not cif_files:
        print(f"ERROR: No *_model_*.cif files found in {pred_dir}")
        sys.exit(1)

    print(f"Found {len(cif_files)} CIF files in {pred_dir}")
    prefix = os.path.basename(cif_files[0]).rsplit("_model_", 1)[0]

    # ── Load metadata for all models ──
    models = []
    for cif_path in cif_files:
        basename  = os.path.basename(cif_path).replace(".cif", "")
        model_num = int(basename.split("_model_")[-1])
        metrics   = load_confidence_metrics(cif_path)
        conf      = metrics.get("confidence_score")
        models.append({
            "cif_path":  cif_path,
            "model_num": model_num,
            "metrics":   metrics,
            "conf":      conf if conf is not None else -1,
            "pae":       None,
            "chains":    None,
        })

    models.sort(key=lambda m: m["conf"], reverse=True)

    print(f"\nTop {min(top_n, len(models))} models by confidence score:")
    for i, m in enumerate(models[:top_n]):
        print(f"  #{i+1}  model_{m['model_num']:>2}  conf={m['conf']:.4f}")

    # ── Load PAE for all models ──
    print(f"\nLoading PAE matrices for all {len(models)} models...")
    for m in models:
        chains = parse_chains_from_cif(m["cif_path"])
        pae    = load_pae(m["cif_path"])
        if pae is None:
            print(f"  WARNING: No PAE for model_{m['model_num']}, skipping")
            continue
        pae, chains = reorder_to_canonical(pae, chains)
        m["pae"]    = pae
        m["chains"] = chains
        print(f"  model_{m['model_num']:>2}  PAE={pae.shape}")

    all_loaded = [m for m in models if m["pae"] is not None]

    # ── Output ──
    grid_path    = os.path.join(pred_dir, f"{prefix}_all_pae_grid.png")
    summary_path = os.path.join(pred_dir, f"{prefix}_metrics_summary.png")

    print("\nGenerating PAE grid...")
    make_grid(all_loaded, vmax=vmax, output_path=grid_path, prefix=prefix)

    print("Generating metrics summary...")
    make_metrics_summary(models, output_path=summary_path,
                         prefix=prefix, top_n=top_n)

    print(f"\nDone!")
    print(f"  {grid_path}")
    print(f"  {summary_path}")


if __name__ == "__main__":
    main()
