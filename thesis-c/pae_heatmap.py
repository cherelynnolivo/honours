# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "numpy",
#     "matplotlib",
# ]
# ///
"""
PAE Heatmap Generator for Boltz-2 outputs.

Usage:
    uv run pae_heatmap.py <path_to_output_cif>

The script auto-detects chain order and residue counts from the CIF,
finds the corresponding PAE NPZ and confidence JSON files, and generates
a labelled heatmap with confidence metrics in the title.

Example:
    uv run pae_heatmap.py base_FG_fixed/5dk3_base_FG_fixed_outputs/boltz_results_5dk3_base_FG_fixed/predictions/5dk3_base_FG_fixed/5dk3_base_FG_fixed_model_0.cif
"""

import sys
import os
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import OrderedDict


def parse_chains_from_cif(cif_path):
    """
    Parse chain IDs and residue counts from a Boltz-2 output CIF file.
    Returns an OrderedDict of {chain_id: residue_count} in file order.

    Uses _atom_site loop, counting unique (label_asym_id, label_seq_id) pairs
    for ATOM records only (skips HETATM).
    """
    chains = OrderedDict()
    seen_residues = {}  # (chain_id, seq_id) -> True

    in_atom_site = False
    field_indices = {}
    field_count = 0

    with open(cif_path, "r") as f:
        for line in f:
            line_stripped = line.strip()

            # Detect start of _atom_site loop
            if line_stripped == "loop_":
                in_atom_site = False
                field_indices = {}
                field_count = 0
                continue

            if line_stripped.startswith("_atom_site."):
                in_atom_site = True
                field_name = line_stripped.split(".")[1]
                field_indices[field_name] = field_count
                field_count += 1
                continue

            # If we were in _atom_site and hit a non-data line, we're done
            if in_atom_site and (
                line_stripped.startswith("_")
                or line_stripped.startswith("#")
                or line_stripped == ""
            ):
                if line_stripped.startswith("_") and not line_stripped.startswith(
                    "_atom_site."
                ):
                    in_atom_site = False
                if line_stripped.startswith("#") or line_stripped == "":
                    in_atom_site = False
                continue

            if not in_atom_site:
                continue

            # Parse data line
            tokens = line_stripped.split()
            if len(tokens) < field_count:
                continue

            # Get record type - skip HETATM
            group_pdb_idx = field_indices.get("group_PDB")
            if group_pdb_idx is not None:
                if tokens[group_pdb_idx] != "ATOM":
                    continue

            # Get chain and residue info
            chain_idx = field_indices.get("label_asym_id")
            seq_idx = field_indices.get("label_seq_id")

            if chain_idx is None or seq_idx is None:
                continue

            chain_id = tokens[chain_idx]
            seq_id = tokens[seq_idx]

            key = (chain_id, seq_id)
            if key not in seen_residues:
                seen_residues[key] = True
                if chain_id not in chains:
                    chains[chain_id] = 0
                chains[chain_id] += 1

    return chains


def find_pae_npz(cif_path):
    """
    Given a CIF path like .../predictions/name/name_model_0.cif,
    find the corresponding PAE NPZ file.
    """
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    cif_basename = os.path.basename(cif_path).replace(".cif", "")

    # Primary: pae_<name>_model_N.npz in same directory as CIF
    pae_npz = os.path.join(cif_dir, f"pae_{cif_basename}.npz")
    if os.path.exists(pae_npz):
        return pae_npz

    # Extract model number for fallback searches
    if "_model_" in cif_basename:
        model_num = cif_basename.split("_model_")[-1]
        name = cif_basename.rsplit("_model_", 1)[0]
    else:
        model_num = "0"
        name = cif_basename

    pred_name_dir = cif_dir
    predictions_dir = os.path.dirname(pred_name_dir)
    results_dir = os.path.dirname(predictions_dir)

    fallback_patterns = [
        os.path.join(
            results_dir, "confidence", name, f"confidence_model_{model_num}.npz"
        ),
        os.path.join(cif_dir, f"*pae*model_{model_num}.npz"),
        os.path.join(results_dir, "**", f"*pae*model_{model_num}.npz"),
    ]
    for pattern in fallback_patterns:
        matches = glob.glob(pattern, recursive=True)
        matches = [m for m in matches if not os.path.basename(m).startswith("._")]
        if matches:
            return matches[0]

    return None


def find_confidence_json(cif_path):
    """
    Given a CIF path like .../predictions/name/name_model_0.cif,
    find the corresponding confidence JSON file.

    Boltz-2 writes: confidence_<name>_model_N.json in same directory.
    """
    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    cif_basename = os.path.basename(cif_path).replace(".cif", "")

    json_path = os.path.join(cif_dir, f"confidence_{cif_basename}.json")
    if os.path.exists(json_path):
        return json_path

    # Fallback glob
    matches = glob.glob(os.path.join(cif_dir, f"confidence*{cif_basename}*.json"))
    matches = [m for m in matches if not os.path.basename(m).startswith("._")]
    if matches:
        return matches[0]

    return None


def load_confidence_metrics(json_path):
    """
    Load confidence metrics from a Boltz-2 confidence JSON file.
    Returns a dict with keys: confidence_score, ptm, iptm, avg_plddt.
    Returns None if file not found or parsing fails.
    """
    if json_path is None or not os.path.exists(json_path):
        return None
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        # Debug: print keys from first JSON so you can verify key names
        if not hasattr(load_confidence_metrics, "_printed_keys"):
            print(f"  [DEBUG] confidence JSON keys: {list(data.keys())}")
            load_confidence_metrics._printed_keys = True
        # Boltz-2 may store pLDDT under different key names
        plddt = (data.get("complex_plddt")
                 or data.get("plddt")
                 or data.get("avg_plddt")
                 or data.get("mean_plddt"))
        return {
            "confidence_score": data.get("confidence_score"),
            "ptm": data.get("ptm"),
            "iptm": data.get("iptm"),
            "avg_plddt": plddt,
        }
    except Exception as e:
        print(f"WARNING: Could not parse confidence JSON: {e}")
        return None


def plot_pae_heatmap(pae, chains, title=None, metrics=None, output_path=None, vmax=32, marks=None):
    """
    Plot a PAE heatmap with viridis colormap and chain boundary labels.

    Args:
        pae: 2D numpy array (N x N) of PAE values
        chains: OrderedDict of {chain_id: residue_count} in matrix order
        title: main plot title (line 1)
        metrics: dict with confidence_score, ptm, iptm, avg_plddt (for line 2)
        output_path: where to save the figure
        vmax: maximum value for colour scale
        marks: list of (chain_id, start_res, end_res) tuples to highlight
    """
    total_residues = sum(chains.values())

    if pae.shape[0] != total_residues:
        print(
            f"WARNING: PAE matrix size ({pae.shape[0]}) != total residues from CIF ({total_residues})"
        )
        print(f"  Chain order and sizes: {dict(chains)}")
        print(f"  This may indicate a chain ordering mismatch.")

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    im = ax.imshow(
        pae,
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        aspect="equal",
        interpolation="nearest",
    )

    # Calculate chain boundaries
    boundaries = []
    pos = 0
    for chain_id, count in chains.items():
        boundaries.append(pos)
        pos += count
    boundaries.append(pos)

    # (chain boundaries indicated by tick labels only — no overlay lines)

    # Primary ticks: chain letters (A, B, F, G) at midpoints
    midpoints = []
    mid_labels = []
    for chain_id, count in chains.items():
        start = boundaries[list(chains.keys()).index(chain_id)]
        midpoints.append(start + count / 2)
        mid_labels.append(chain_id)

    ax.set_xticks(midpoints)
    ax.set_xticklabels(mid_labels, fontsize=11, fontweight="bold")
    ax.set_yticks(midpoints)
    ax.set_yticklabels(mid_labels, fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", which="major", length=0, pad=4)

    # Minor ticks: residue indices at chain boundaries (grey, small)
    boundary_pos = [b - 0.5 for b in boundaries]
    boundary_labels = [str(b) for b in boundaries]
    ax.set_xticks(boundary_pos, minor=True)
    ax.set_xticklabels(boundary_labels, minor=True, fontsize=8, color="grey")
    ax.set_yticks(boundary_pos, minor=True)
    ax.set_yticklabels(boundary_labels, minor=True, fontsize=8, color="grey")
    ax.tick_params(axis="both", which="minor", length=3, pad=1)

    ax.set_xlabel("Scored residue index", fontsize=12)
    ax.set_ylabel("Aligned residue index", fontsize=12)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, label="Predicted Aligned Error (Å)")

    # Draw marker bands for highlighted regions
    if marks:
        # Convert chain-local residue ranges to global PAE matrix indices
        chain_starts = {}
        pos = 0
        for cid, count in chains.items():
            chain_starts[cid] = pos
            pos += count

        for chain_id, res_start, res_end in marks:
            if chain_id not in chain_starts:
                print(f"WARNING: mark chain '{chain_id}' not in structure, skipping")
                continue
            # Global indices in PAE matrix (0-based, res numbers are 1-based)
            g_start = chain_starts[chain_id] + (res_start - 1) - 0.5
            g_end   = chain_starts[chain_id] + res_end - 0.5
            g_mid   = (g_start + g_end) / 2
            width = g_end - g_start

            # Horizontal band (full width)
            ax.add_patch(plt.Rectangle(
                (-0.5, g_start), pae.shape[1], width,
                linewidth=0, edgecolor="none",
                facecolor="red", alpha=0.15,
            ))
            # Vertical band (full height)
            ax.add_patch(plt.Rectangle(
                (g_start, -0.5), width, pae.shape[0],
                linewidth=0, edgecolor="none",
                facecolor="red", alpha=0.15,
            ))
            # Thin border lines for clarity
            for edge in (g_start, g_end):
                ax.axhline(y=edge, color="red", linewidth=0.6, alpha=0.5)
                ax.axvline(x=edge, color="red", linewidth=0.6, alpha=0.5)

            # Label on axes: "C:start-end" at the midpoint of the marked region
            label = f"{chain_id}:{res_start}-{res_end}"
            ax.annotate(label, xy=(0, g_mid), xycoords=("axes fraction", "data"),
                        xytext=(-4, 0), textcoords="offset points",
                        fontsize=6, color="red", fontweight="bold",
                        ha="right", va="center")
            ax.annotate(label, xy=(g_mid, 0), xycoords=("data", "axes fraction"),
                        xytext=(0, -4), textcoords="offset points",
                        fontsize=6, color="red", fontweight="bold",
                        ha="center", va="top", rotation=90)

    # Build title
    if title:
        full_title = title
        if metrics:
            conf = metrics.get("confidence_score")
            ptm = metrics.get("ptm")
            iptm = metrics.get("iptm")
            plddt = metrics.get("avg_plddt")

            parts = []
            if conf is not None:
                parts.append(f"Conf: {conf:.4f}")
            if ptm is not None:
                parts.append(f"pTM: {ptm:.4f}")
            if iptm is not None:
                parts.append(f"ipTM: {iptm:.4f}")
            if plddt is not None:
                parts.append(f"pLDDT: {plddt:.4f}")

            if parts:
                full_title += "\n" + "  |  ".join(parts)

        ax.set_title(full_title, fontsize=13, pad=12, fontweight="bold")

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved: {output_path}")
    else:
        plt.show()

    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python pae_heatmap.py <path_to_output_cif> [--vmax 32] [--mark A:111-117 ...] [--auto-mark 20]")
        sys.exit(1)

    cif_path = sys.argv[1]
    vmax = 32
    if "--vmax" in sys.argv:
        vmax_idx = sys.argv.index("--vmax") + 1
        vmax = float(sys.argv[vmax_idx])

    # Parse --mark CHAIN:START-END flags (can appear multiple times)
    # Example: --mark A:111-117 B:111-117 F:111-117 G:111-117
    marks = []
    if "--mark" in sys.argv:
        mark_idx = sys.argv.index("--mark") + 1
        while mark_idx < len(sys.argv) and not sys.argv[mark_idx].startswith("--"):
            token = sys.argv[mark_idx]
            try:
                chain_id, rng = token.split(":")
                start, end = rng.split("-")
                marks.append((chain_id, int(start), int(end)))
            except ValueError:
                print(f"WARNING: could not parse mark '{token}', expected CHAIN:START-END")
            mark_idx += 1

    # Parse --auto-mark threshold (detect high-PAE rows per chain)
    # Example: --auto-mark 20   (marks rows with mean PAE > 20 Å)
    auto_mark_threshold = None
    if "--auto-mark" in sys.argv:
        ami = sys.argv.index("--auto-mark") + 1
        auto_mark_threshold = float(sys.argv[ami]) if ami < len(sys.argv) else 20.0

    if not os.path.exists(cif_path):
        print(f"ERROR: CIF file not found: {cif_path}")
        sys.exit(1)

    # 1. Parse chain order from CIF
    print(f"Reading CIF: {cif_path}")
    chains = parse_chains_from_cif(cif_path)
    print(f"Chain order (from CIF): {dict(chains)}")
    total = sum(chains.values())
    print(f"Total residues: {total}")

    # 2. Find and load confidence metrics
    json_path = find_confidence_json(cif_path)
    metrics = None
    if json_path:
        print(f"Confidence JSON: {json_path}")
        metrics = load_confidence_metrics(json_path)
        if metrics:
            print(f"Metrics: {metrics}")
    else:
        print("WARNING: No confidence JSON found — title will omit metrics.")

    # 3. Find PAE NPZ
    npz_path = find_pae_npz(cif_path)
    if npz_path is None:
        print("ERROR: Could not find confidence NPZ file.")
        sys.exit(1)
    print(f"PAE file:  {npz_path}")

    # 4. Load PAE
    data = np.load(npz_path)
    print(f"NPZ keys:  {data.files}")

    if "pae" in data:
        pae = data["pae"]
    elif "predicted_aligned_error" in data:
        pae = data["predicted_aligned_error"]
    else:
        print(f"ERROR: No PAE key found. Available keys: {data.files}")
        sys.exit(1)

    if pae.ndim == 3:
        print(f"PAE shape is 3D {pae.shape}, using first sample (index 0)")
        pae = pae[0]

    print(f"PAE shape: {pae.shape}")

    # 5. Reorder PAE matrix to canonical chain order: A, B, F, G
    CANONICAL_ORDER = ["A", "B", "F", "G"]
    cif_order = list(chains.keys())

    if cif_order != CANONICAL_ORDER:
        print(f"Reordering PAE from {cif_order} -> {CANONICAL_ORDER}")
        starts = {}
        pos = 0
        for chain_id in cif_order:
            starts[chain_id] = pos
            pos += chains[chain_id]

        new_indices = []
        for chain_id in CANONICAL_ORDER:
            if chain_id not in chains:
                print(f"WARNING: Chain {chain_id} not found in CIF, skipping")
                continue
            s = starts[chain_id]
            new_indices.extend(range(s, s + chains[chain_id]))

        new_indices = np.array(new_indices)
        pae = pae[np.ix_(new_indices, new_indices)]
        chains = OrderedDict((c, chains[c]) for c in CANONICAL_ORDER if c in chains)

    # 6. Auto-detect high-PAE rows per chain
    if auto_mark_threshold is not None:
        print(f"\nAuto-detecting high-PAE regions (threshold: {auto_mark_threshold} Å)...")
        row_means = np.mean(pae, axis=1)  # mean PAE for each row

        chain_start = 0
        for chain_id, count in chains.items():
            chain_rows = row_means[chain_start : chain_start + count]
            above = chain_rows > auto_mark_threshold

            # Find contiguous stretches of True
            in_region = False
            region_start = None
            for i in range(count):
                if above[i] and not in_region:
                    in_region = True
                    region_start = i
                elif not above[i] and in_region:
                    in_region = False
                    # Convert to 1-based residue numbers
                    marks.append((chain_id, region_start + 1, i))
                    print(f"  {chain_id}:{region_start + 1}-{i}  "
                          f"(mean PAE {np.mean(chain_rows[region_start:i]):.1f} Å)")
            if in_region:
                marks.append((chain_id, region_start + 1, count))
                print(f"  {chain_id}:{region_start + 1}-{count}  "
                      f"(mean PAE {np.mean(chain_rows[region_start:]):.1f} Å)")

            chain_start += count

        if not any(m for m in marks):
            print("  No regions above threshold found.")

    # 7. Build title
    cif_basename = os.path.basename(cif_path).replace(".cif", "")
    # Pretty-print the name: replace underscores with spaces for display
    display_name = cif_basename.replace("_", " ")
    title = f"Boltz-2 PAE — {display_name}"

    # 7. Output path
    output_dir = os.path.dirname(cif_path)
    output_png = os.path.join(output_dir, f"{cif_basename}_pae.png")

    # 8. Plot
    plot_pae_heatmap(
        pae, chains, title=title, metrics=metrics, output_path=output_png, vmax=vmax,
        marks=marks if marks else None,
    )

    # Summary
    chain_list = list(chains.keys())
    print(f"\nChain order used for plot: {' -> '.join(chain_list)}")
    pos = 0
    for chain_id, count in chains.items():
        print(f"  {chain_id}: residues {pos}–{pos + count - 1} ({count} residues)")
        pos += count


if __name__ == "__main__":
    main()
