#!/usr/bin/env python3
"""
plddt_overlay.py — Insertion-centred pLDDT overlay across all Boltz-2 variants.

Generates a figure showing mean ± std pLDDT traces for each variant × strategy,
aligned to the linker insertion point, so dips at the insertion site
are directly comparable across linker lengths / types / strategies.

Run from the root 5DK3 directory:
    uv run --with numpy --with matplotlib python3 plddt_overlay.py

Output: plddt_insertion_overlay.png, plddt_fullchain_overlay.png
"""

import sys, os, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
from collections import defaultdict

# ============================================================
# CONFIG
# ============================================================

# Variant colours — same colour for all strategies of a given variant
VARIANT_COLORS = {
    "base":    "#1a1a1a",
    "15mer":   "#06b6d4",
    "18mer":   "#22c55e",
    "20mer":   "#3b82f6",
    "charged": "#eab308",
    "1al1":    "#ec4899",
    "2a3d":    "#a855f7",
}

VARIANT_LABELS = {
    "base":    "Base (no linker)",
    "15mer":   "G₄S ×3 (15-mer)",
    "18mer":   "G₄S ×4.5 (18-mer)",
    "20mer":   "G₄S ×5 (20-mer)",
    "charged": "Charged linker",
    "1al1":    "1AL1 helix",
    "2a3d":    "2A3D helix",
}

# Line style per strategy
STRATEGY_LS = {
    "replacement":    "-",
    "insertion_left":  "--",
    "insertion_right": ":",
}

STRATEGY_LABELS = {
    "replacement":    "repl",
    "insertion_left":  "ins_L",
    "insertion_right": "ins_R",
}

STRATEGIES = ["replacement", "insertion_left", "insertion_right"]

# Chain lengths in the BASE construct
BASE_CHAIN_A_LEN = 218
BASE_CHAIN_B_LEN = 449

# 0-indexed position of the first replaced/inserted residue in the BASE sequence
INSERTION_POS_A = 107
INSERTION_POS_B = 112

# Window around insertion point to display
WINDOW = 60

# ============================================================
# HELPERS
# ============================================================

def find_plddt_files_for_strategy(variant_dir, strategy):
    """Find all pLDDT .npz files for a specific strategy within a variant."""
    if strategy is None:
        # Base case — no strategy subdirectory
        patterns = [
            os.path.join(variant_dir, "*_outputs", "boltz_results_*", "predictions", "*", "plddt_*.npz"),
        ]
    else:
        patterns = [
            os.path.join(variant_dir, strategy, "*_outputs", "boltz_results_*", "predictions", "*", "plddt_*.npz"),
            os.path.join(variant_dir, strategy, "*", "boltz_results_*", "predictions", "*", "plddt_*.npz"),
        ]
    files = []
    for p in patterns:
        found = glob.glob(p)
        found = [f for f in found if not os.path.basename(f).startswith("._")]
        files.extend(found)
    return sorted(set(files))


def parse_chains_from_cif(cif_path):
    """Parse chain names and lengths from a CIF file (CA atoms only)."""
    chains = []
    try:
        with open(cif_path, 'r') as f:
            lines = f.readlines()
        col_names = []
        in_atom_site = False
        data_start = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('_atom_site.'):
                in_atom_site = True
                col_names.append(line.strip())
            elif in_atom_site and not line.strip().startswith('_atom_site.'):
                data_start = i
                break

        atom_col = chain_col = None
        for j, name in enumerate(col_names):
            if name == '_atom_site.label_atom_id':
                atom_col = j
            if name == '_atom_site.auth_asym_id':
                chain_col = j

        if atom_col is None or chain_col is None:
            return []

        current_chain = None
        count = 0
        for line in lines[data_start:]:
            if not line.startswith('ATOM'):
                if line.startswith('#') or line.strip() == '':
                    break
                continue
            parts = line.split()
            if len(parts) <= max(atom_col, chain_col):
                continue
            if parts[atom_col] != 'CA':
                continue
            cid = parts[chain_col]
            if cid != current_chain:
                if current_chain is not None:
                    chains.append((current_chain, count))
                current_chain = cid
                count = 1
            else:
                count += 1
        if current_chain is not None:
            chains.append((current_chain, count))
    except Exception:
        pass
    return chains


def get_chain_plddt(plddt_arr, chains, target_chain):
    """Extract pLDDT values for a specific chain, given chain ordering from CIF."""
    offset = 0
    for name, length in chains:
        if name == target_chain:
            return plddt_arr[offset:offset + length]
        offset += length
    return None


def load_plddt_for_chain_and_strategy(variant_dir, strategy, target_chain):
    """Load pLDDT arrays for a specific chain and strategy."""
    plddt_files = find_plddt_files_for_strategy(variant_dir, strategy)
    if not plddt_files:
        return []

    all_plddt = []
    for pf in plddt_files:
        try:
            d = np.load(pf)
            arr = d[d.files[0]]
            if arr.ndim == 2:
                arr = arr[0]

            # Find matching CIF
            pred_dir = os.path.dirname(pf)
            model_match = re.search(r'model_(\d+)', os.path.basename(pf))
            if model_match:
                model_id = f"model_{model_match.group(1)}"
                cif_files = glob.glob(os.path.join(pred_dir, f"*{model_id}.cif"))
                cif_files = [c for c in cif_files if not os.path.basename(c).startswith("._")]
            else:
                cif_files = glob.glob(os.path.join(pred_dir, "*.cif"))
                cif_files = [c for c in cif_files if not os.path.basename(c).startswith("._")]

            if not cif_files:
                continue

            chains = parse_chains_from_cif(cif_files[0])
            if not chains:
                continue

            total_cif = sum(l for _, l in chains)
            if total_cif != len(arr):
                continue

            chain_plddt = get_chain_plddt(arr, chains, target_chain)
            if chain_plddt is not None:
                all_plddt.append(chain_plddt)
        except Exception:
            continue

    return all_plddt


# ============================================================
# MAIN
# ============================================================

def main():
    root = "."

    # results key: (variant, strategy, chain) -> {"mean", "std", "chain_len", "n"}
    results = {}

    variant_names = list(VARIANT_COLORS.keys())

    for vname in variant_names:
        vdir = os.path.join(root, vname)
        if not os.path.isdir(vdir):
            print(f"  Skipping {vname}: directory not found")
            continue

        if vname == "base":
            strategies_to_check = [None]
        else:
            strategies_to_check = STRATEGIES

        for strategy in strategies_to_check:
            strat_label = strategy if strategy else "control"
            for chain in ["A", "B"]:
                print(f"  Loading {vname}/{strat_label} chain {chain}...")
                all_plddt = load_plddt_for_chain_and_strategy(vdir, strategy, chain)
                if not all_plddt:
                    print(f"    No data found")
                    continue

                # Group by chain length (should be uniform within a strategy)
                by_len = defaultdict(list)
                for arr in all_plddt:
                    by_len[len(arr)].append(arr)

                most_common_len = max(by_len.keys(), key=lambda k: len(by_len[k]))
                arrays = by_len[most_common_len]

                stacked = np.stack(arrays)
                mean = np.mean(stacked, axis=0)
                std = np.std(stacked, axis=0)

                key = (vname, strategy, chain)
                results[key] = {
                    "mean": mean,
                    "std": std,
                    "chain_len": most_common_len,
                    "n": len(arrays),
                }
                print(f"    Got {len(arrays)} predictions, chain length {most_common_len}")

    if not results:
        print("ERROR: No data loaded.")
        sys.exit(1)

    # ============================================================
    # Compute per-chain x-axis ranges (tight fit)
    # ============================================================
    chain_xmin = {"A": 0, "B": 0}
    chain_xmax = {"A": 0, "B": 0}
    for (vname, strategy, chain), data in results.items():
        ins_pos = INSERTION_POS_A if chain == "A" else INSERTION_POS_B
        xmin = -ins_pos
        xmax = data["chain_len"] - ins_pos - 1
        if xmin < chain_xmin[chain]:
            chain_xmin[chain] = xmin
        if xmax > chain_xmax[chain]:
            chain_xmax[chain] = xmax
    for c in ["A", "B"]:
        chain_xmin[c] -= 5
        chain_xmax[c] += 5

    # ============================================================
    # Helper to make one figure per chain
    # ============================================================
    def plot_single_chain(chain, ins_pos, strategies_to_plot, title_suffix, filename):
        fig, ax = plt.subplots(figsize=(18, 6))

        for vname in variant_names:
            color = VARIANT_COLORS[vname]

            if vname == "base":
                key = (vname, None, chain)
                if key not in results:
                    continue
                data = results[key]
                x = np.arange(data["chain_len"]) - ins_pos
                label = f"{VARIANT_LABELS[vname]} (n={data['n']})"
                ax.plot(x, data["mean"], color=color, linestyle="-",
                        linewidth=2.5, label=label, zorder=5)
                ax.fill_between(x, data["mean"] - data["std"],
                                data["mean"] + data["std"],
                                alpha=0.15, color=color, zorder=1)
            else:
                for strategy in strategies_to_plot:
                    key = (vname, strategy, chain)
                    if key not in results:
                        continue
                    data = results[key]
                    x = np.arange(data["chain_len"]) - ins_pos
                    ls = STRATEGY_LS[strategy]
                    slabel = STRATEGY_LABELS[strategy]
                    label = f"{VARIANT_LABELS[vname]} — {slabel} (n={data['n']})"
                    lw = 1.5 if strategy == "replacement" else 1.3
                    ax.plot(x, data["mean"], color=color, linestyle=ls,
                            linewidth=lw, label=label, zorder=3)
                    ax.fill_between(x, data["mean"] - data["std"],
                                    data["mean"] + data["std"],
                                    alpha=0.08, color=color, zorder=1)

        ax.axvline(x=0, color="red", linewidth=1.2, linestyle="-", alpha=0.6, zorder=6)
        ax.axhline(y=0.7, color="gray", linewidth=0.8, linestyle="--", alpha=0.4)

        ax.set_xlim(chain_xmin[chain], chain_xmax[chain])
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Residue position relative to insertion site", fontsize=11)
        ax.set_ylabel("pLDDT (mean ± std)", fontsize=11)
        ax.grid(alpha=0.15)

        # Title on the left, legend on the right — both at the top
        ax.set_title(
            f"Boltz-2 pLDDT — Chain {chain} — {title_suffix}\n"
            "Pembrolizumab (5DK3) | Mean ± std across all predictions",
            fontsize=12, fontweight="bold", loc="left"
        )
        ax.legend(fontsize=6.5, loc="upper left", framealpha=0.9, ncol=1,
                  bbox_to_anchor=(1.01, 1.0))

        plt.tight_layout()
        outpath = os.path.join(root, filename)
        plt.savefig(outpath, dpi=250, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  Saved: {outpath}")

    # ============================================================
    # Generate separate figures: per chain × per strategy set
    # ============================================================

    for focus_strategy in ["insertion_left", "insertion_right"]:
        focus_label = focus_strategy.replace("_", " ").title()
        strats = ["replacement", focus_strategy]
        ls_note = f"solid: replacement | dashed: {focus_strategy}"

        for chain, ins_pos in [("A", INSERTION_POS_A), ("B", INSERTION_POS_B)]:
            plot_single_chain(
                chain, ins_pos, strats,
                f"Replacement vs {focus_label}",
                f"plddt_{focus_strategy}_chain{chain}.png"
            )

    # Full chain overlay (all strategies)
    for chain, ins_pos in [("A", INSERTION_POS_A), ("B", INSERTION_POS_B)]:
        plot_single_chain(
            chain, ins_pos, STRATEGIES,
            "All Strategies",
            f"plddt_fullchain_chain{chain}.png"
        )

    # ============================================================
    # Summary stats
    # ============================================================

    print("\n=== Summary ===")
    print(f"{'Variant':<12} {'Strategy':<18} {'Chain':>5} {'Len':>5} {'N':>5} "
          f"{'Mean pLDDT':>11} {'Min pLDDT':>10} {'@ Insertion':>12}")

    for vname in variant_names:
        if vname == "base":
            strats = [None]
        else:
            strats = STRATEGIES

        for strategy in strats:
            for chain, base_len, ins_pos in [("A", BASE_CHAIN_A_LEN, INSERTION_POS_A),
                                              ("B", BASE_CHAIN_B_LEN, INSERTION_POS_B)]:
                key = (vname, strategy, chain)
                if key not in results:
                    continue
                data = results[key]
                mean = data["mean"]
                window_start = max(0, ins_pos - 5)
                window_end = min(len(mean), ins_pos + 20)
                local_mean = np.mean(mean[window_start:window_end])
                strat_str = strategy if strategy else "control"
                print(f"{vname:<12} {strat_str:<18} {chain:>5} {data['chain_len']:>5} "
                      f"{data['n']:>5} {np.mean(mean):>11.4f} {np.min(mean):>10.4f} "
                      f"{local_mean:>12.4f}")

    print("\nDone!")


if __name__ == "__main__":
    main()
