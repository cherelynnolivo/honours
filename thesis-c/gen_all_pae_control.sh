#!/bin/bash
# generate_all_control_pae.sh — Run from root dir containing base/, base_FG_fixed/
# Generates PAE+pLDDT plots for EVERY model in both control conditions.
# Outputs go into control_pae_plots/

mkdir -p control_pae_plots
TMPFILE=$(mktemp)

for control in base base_FG_fixed; do
  [ ! -d "$control" ] && echo "Skipping $control (not found)" && continue

  # Find all confidence JSON files (one per model)
  for conf_file in ${control}/*/boltz_results_*/predictions/*/confidence_*.json; do
    [ ! -f "$conf_file" ] && continue

    dir=$(dirname "$conf_file")
    model=$(echo "$conf_file" | sed 's/.*\(model_[0-9]*\).*/\1/')
    score=$(jq -r '.confidence_score' "$conf_file" 2>/dev/null || echo "0")

    pae=$(find "$dir" -name "pae_*${model}.npz" ! -name "._*" | head -1)
    cif=$(find "$dir" -name "*${model}.cif" ! -name "._*" | head -1)
    plddt_f=$(find "$dir" -name "plddt_*${model}.npz" ! -name "._*" | head -1)

    if [ -n "$pae" ]; then
      echo "${control}|${model}|${pae}|${conf_file}|${plddt_f}|${cif}|${score}" >> "$TMPFILE"
    else
      echo "  WARNING: No PAE file for $control $model"
    fi
  done
done

total=$(wc -l < "$TMPFILE")
echo "=== Found $total models across controls ==="
cat "$TMPFILE"
echo ""
echo "Generating plots..."

uv run --with numpy --with matplotlib python3 - "$TMPFILE" << 'PYEOF'
import sys, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

tmpfile = sys.argv[1]
palette = ["#06b6d4", "#22c55e", "#eab308", "#ec4899", "#a855f7", "#f97316"]

with open(tmpfile) as f:
    entries = [line.strip().split('|') for line in f if line.strip()]

all_stats = []

def parse_chains(cif_path):
    chains = []
    try:
        with open(cif_path, 'r') as cf:
            lines = cf.readlines()
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
            if name == '_atom_site.label_atom_id': atom_col = j
            if name == '_atom_site.auth_asym_id': chain_col = j
        if atom_col is None or chain_col is None:
            return []
        current_chain = None
        count = 0
        for line in lines[data_start:]:
            if not line.startswith('ATOM'):
                if line.startswith('#') or line.strip() == '': break
                continue
            parts = line.split()
            if len(parts) <= max(atom_col, chain_col): continue
            if parts[atom_col] != 'CA': continue
            cid = parts[chain_col]
            if cid != current_chain:
                if current_chain is not None:
                    chains.append({"name": current_chain, "length": count})
                current_chain = cid
                count = 1
            else:
                count += 1
        if current_chain is not None:
            chains.append({"name": current_chain, "length": count})
    except:
        pass
    return chains

for entry in entries:
    variant, model, pae_path, conf_path, plddt_path, cif_path, _ = entry

    with open(conf_path) as cf:
        conf = json.load(cf)

    pae = np.load(pae_path)['pae']
    if pae.ndim == 3: pae = pae[0]
    total = pae.shape[0]

    # Load pLDDT
    plddt_arr = None
    if plddt_path:
        try:
            d = np.load(plddt_path)
            plddt_arr = d[d.files[0]]
            if plddt_arr.ndim == 2: plddt_arr = plddt_arr[0]
        except:
            pass

    # Parse chains
    chains = parse_chains(cif_path) if cif_path else []
    if chains and sum(c["length"] for c in chains) != (len(plddt_arr) if plddt_arr is not None else total):
        chains = []

    stats_line = (f"Conf: {conf['confidence_score']:.4f}  |  "
                  f"pTM: {conf['ptm']:.4f}  |  "
                  f"ipTM: {conf['iptm']:.4f}  |  "
                  f"pLDDT: {conf['complex_plddt']:.4f}")

    # === Combined figure: PAE left, pLDDT right ===
    if plddt_arr is not None:
        fig, (ax_pae, ax_plddt) = plt.subplots(1, 2, figsize=(18, 8),
            gridspec_kw={'width_ratios': [3, 1.2]})
    else:
        fig, ax_pae = plt.subplots(figsize=(10, 8))
        ax_plddt = None

    # PAE
    im = ax_pae.imshow(pae, cmap="viridis", aspect="equal", vmin=0, vmax=31,
                       interpolation="nearest", origin="upper")
    cbar = plt.colorbar(im, ax=ax_pae, shrink=0.82, pad=0.02)
    cbar.set_label("Predicted Aligned Error (Å)", fontsize=11)
    ax_pae.set_xlabel("Scored residue index", fontsize=11)
    ax_pae.set_ylabel("Aligned residue index", fontsize=11)
    ax_pae.set_title("PAE", fontsize=12, fontweight="bold")

    # pLDDT
    if ax_plddt is not None and plddt_arr is not None:
        x = np.arange(len(plddt_arr))
        if chains:
            offset = 0
            for ci, c in enumerate(chains):
                color = palette[ci % len(palette)]
                mask = (x >= offset) & (x < offset + c["length"])
                ax_plddt.fill_between(x[mask], plddt_arr[mask], alpha=0.3, color=color)
                ax_plddt.plot(x[mask], plddt_arr[mask], color=color, linewidth=0.8,
                              label=f"Chain {c['name']} ({c['length']})")
                ax_plddt.axvline(x=offset, color="gray", linewidth=0.5, linestyle="--", alpha=0.4)
                offset += c["length"]
        else:
            ax_plddt.fill_between(x, plddt_arr, alpha=0.3, color="#06b6d4")
            ax_plddt.plot(x, plddt_arr, color="#06b6d4", linewidth=0.8)

        ax_plddt.set_ylim(0, 1)
        ax_plddt.set_xlim(0, len(plddt_arr))
        ax_plddt.axhline(y=0.7, color="red", linewidth=0.8, linestyle="--", alpha=0.5, label="0.7 threshold")
        ax_plddt.set_xlabel("Residue index", fontsize=11)
        ax_plddt.set_ylabel("pLDDT", fontsize=11)
        ax_plddt.set_title("pLDDT", fontsize=12, fontweight="bold")
        ax_plddt.legend(fontsize=7, loc="lower left")
        ax_plddt.grid(alpha=0.2)

    fig.suptitle(f"Boltz-2 — {variant} {model}\n{stats_line}",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    outname = f"control_pae_plots/{variant}_{model}.png"
    plt.savefig(outname, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {outname}")

    # Extract model number for sorting
    m_num = re.search(r'model_(\d+)', model)
    model_num = int(m_num.group(1)) if m_num else 0

    all_stats.append({
        "variant": variant, "model": model, "model_num": model_num,
        "conf": conf["confidence_score"], "ptm": conf["ptm"],
        "iptm": conf["iptm"], "plddt": conf["complex_plddt"],
        "pair_iptm": conf.get("pair_chains_iptm", {}),
        "chains": [(c["name"], c["length"]) for c in chains],
        "img": f"{variant}_{model}.png",
    })

# === Summary markdown ===
md = []
md.append("# Boltz-2 Control Predictions — All Models")
md.append("")
md.append("Pembrolizumab (5DK3) control conditions: `base` (no constraints) and `base_FG_fixed` (chains F/G as template anchors).")
md.append("")

for control_name in ["base", "base_FG_fixed"]:
    ctrl_stats = sorted([d for d in all_stats if d["variant"] == control_name],
                        key=lambda x: x["model_num"])
    if not ctrl_stats:
        continue

    md.append(f"## {control_name}")
    md.append("")

    # Summary stats
    confs = [d["conf"] for d in ctrl_stats]
    ptms = [d["ptm"] for d in ctrl_stats]
    iptms = [d["iptm"] for d in ctrl_stats]
    plddts = [d["plddt"] for d in ctrl_stats]
    best = max(ctrl_stats, key=lambda x: x["conf"])

    md.append(f"**{len(ctrl_stats)} models** | "
              f"Best: {best['model']} (conf={best['conf']:.4f}) | "
              f"Mean conf: {sum(confs)/len(confs):.4f} | "
              f"Range: [{min(confs):.4f}, {max(confs):.4f}]")
    md.append("")

    # Table of all models
    md.append("| Model | Confidence | pTM | ipTM | pLDDT |")
    md.append("|-------|-----------|------|------|-------|")
    for d in ctrl_stats:
        bold = "**" if d["model"] == best["model"] else ""
        md.append(f"| {d['model']} | {bold}{d['conf']:.4f}{bold} | {d['ptm']:.4f} | {d['iptm']:.4f} | {d['plddt']:.4f} |")
    md.append("")

    # Chain info (from first model — should be same for all)
    if ctrl_stats[0]["chains"]:
        md.append("**Chains:** " + ", ".join(f"{n}={l}" for n, l in ctrl_stats[0]["chains"])
                  + f" (total: {sum(l for _,l in ctrl_stats[0]['chains'])})")
        md.append("")

    # Inter-chain iPTM for best model
    pair = best["pair_iptm"]
    if pair and best["chains"]:
        chain_names = [c[0] for c in best["chains"]]
        md.append(f"**Inter-chain iPTM (best model: {best['model']}):**")
        md.append("")
        md.append("|  | " + " | ".join(chain_names) + " |")
        md.append("|--|" + "|".join(["---"] * len(chain_names)) + "|")
        for i, (ki, row) in enumerate(sorted(pair.items(), key=lambda x: int(x[0]))):
            name = chain_names[i] if i < len(chain_names) else ki
            vals = " | ".join(f"{float(row[kj]):.3f}" for kj in sorted(row.keys(), key=int))
            md.append(f"| **{name}** | {vals} |")
        md.append("")

    # Individual plots
    md.append("### Individual Model Plots")
    md.append("")
    for d in ctrl_stats:
        md.append(f"#### {d['model']} (conf={d['conf']:.4f})")
        md.append("")
        md.append(f"![{d['variant']} {d['model']}]({d['img']})")
        md.append("")

    md.append("---")
    md.append("")

with open("control_pae_plots/summary.md", "w") as f:
    f.write("\n".join(md))
print("\n  Saved: control_pae_plots/summary.md")
print("\nDone!")
PYEOF

rm "$TMPFILE"