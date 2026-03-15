#!/bin/bash
# generate_best_pae_controls.sh — Run from root dir containing base/, base_FG_fixed/
# Outputs go into the existing best_pae_plots/ folder alongside variant results.

mkdir -p best_pae_plots
TMPFILE=$(mktemp)

find_best_conf() {
  local pattern="$1"
  local best_score="0"
  local best_file=""
  for f in $pattern; do
    [ ! -f "$f" ] && continue
    score=$(jq -r '.confidence_score' "$f" 2>/dev/null || echo "0")
    is_better=$(python3 -c "print(1 if float('${score}') > float('${best_score}') else 0)" 2>/dev/null || echo "0")
    if [ "$is_better" -eq 1 ]; then
      best_score="$score"
      best_file="$f"
    fi
  done
  echo "$best_file"
}

# Process each control directory
for control in base base_FG_fixed; do
  [ ! -d "$control" ] && echo "Skipping $control (not found)" && continue

  # Flat structure: control/5dk3_*_outputs/boltz_results_*/predictions/*/confidence_*.json
  file=$(find_best_conf "${control}"/*/boltz_results_*/predictions/*/confidence_*.json)
  if [ -n "$file" ] && [ -f "$file" ]; then
    dir=$(dirname "$file")
    model=$(echo "$file" | sed 's/.*\(model_[0-9]*\).*/\1/')
    pae=$(find "$dir" -name "pae_*${model}.npz" ! -name "._*" | head -1)
    cif=$(find "$dir" -name "*${model}.cif" ! -name "._*" | head -1)
    plddt_f=$(find "$dir" -name "plddt_*${model}.npz" ! -name "._*" | head -1)
    [ -n "$pae" ] && echo "${control}|control|${pae}|${file}|${plddt_f}|${cif}" >> "$TMPFILE"
  else
    echo "No confidence files found for $control"
  fi
done

echo "=== Best control models found ==="
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
    variant, strategy, pae_path, conf_path, plddt_path, cif_path = entry

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

    fig.suptitle(f"Boltz-2 — {variant} {strategy}\n{stats_line}",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    outname = f"best_pae_plots/{variant}_{strategy}.png"
    plt.savefig(outname, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {outname}")

    # Collect stats
    m_model = re.search(r'model_(\d+)', conf_path)
    model = f"model_{m_model.group(1)}" if m_model else "?"

    all_stats.append({
        "variant": variant, "strategy": strategy,
        "conf": conf["confidence_score"], "ptm": conf["ptm"],
        "iptm": conf["iptm"], "plddt": conf["complex_plddt"],
        "pair_iptm": conf.get("pair_chains_iptm", {}),
        "chains": [(c["name"], c["length"]) for c in chains],
        "source": model,
        "img": f"{variant}_{strategy}.png",
    })

# === Summary markdown ===
md = []
md.append("# Boltz-2 Control Predictions — Best Models Summary")
md.append("")
md.append("Pembrolizumab (5DK3) control conditions. Best model per control from all predictions.")
md.append("")
md.append("## Overview")
md.append("")
md.append("| Control | Confidence | pTM | ipTM | pLDDT | Source |")
md.append("|---------|-----------|------|------|-------|--------|")
for d in sorted(all_stats, key=lambda x: -x["conf"]):
    md.append(f"| {d['variant']} | **{d['conf']:.4f}** | {d['ptm']:.4f} | {d['iptm']:.4f} | {d['plddt']:.4f} | {d['source']} |")

md.append("")
md.append("## Per-Control Details")
md.append("")

for d in all_stats:
    md.append(f"### {d['variant']}")
    md.append("")
    md.append(f"**Source**: {d['source']} | **Confidence**: {d['conf']:.4f} | **pTM**: {d['ptm']:.4f} | **ipTM**: {d['iptm']:.4f} | **pLDDT**: {d['plddt']:.4f}")
    md.append("")
    if d["chains"]:
        md.append("**Chains:** " + ", ".join(f"{n}={l}" for n, l in d["chains"]) + f" (total: {sum(l for _,l in d['chains'])})")
        md.append("")
    pair = d["pair_iptm"]
    if pair and d["chains"]:
        chain_names = [c[0] for c in d["chains"]]
        md.append("**Inter-chain iPTM:**")
        md.append("")
        md.append("|  | " + " | ".join(chain_names) + " |")
        md.append("|--|" + "|".join(["---"] * len(chain_names)) + "|")
        for i, (ki, row) in enumerate(sorted(pair.items(), key=lambda x: int(x[0]))):
            name = chain_names[i] if i < len(chain_names) else ki
            vals = " | ".join(f"{float(row[kj]):.3f}" for kj in sorted(row.keys(), key=int))
            md.append(f"| **{name}** | {vals} |")
        md.append("")
    md.append(f"![{d['variant']}]({d['img']})")
    md.append("")
    md.append("---")
    md.append("")

with open("best_pae_plots/controls_summary.md", "w") as f:
    f.write("\n".join(md))
print("\n  Saved: best_pae_plots/controls_summary.md")
print("\nDone!")
PYEOF

rm "$TMPFILE"