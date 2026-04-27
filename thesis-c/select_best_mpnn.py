#!/usr/bin/env python3
"""
select_best_mpnn.py

Parse ProteinMPNN output FASTA and select the best designed sample.

Ranking logic (rank-based composite, most positive = best):
  composite = -rank_score - rank_global + 0.5 * rank_recovery

  - rank_score:    rank by score ascending (lowest = rank 1)
  - rank_global:   rank by global_score ascending (lowest = rank 1)
  - rank_recovery: rank by seq_recovery descending (highest = rank 1)

Score and global_score are weighted equally and dominate; seq_recovery
contributes at half weight as a bonus for preserving WT identity.

Usage:
  Single file:   python3 select_best_mpnn.py run_0.fa
  Batch (root):  python3 select_best_mpnn.py .          # finds all .fa recursively
  Batch (dir):   python3 select_best_mpnn.py 10A/       # finds all .fa under 10A/
"""

import argparse
import re
import sys
from pathlib import Path


def parse_mpnn_fasta(fasta_path: str) -> list[dict]:
    """Parse ProteinMPNN .fa output into a list of sample dicts."""
    samples = []
    current_header = None
    current_seq_lines = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    entry = _parse_header(current_header)
                    if entry is not None:
                        entry["sequence"] = "".join(current_seq_lines)
                        samples.append(entry)
                current_header = line[1:]
                current_seq_lines = []
            else:
                current_seq_lines.append(line)

        if current_header is not None:
            entry = _parse_header(current_header)
            if entry is not None:
                entry["sequence"] = "".join(current_seq_lines)
                samples.append(entry)

    return samples


def _parse_header(header: str) -> dict | None:
    """Extract metrics from a FASTA header. Skips the run_0 reference line."""
    if header.startswith("run_0") or "sample=" not in header:
        return None

    info = {}

    m = re.search(r"sample=(\d+)", header)
    if m:
        info["sample"] = int(m.group(1))

    m = re.search(r"T=([\d.]+)", header)
    if m:
        info["temperature"] = float(m.group(1))

    m = re.search(r"(?<![a-z_])score=([\d.]+)", header)
    if m:
        info["score"] = float(m.group(1))

    m = re.search(r"global_score=([\d.]+)", header)
    if m:
        info["global_score"] = float(m.group(1))

    m = re.search(r"seq_recovery=([\d.]+)", header)
    if m:
        info["seq_recovery"] = float(m.group(1))

    if not all(k in info for k in ("score", "global_score", "seq_recovery")):
        return None

    return info


def rank_samples(samples: list[dict]) -> list[dict]:
    """Rank samples by rank-based composite (most positive = best).

    composite = -rank_score - rank_global + 0.5 * rank_recovery
    Highest composite wins.
    """
    n = len(samples)

    by_score = sorted(range(n), key=lambda i: samples[i]["score"])
    by_global = sorted(range(n), key=lambda i: samples[i]["global_score"])
    by_recovery = sorted(range(n), key=lambda i: -samples[i]["seq_recovery"])

    for rank, idx in enumerate(by_score, 1):
        samples[idx]["rank_score"] = rank
    for rank, idx in enumerate(by_global, 1):
        samples[idx]["rank_global"] = rank
    for rank, idx in enumerate(by_recovery, 1):
        samples[idx]["rank_recovery"] = rank

    for s in samples:
        s["composite"] = (
            -s["rank_score"]
            - s["rank_global"]
            + 0.5 * s["rank_recovery"]
        )

    return sorted(samples, key=lambda x: -x["composite"])


def write_best_fasta(sample: dict, output_path: str):
    """Write the best sample to a FASTA file."""
    header = (
        f">best_sample={sample['sample']}, "
        f"score={sample['score']:.4f}, "
        f"global_score={sample['global_score']:.4f}, "
        f"seq_recovery={sample['seq_recovery']:.4f}, "
        f"composite={sample['composite']:.2f}"
    )
    seq = sample["sequence"]
    with open(output_path, "w") as f:
        f.write(header + "\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def process_single(fasta_path: str, output: str | None, top_n: int):
    """Process a single FASTA file with detailed output."""
    samples = parse_mpnn_fasta(fasta_path)
    if not samples:
        print(f"  No valid samples in {fasta_path}", file=sys.stderr)
        return None

    ranked = rank_samples(samples)

    top_n = min(top_n, len(ranked))
    print(f"\n{'Rank':<5} {'Sample':<8} {'Score':<10} {'Global':<10} {'SeqRec':<10} {'rS':<4} {'rG':<4} {'rR':<4} {'Composite':<10}")
    print("-" * 65)
    for i, s in enumerate(ranked[:top_n]):
        print(
            f"{i+1:<5} {s['sample']:<8} {s['score']:<10.4f} "
            f"{s['global_score']:<10.4f} {s['seq_recovery']:<10.4f} "
            f"{s['rank_score']:<4} {s['rank_global']:<4} {s['rank_recovery']:<4} "
            f"{s['composite']:<+10.2f}"
        )

    best = ranked[0]
    out_path = output or Path(fasta_path).stem + "_best.fa"
    write_best_fasta(best, out_path)
    print(f"\n  Best: sample {best['sample']} -> {out_path}")
    return best


def process_batch(root_dir: str, output_dir: str | None):
    """Find all .fa files under root_dir, rank each, print summary table."""
    root = Path(root_dir)
    fa_files = sorted(root.rglob("*.fa"))

    if not fa_files:
        print(f"ERROR: No .fa files found under {root_dir}", file=sys.stderr)
        sys.exit(1)

    # Where to write best fastas
    out_dir = Path(output_dir) if output_dir else root / "best_samples"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    skipped = []

    for fa in fa_files:
        # Skip any *_best.fa we previously generated
        if fa.stem.endswith("_best"):
            continue

        samples = parse_mpnn_fasta(str(fa))
        if not samples:
            skipped.append(str(fa))
            continue

        ranked = rank_samples(samples)
        best = ranked[0]

        # Label from relative path
        rel = fa.relative_to(root)
        label = str(rel)

        best["source_file"] = str(fa)
        best["label"] = label
        best["n_samples"] = len(samples)

        # Write best fasta — mirror directory structure
        out_fa = out_dir / rel.parent / (fa.stem + "_best.fa")
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        write_best_fasta(best, str(out_fa))
        best["output_file"] = str(out_fa)

        results.append(best)

    if not results:
        print("ERROR: No valid samples found in any .fa file.", file=sys.stderr)
        sys.exit(1)

    # Print summary table
    print(f"\n{'File':<50} {'Samp':<5} {'Score':<9} {'Global':<9} {'SeqRec':<8} {'Comp':<8}")
    print("=" * 89)
    for r in results:
        print(
            f"{r['label']:<50} {r['sample']:<5} {r['score']:<9.4f} "
            f"{r['global_score']:<9.4f} {r['seq_recovery']:<8.4f} {r['composite']:<+8.2f}"
        )

    print(f"\n  Processed {len(results)} files ({len(skipped)} skipped)")
    print(f"  Best FASTAs written to: {out_dir}/")

    if skipped:
        print(f"  Skipped (no valid samples): {', '.join(skipped)}")


def main():
    parser = argparse.ArgumentParser(
        description="Select best ProteinMPNN sample(s) for Boltz prediction."
    )
    parser.add_argument(
        "input",
        help="Single .fa file OR directory to scan recursively for .fa files",
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output path: file (single mode) or directory (batch mode)",
    )
    parser.add_argument(
        "-n", "--top-n", type=int, default=5,
        help="Show top N samples per file in single mode (default: 5)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)

    if input_path.is_file():
        process_single(str(input_path), args.output, args.top_n)
    elif input_path.is_dir():
        process_batch(str(input_path), args.output)
    else:
        print(f"ERROR: {args.input} is not a file or directory", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()