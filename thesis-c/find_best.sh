#!/bin/bash
# find_best.sh — Run from root dir containing 15mer/, 18mer/, etc.

printf "%-12s %-18s %-10s %-8s %-8s %-8s\n" "VARIANT" "STRATEGY" "CONF" "PTM" "IPTM" "PLDDT"
printf '%0.s=' {1..80}; echo

for variant_dir in */; do
  variant=$(basename "$variant_dir")
  [ ! -d "$variant_dir" ] && continue

  # Replacement
  for f in "$variant_dir"replacement/*/boltz_results_*/predictions/*/confidence_*.json; do
    [ ! -f "$f" ] && continue
    printf "%s %s\n" "$(jq -r '[.confidence_score, .ptm, .iptm, .complex_plddt] | map(tostring) | join(" ")' "$f")" "$f"
  done | sort -rn | head -1 | while read -r conf ptm iptm plddt file; do
    [ -z "$conf" ] && continue
    run=$(echo "$file" | grep -o 'run_[0-9]*')
    model=$(echo "$file" | grep -o 'model_[0-9]*')
    printf "%-10s %-14s %-8.4f %-7.4f %-7.4f %-7.4f %s / %s\n" "$variant" "replacement" "$conf" "$ptm" "$iptm" "$plddt" "$run" "$model"
  done

  # Insertion (left + right combined)
  for f in "$variant_dir"insertion_*/*/boltz_results_*/predictions/*/confidence_*.json; do
    [ ! -f "$f" ] && continue
    printf "%s %s\n" "$(jq -r '[.confidence_score, .ptm, .iptm, .complex_plddt] | map(tostring) | join(" ")' "$f")" "$f"
  done | sort -rn | head -1 | while read -r conf ptm iptm plddt file; do
    [ -z "$conf" ] && continue
    run=$(echo "$file" | grep -o 'run_[0-9]*')
    model=$(echo "$file" | grep -o 'model_[0-9]*')
    side=$(echo "$file" | grep -o 'insertion_[a-z]*')
    printf "%-10s %-14s %-8.4f %-7.4f %-7.4f %-7.4f %s / %s / %s\n" "$variant" "insertion" "$conf" "$ptm" "$iptm" "$plddt" "$side" "$run" "$model"
  done
done
