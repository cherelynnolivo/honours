# Thesis C: Current Work - Analysis & Results

## Overview

This directory contains the scripts for analyzing Boltz predictions and comparing designed linker structures. These tools were used to evaluate the quality of designs from Thesis B.

## Directory Contents

### Analysis Scripts

| Script | Description |
|--------|-------------|
| `find_best.sh` | Finds the best performing designs from Boltz outputs |
| `gen_best_pae.sh` | Generates PAE (predicted aligned error) plots for best designs |
| `gen_best_pae_control.sh` | Generates control PAE plots for comparison |
| `gen_all_pae_control.sh` | Generates PAE plots for all control structures |

### Visualization Scripts

| Script | Description |
|--------|-------------|
| `pae_compare.py` | Compares PAE matrices between designs |
| `pae_heatmap.py` | Creates heatmap visualizations of PAE data |
| `pae_pymol.py` | PyMOL script for PAE visualization |
| `pae_to_A.pml` | PyMOL script to map PAE to chain A |
| `pae_to_B.pml` | PyMOL script to map PAE to chain B |
| `plddt_overlay.py` | Creates overlay plots of pLDDT confidence scores |
| `plot_model.py` | General model visualization script |
| `residue_count.py` | Counts residues in designed linkers |

### Data Files

| File | Description |
|------|-------------|
| `5dk3_cleaned.cif` | Cleaned 5DK3 structure in mmCIF format |

## Usage Examples

### Find best designs
```bash
./find_best.sh
```

### Generate PAE plots
```bash
./gen_best_pae.sh
```

### Create pLDDT overlay
```bash
python plddt_overlay.py
```

### Compare PAE matrices
```bash
python pae_compare.py --design1 <path> --design2 <path>
```

## Output Folders (in Thesis C parent)

Note: The following folders contain Boltz prediction outputs and are NOT copied here (too large):

- `15mer/` - 15-residue linker predictions
- `18mer/` - 18-residue linker predictions
- `20mer/` - 20-residue linker predictions
- `1al1/` - 1AL1 helix-based designs
- `2a3d/` - 2A3D helix-based designs
- `charged/` - Charged residue experiments
- `base/` - Base configuration results
- `base_FG_fixed/` - FG fixed experiments
- `best_pae_plots/` - Best PAE visualizations
- `control_pae_plots/` - Control PAE visualizations

## Dependencies

- **Python 3** with packages: numpy, matplotlib, pandas
- **PyMOL** - For molecular visualization
- **Boltz** - For PAE data generation

## Notes

- These scripts analyze output from the Thesis B workflow
- PAE (Predicted Aligned Error) indicates prediction confidence
- pLDDT (predicted LDDT) indicates per-residue confidence
- Lower PAE values = higher confidence in relative positioning
- Higher pLDDT values = higher confidence in structure
