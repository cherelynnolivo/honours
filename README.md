# Honours Project: 5DK3 Antibody Linker Design

## Overview

This project explores computational approaches for designing novel linker regions in **Pembrolizumab (5DK3)**, a full-length IgG4 therapeutic antibody. The work investigates the use of deep learning-based protein design tools to generate and optimize linker sequences connecting antibody domains.

**PDB ID:** [5DK3](https://www.rcsb.org/structure/5DK3) - Crystal Structure of Pembrolizumab

## Project Structure

```
honours/
├── README.md                 # This file
├── thesis-b/                # Systematic spacing experiments (completed)
│   ├── structures/          # Input PDB structures at various spacings
│   ├── rfdiffusion/        # RFdiffusion backbone generation scripts
│   ├── proteinmpnn/        # ProteinMPNN sequence design
│   └── boltz/              # Boltz structure prediction/validation
└── thesis-c                # Current/ongoing work
    └── (experiments go here)
```

## Thesis Directories

### Thesis B - Systematic Experiments

Contains the systematic exploration of linker design at various spacing distances (10Å-50Å) with different linker lengths (10-50 residues).

**See:** [thesis-b/README.md](thesis-b/README.md) for detailed documentation

- **structures/** - 5DK3 PDB files at various spacings
- **rfdiffusion/** - Backbone generation scripts (distance × length matrix)
- **proteinmpnn/** - Sequence design organized by distance
- **boltz/** - Structure validation (scratch experiments + systematic runs)

### Thesis C - Current Work

Active experiments and new directions. See [thesis-c/README.md](thesis-c/README.md) for details.

## Methodology Summary

1. **Structure Preparation** - Artificially space Fab regions from constant regions (10Å-50Å gaps)
2. **RFdiffusion** - Generate novel backbone conformations for linker regions
3. **ProteinMPNN** - Design amino acid sequences for generated backbones
4. **Boltz** - Validate designed sequences through structure prediction

## Key Parameters

| Tool | Parameter | Value |
|------|-----------|-------|
| RFdiffusion | num_designs | 9-10 per configuration |
| ProteinMPNN | num_seq_per_target | 8 |
| ProteinMPNN | sampling_temp | 0.1 |
| Boltz | num_models | 25 per prediction |

## Dependencies

- **RFdiffusion**: Deep learning model for protein backbone generation
- **ProteinMPNN**: Neural network for protein sequence design
- **Boltz**: Structure prediction for validation
- **PyMOL**: Structure visualization and preparation
- **Biopython**: Sequence extraction utilities

## Author

Honours research project - UNSW
