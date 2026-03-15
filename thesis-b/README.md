# Thesis B: 5DK3 Antibody Linker Design - Systematic Spacing Experiments

## Overview

This directory contains the systematic exploration of linker design for the 5DK3 antibody (Pembrolizumab) at various spacing distances. The work tests different gap distances (10Å-50Å) and linker lengths (10-50 residues) to find optimal linker designs.

## Directory Structure

```
thesis-b/
├── structures/           # Input PDB structures at various spacings
├── rfdiffusion/         # RFdiffusion backbone generation scripts
├── proteinmpnn/         # ProteinMPNN sequence design (by distance)
└── boltz/              # Boltz structure validation
    ├── scratch/         # Experimental runs (15mer, 18mer, 20mer, helix designs)
    ├── base/            # Base configuration
    └── 10A/ - 50A/      # Distance-specific validation runs
```

## Structures (`structures/`)

Input structures for RFdiffusion experiments, created by artificially spacing Fab regions from constant regions.

| File | Description |
|------|-------------|
| `5dk3_original.pdb` | Original crystal structure from PDB |
| `5dk3_cleaned.pdb` | Cleaned structure (chains A, B, F, G) |
| `5dk3_base.cif` | Base structure in mmCIF format |
| `5dk3_10A.pdb` - `5dk3_50A.pdb` | Structures with 10-50Å gaps |
| `5dk3_chains_AF.yaml` | YAML config for chains A & F |
| `5dk3_original.yaml` | YAML config for original structure |

### Chain Mapping
- **A, F**: Light chains
- **B, G**: Heavy chains (linker design targets)

## RFdiffusion (`rfdiffusion/`)

Scripts for generating novel backbone conformations for linker regions.

### Script Naming Convention
`{distance}A_{length}AA_random.{sh|pbs}`

- Distance: 10A, 15A, 20A, 30A, 40A, 50A
- Length: 10AA, 20AA, 30AA, 40AA, 50AA

Example: `20A_30AA_random.sh` generates 30-residue linkers for a 20Å gap.

### Key Scripts

| Script | Description |
|--------|-------------|
| `both_fab_20aa_unconditional.sh` | Unconditional design on both Fab regions |
| `heavy_light_20aa_unconditional.sh` | Heavy + light chain linker design |
| `heavy_50aa_unconditional.sh` | 50-residue unconditional on heavy chain |
| `heavy_5aa_unconditional.sh` | 5-residue short linker |
| `light_20aa_unconditional.sh` | Light chain specific design |
| `light_1AL1_helix.sh` | Helix-based linker design using 1AL1 motif |
| `extract.py` | Extracts sequences from PDB to FASTA |
| `playground.sh` | Experimental/test scripts |

### Contig Specification Example
```python
contigmap.contigs='[A1-218/0 B1-226/B233-238/30/B239-444/0 F1-218/0 G1-229/G236-238/30/G239-444/0]'
```
This specifies fixed regions (1-218, etc.) with a 30-residue insertion point.

### Parameters
- `num_designs`: 9-10 designs per configuration
- Output: `{distance}A_{length}AA_outputs/`

## ProteinMPNN (`proteinmpnn/`)

Sequence design for RFdiffusion-generated backbones, organized by distance.

```
proteinmpnn/
├── 10A/
│   ├── 10A_10AA_script.sh
│   ├── 10A_10AA_outputs/
│   └── ...
├── 15A/
├── 20A/
├── 30A/
├── 40A/
└── 50A/
```

### Workflow
1. Parse RFdiffusion output PDBs into JSONL format
2. Assign chains B and G (heavy chains) for design
3. Specify fixed positions (keep original residues, design linker)
4. Generate 8 sequences per structure

### Parameters
- `num_seq_per_target`: 8
- `sampling_temp`: 0.1

## Boltz (`boltz/`)

Structure prediction for validating designed sequences.

```
boltz/
├── scratch/           # Experimental/preliminary runs
│   ├── 5dk3_15mer.pbs/yaml   # 15-residue linker predictions
│   ├── 5dk3_18mer.pbs/yaml   # 18-residue linker predictions
│   ├── 5dk3_20mer.pbs/yaml   # 20-residue linker predictions
│   ├── 5dk3_charged.pbs/yaml # Charged residue experiments
│   ├── 5dk3_1al1_heavy*.pbs/yaml # 1AL1 helix-based designs
│   └── 5dk3_2a3d_heavy_light.pbs/yaml # 2A3D helix-based designs
├── base/              # Base configuration
├── 10A/ - 50A/       # Distance-specific validation runs
```

### Experiments
- **Linker lengths**: 15mer, 18mer, 20mer designs
- **Helix-based linkers**: Using 1AL1 and 2A3D helical motifs
- **Charged residues**: Testing charged linker sequences
- **Systematic validation**: All RFdiffusion/ProteinMPNN outputs by distance

### Parameters
- `num_models`: 25 per prediction

## Running the Workflow

### 1. RFdiffusion
```bash
cd thesis-b/rfdiffusion
bash 20A_30AA_random.sh
```

### 2. ProteinMPNN
```bash
cd thesis-b/proteinmpnn/20A
bash 20A_30AA_script.sh
```

### 3. Boltz (HPC)
```bash
cd thesis-b/boltz/20A
qsub run_boltz_30AA_run_0.pbs
```

## Dependencies

- **RFdiffusion**: Protein backbone generation
- **ProteinMPNN**: Sequence design
- **Boltz**: Structure prediction
- **PyMOL**: Structure visualization
- **Biopython**: Sequence extraction

## Notes

- Hinge region (residues ~227-238) has missing coordinates due to flexibility
- Heavy chains (B, G) were primary targets for linker design
- Systematic exploration covered all distance × length combinations
