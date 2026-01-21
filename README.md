# Honours Project: 5DK3 Antibody Linker Design

## Overview

This project explores computational approaches for designing novel linker regions in **Pembrolizumab (5DK3)**, a full-length IgG4 therapeutic antibody. The work investigates the use of deep learning-based protein design tools to generate and optimize linker sequences connecting antibody domains.

**PDB ID:** [5DK3](https://www.rcsb.org/structure/5DK3) - Crystal Structure of Pembrolizumab

## Project Structure

```
honours/
├── README.md
├── structures/              # Input PDB structures at various spacings
├── rfdiffusion/             # RFdiffusion backbone generation scripts (.sh + .pbs)
├── proteinmpnn/             # ProteinMPNN sequence design scripts
└── boltz/                   # Boltz structure prediction scripts
    ├── scratch/             # Experimental runs (15mer, 18mer, 20mer, helix designs)
    ├── base/                # Base configuration
    ├── 10A/ - 50A/          # Distance-specific validation runs
    └── prep_boltz_simple.py # Input preparation script
```

## Methodology

### 1. Structure Preparation (`structures/`)

The antibody structure was prepared by artificially spacing the Fab variable regions from the constant regions at different distances (10Å, 15Å, 20Å, 30Å, 40Å, 50Å). This creates gaps that the design tools must fill with new linker backbones.

**Key files:**
- `5dk3_original.pdb` - Original crystal structure
- `5dk3_cleaned.pdb` - Cleaned structure (chains A, B, F, G)
- `5dk3_10A.pdb` to `5dk3_50A.pdb` - Structures with artificial gaps

**Chains:**
- **A, F**: Light chains
- **B, G**: Heavy chains
- Hinge region (residues 227-238) was noted to have missing coordinates in the original structure

### 2. Backbone Generation with RFdiffusion (`rfdiffusion/`)

RFdiffusion was used to generate novel backbone conformations for linker regions. The experiments systematically tested:

- **Gap distances**: 10Å, 15Å, 20Å, 30Å, 40Å, 50Å
- **Linker lengths**: 10, 20, 30, 40, 50 amino acids

**Script naming convention:** `{distance}A_{length}AA_random.sh` (and corresponding `.pbs` for HPC submission)

Example: `20A_30AA_random.sh` generates 30-residue linkers for a 20Å gap.

**Files:**
- `*.sh` - Shell scripts for local/interactive runs
- `*.pbs` - PBS job scripts for HPC cluster submission

**Additional experiments:**
- `heavy_light_20aa_unconditional.sh` - Unconditional linker design on both heavy and light chains
- `heavy_50aa_unconditional.sh` - 50-residue unconditional design on heavy chain
- `light_20aa_unconditional.sh` - Light chain specific design

**Contig specification example:**
```
contigmap.contigs='[A1-218/0 B1-226/B233-238/10/B239-444/0 F1-218/0 G1-229/G236-238/10/G239-444/0]'
```
This specifies fixed regions with a 10-residue insertion point.

### 3. Sequence Design with ProteinMPNN (`proteinmpnn/`)

ProteinMPNN was used to design amino acid sequences for the RFdiffusion-generated backbones.

**Workflow:**
1. Parse RFdiffusion output PDBs into JSONL format
2. Assign chains B and G (heavy chains) for design
3. Specify fixed positions (keeping original residues except linker)
4. Run ProteinMPNN to generate 8 sequences per structure

**Directory structure:**
```
proteinmpnn/
├── 10A/
│   ├── 10A_10AA_script.sh
│   ├── 10A_20AA_script.sh
│   └── ...
├── 15A/
├── 20A/
├── 30A/
├── 40A/
└── 50A/
```

### 4. Structure Validation with Boltz (`boltz/`)

Boltz was used for structure prediction and validation of designed sequences from ProteinMPNN.

**Directory structure:**
```
boltz/
├── scratch/                 # Experimental/preliminary runs
│   ├── 5dk3_15mer.pbs/yaml  # 15-residue linker predictions
│   ├── 5dk3_18mer.pbs/yaml  # 18-residue linker predictions
│   ├── 5dk3_20mer.pbs/yaml  # 20-residue linker predictions
│   ├── 5dk3_charged.pbs/yaml # Charged residue experiments
│   ├── 5dk3_1al1_heavy*.pbs/yaml # 1AL1 helix-based designs
│   └── 5dk3_2a3d_heavy_light.pbs/yaml # 2A3D helix-based designs
├── base/                    # Base configuration PBS script
├── 10A/ - 50A/              # Systematic validation by distance
│   ├── prep_boltz_simple.py # Prepares Boltz input from ProteinMPNN outputs
│   └── run_boltz_*AA_run_*.pbs # PBS scripts for each linker length/run
```

**Experiments include:**
- **Linker length variations**: 15mer, 18mer, 20mer designs
- **Helix-based linkers**: Using 1AL1 and 2A3D helical motifs
- **Charged residues**: Testing charged linker sequences
- **Systematic validation**: Running Boltz on all RFdiffusion/ProteinMPNN outputs organized by distance

## Key Parameters

| Tool | Parameter | Value |
|------|-----------|-------|
| RFdiffusion | num_designs | 9 per configuration |
| ProteinMPNN | num_seq_per_target | 8 |
| ProteinMPNN | sampling_temp | 0.1 |
| Boltz | num_models | 25 per prediction |

## Utilities

- `extract.py` - Extracts sequences from PDB files in FASTA format
- `prep_boltz_simple.py` - Prepares Boltz input YAML files from ProteinMPNN-designed sequences

## Running the Scripts

### RFdiffusion (local)
```bash
cd honours/rfdiffusion
bash 20A_30AA_random.sh
```

### RFdiffusion (HPC)
```bash
cd honours/rfdiffusion
qsub 20A_30AA_random.pbs
```

### ProteinMPNN
```bash
cd honours/proteinmpnn/20A
bash 20A_30AA_script.sh
```

### Boltz (HPC)
```bash
cd honours/boltz/20A
qsub run_boltz_30AA_run_0.pbs
```

## Dependencies

- **RFdiffusion**: Deep learning model for protein backbone generation
- **ProteinMPNN**: Neural network for protein sequence design
- **Boltz**: Structure prediction for validation
- **PyMOL**: Structure visualization and preparation
- **Biopython**: Sequence extraction utilities

## Notes

- The hinge region (residues ~227-238) in 5DK3 has missing coordinates due to flexibility, requiring special handling in contig maps
- Heavy chains (B, G) were the primary targets for linker design
- Multiple gap distances and linker lengths were tested to explore the design space systematically
- PBS scripts are provided for running on HPC clusters
- Scratch directory contains exploratory experiments with alternative linker motifs

## Author

Honours research project - UNSW
