# ENCORE Pipeline - Module Reference Guide

## Overview

This directory contains all ENCORE individual modules that can be executed standalone or as part of integrated pipelines.
The examples below run from the `Example/` folder, which has its own `ENCORE.sh` and `ENCORE_config.yaml`.

## Available Datasets

The pipeline includes a toy dataset for testing and validation:
- **Input Dataset**: `Toy_Dataset` (three samples: Sample1, Sample2, Sample3, one sub-folder per sample with `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`)
- **Location**: `Example/Toy_Dataset` (the reads are not stored on GitHub because of their size, see the main README for the download link)

## Conda environments

The modules run in the `metagem` conda environment. The metaWRAP modules (`--refinement`, `--reassembly`) run in the `metawrap` environment, which `ENCORE.sh` stacks on top automatically when its path is set under `envs:` and `module_envs:` in `ENCORE_config.yaml`. Otherwise, activate it yourself with `conda activate --stack metawrap`.

## Quick Start Examples

### From raw reads to the mandatory exchange reactions

Execute all steps from quality control to the ecGEMs and their FVA results:

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --all
```

Then run `--reporter` with the sample metadata (Step 14).

### Running Joint Modules

#### Preprocessing
Executes: Quality Check → Trimming → Assembly
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --preprocessing
```

#### Downstream Analysis
Executes: Crossmapping → Binning (MaxBin2/MetaBAT2/CONCOCT) → Refinement → Reassembly → Abundance → Extraction → Taxonomy → GEMs → ecGEMs
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --downstream
```

---

## Step-by-Step Module Execution

### Stage 1: Quality Control & Assembly (Preprocessing)

#### Step 1: Quality Check
Analyze read quality metrics using FastQC.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --quality
```

#### Step 2: Trim Sequences
Filter and trim low-quality reads using fastp.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --trim
```

#### Step 3: Assembly to Contigs
Assemble trimmed reads into contigs using MEGAHIT.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --assembly
```
---

### Stage 2 (Step 4): Crossmapping for Depth and Coverage Preparation
Prepare depth files for binning by mapping all samples to the assembly.
> **Important**: Complete Steps 1-3 before running crossmapping, as it requires assembled contigs from all samples.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --crossmap
```

---

### Stage 3: Binning, Abundance and Taxonomy
Choose one or more binning algorithms. These steps can be run in parallel:
#### Step 5.1: Binning with MaxBin2
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --maxbin
```

#### Step 5.2: Binning with MetaBAT2
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --metabat
```

#### Step 5.3: Binning with CONCOCT
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --concoct
```

#### Step 6: Bin Refinement
Refine and improve bin quality using metaWRAP.
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --refinement
```

#### Step 7: Bin Reassembly
Reassemble refined bins to improve completeness.
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --reassembly
```

#### Step 8: Calculate Bin Abundance
Calculate the abundance of each reassembled bin by mapping the quality-filtered reads back to the bins (bwa and samtools).

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --abundance
```

#### Step 9: Extract Genomic and Proteomic Sequences
Extract genes and proteins from the reassembled bins.
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --extraction
```

#### Step 10: Bin Taxonomy Classification
Assign taxonomic classifications to bins using GTDB-Tk.
> **Important**: Run Step 9 first, GTDB-Tk reads the DNA sequences of the bins.
```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --taxonomy
```

---

### Stage 4: Metabolic Modeling

#### Step 11: Reconstruct Genome-Scale Metabolic Models (GEMs)
Reconstruct metabolic models from the protein sequences using CarveMe.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --gem --keep-incomplete
```

#### Step 12: Reconstruct Enzyme-Constrained GEMs (ecGEMs)
Reconstruct enzyme-constrained genome-scale metabolic models from the GEMs with GECKO 3.0, and run flux variability analysis (FVA) on their exchange reactions.

For each sample, `--ecgem` prepares one GECKO folder per bin in `output/ecGEMs/<sample>/` and runs `meta-gecko.m` in MATLAB batch mode:

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --ecgem
```

It needs the `matlab` command on PATH, GECKO 3, RAVEN and Gurobi (paths under `matlab:` in `ENCORE_config.yaml`), Docker running (DLKcat) and internet access (UniProt, KEGG and PubChem). The MATLAB output of each sample is written to `output/ecGEMs/<sample>/matlab.log`. The FVA results of each bin are written to `output/ecGEMs/<sample>/<bin>/output/<sample>.<binID>FVA.txt`.

---

### Stage 5: Microbe-Metabolite Network & Reporter Metabolites

#### Step 13: Collect the FVA Results
Copy the FVA results of all bins to `output/fva/` and combine the mandatory exchange reactions (flux range does not include zero: `UB_exc < 0` Uptake, `LB_exc > 0` Secretion) of all samples and bins into `output/fva/fva_mandatory_exchange.txt`, with `Sample` and `Bin` columns.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --fva
```

#### Step 14: DAA, Microbe-Metabolite Network and Reporter Metabolites
`--reporter` runs `Scripts/RScripts/DAA_RM.R`:
1. Differential abundance analysis (DAA) of the species with MaAsLin2, on the relative abundance of the bins (Step 8) summed per species (Step 10). Only bins with CheckM completeness > `reporterComp` (default 70) and contamination < `reporterCont` (default 10) are used (`params:` in `ENCORE_config.yaml`, CheckM from the reassembly)
2. Microbe-metabolite network from the mandatory exchange reactions (Step 13)
3. Directional reporter metabolites (PredCMB-style, with species in place of gene families, weighted by the mean relative abundance of each species over all samples), in the entire network and in the uptake subnetwork. A metabolite linked to all species with only uptake or only secretion has no null variance: its corrected Z-score is NA and its P = 1, and BH-FDR is applied to all metabolites

Required files:
- Sample metadata (`--metadata`): tab-separated, the first column is the sample ID (as in the data folder), plus the column with the two groups compared (`--fixed-effect`, default `Group`). MaAsLin2 uses this column only, without covariates or a paired design.

```bash
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --reporter \
    --metadata metadata.tsv --fixed-effect Group
```

Results in `output/reporter_metabolite/`:

| File | Content |
|------|---------|
| `daa.txt` | DAA results (`classification`, `logFC`, `PValue`, `FDR`) |
| `maaslin2/` | MaAsLin2 output |
| `taxoID.txt` | Species ID used in MaAsLin2 and its GTDB-Tk classification |
| `MMNetwork.txt` | Microbe-metabolite network |
| `reporter_metabolites.txt` | Reporter metabolites in the entire network |
| `reporter_metabolites_uptake.txt` | Reporter metabolites in the uptake subnetwork |
| `reporter_edge_scores.txt` | Signed Z-score and abundance weight of each microbe-metabolite edge |
| `reporter_microbial_background.txt` | Background species and their abundance weights |
| `species_abundance.txt` | Relative abundance of each species in each sample |
| `abundance.summary.txt`, `checkm.summary.txt` | Mapped reads and CheckM quality of all bins (inputs of `DAA_RM.R`) |

---

## Module Files Reference

| Module | File | Purpose |
|--------|------|---------|
| Quality Check | `fastqc.smk` | FastQC quality metrics |
| Trimming | `fastp.smk` | Read quality filtering |
| Assembly | `megahit.smk` | Contig assembly |
| Crossmapping | `crossmapping.smk` | Depth file preparation |
| Binning (MaxBin2) | `maxbin.smk` | MaxBin2 binning |
| Binning (CONCOCT) | `concoct.smk` | CONCOCT binning |
| Binning (MetaBAT2) | `metabat.smk` | MetaBAT2 binning |
| Refinement | `refinement.smk` | metaWRAP bin refinement |
| Reassembly | `reassembly.smk` | metaWRAP bin reassembly |
| Abundance | `AbundanceCalculation.smk` | Calculate bin abundances |
| Extraction | `extraction.smk` | Gene/protein extraction |
| Taxonomy | `taxonomy.smk` | GTDB-Tk classification |
| GEM Reconstruction | `carveme.smk` | CarveMe GEM generation |
| ecGEM Reconstruction | `GEMtoECGEM.smk` + `Scripts/MATLABScripts/meta-gecko.m` | Enzyme-constrained GEMs and FVA of their exchange reactions |
| FVA results | `ENCORE.sh` (bash) | Mandatory exchange reactions of all samples and bins |
| Reporter Metabolites | `ENCORE.sh` (bash) + `Scripts/RScripts/DAA_RM.R` | DAA, microbe-metabolite network and reporter metabolites |
