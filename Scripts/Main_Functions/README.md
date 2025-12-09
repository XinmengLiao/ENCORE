# ENCORE Pipeline - Module Reference Guide

## Overview

This directory contains all ENCORE individual modules that can be executed standalone or as part of integrated pipelines.

## Available Datasets

The pipeline includes a toy dataset for testing and validation:
- **Input Dataset**: `Toy_Dataset` (containing Sample1, Sample2, Sample3)
- **Location**: Root directory of ENCORE

## Quick Start Examples

### Running the Complete Pipeline

Execute all steps from quality control to reporter metabolite analysis:

```bash
conda activate metagem
conda activate --stack metawrap

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --all
```

### Running Joint Modules

#### Preprocessing
Executes: Quality Check → Trimming → Assembly
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --preprocessing
```

#### Downstream Analysis
Executes: Crossmapping → Binning (MaxBin2/MetaBat2/CONCOCT) → Refinement → Reassembly → Extraction → Metabolic Modeling
```bash
conda activate metagem
conda activate --stack metawrap

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --downstream
```

#### Preprocessing to Downstream Analysis
```bash
conda activate metagem
conda activate --stack metawrap

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --preprocessing --crossmap --downstream
```

---

## Step-by-Step Module Execution

### Stage 1: Quality Control & Assembly (Preprocessing)

#### Step 1: Quality Check
Analyze read quality metrics using FastQC.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --quality
```

#### Step 2: Trim Sequences
Filter and trim low-quality reads using fastp.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --trim
```

#### Step 3: Assembly to Contigs
Assemble trimmed reads into contigs using MEGAHIT.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --assembly
```
---

### Stage 2 (Step 4): Crossmapping for Depth and Coverage Preparation
Prepare depth files for binning by mapping all samples to the assembly.
> **Important**: Complete Steps 1-3 before running crossmapping, as it requires assembled contigs from all samples.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --crossmap
```

---

### Stage 3: Binning and Reconstructing ecGEMs (Downstream Analysis)
Choose one or more binning algorithms. These steps can be run in parallel:
#### Step 5.1: Binning with MaxBin2
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --maxbin
```

#### Step 5.2: Binning with MetaBat2
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --metabat
```

#### Step 5.3: Binning with CONCOCT
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --concoct
```

#### Step 6: Bin Refinement
Refine and improve bin quality using metaWRAP.
```bash
conda activate metagem
conda activate --stack metawrap

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --refinement
```

#### Step 7: Bin Reassembly
Reassemble refined bins to improve completeness.
```bash
conda activate metagem
conda activate --stack metawrap

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --reassembly
```

#### Step 8: Extract Genomic and Proteomic Sequences
Extract genes and proteins from refined bins.
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --extraction
```

#### Step 9: Bin Taxonomy Classification
Assign taxonomic classifications to bins using GTDB-Tk.
```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --taxonomy
```

#### Step 10: Reconstruct Genome-Scale Metabolic Models (GEMs)
Reconstruct metabolic models from genomic sequences using CarveMe.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --gem
```

#### Step 11: Reconstruct Ecological GEMs (ecGEMs)
Reconstruct enzyme-constrained genome-scale metabolic models from conventional models.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --ecgem
```

---

### Stage 4: Network Analysis & Reporter Metabolites

#### Step 12: Create Microbe-Metabolite Network
Construct microbial community ecological enzyme-constrained genome-scale metabolic models networks.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --network
```

#### Step 13: Identify Reporter Metabolites
Identify reporter metabolites from the community network.

```bash
conda activate metagem

bash ENCORE.sh -o ./output -c test -d Toy_Dataset --reporter
```

---

## Module Files Reference

| Module | File | Purpose |
|--------|------|---------|
| Quality Check | `fastqc.smk` | FastQC quality metrics |
| Trimming | `fastp.smk` | Read quality filtering |
| Assembly | `megahit.smk` | Contig assembly |
| **Preprocessing** | `preprocessing.smk` | Integrated: Quality → Trim → Assembly |
| Crossmapping | `crossmapping.smk` | Depth file preparation |
| Binning (MaxBin2) | `maxbin.smk` | MaxBin2 binning |
| Binning (CONCOCT) | `concoct.smk` | CONCOCT binning |
| Binning (MetaBat2) | `metabat.smk` | MetaBat2 binning |
| Refinement | `refinement.smk` | metaWRAP bin refinement |
| Reassembly | `reassembly.smk` | metaWRAP bin reassembly |
| Abundance | `AbundanceCalculation.smk` | Calculate bin abundances |
| Taxonomy | `taxonomy.smk` | GTDB-Tk classification |
| Extraction | `extraction.smk` | Gene/protein extraction |
| GEM Reconstruction | `carveme.smk` | CarveMe GEM generation |
| ecGEM Reconstruction | `GEMtoECGEM.smk` | Enzyme-constrained GEM generation |
| **Downstream** | `downstream.smk` | Integrated: Crossmap → Binning → Refinement → Extraction → Metabolic |
| Network Analysis | `CoGEMNetwork.smk` | Community network construction |
| Reporter Metabolites | `ReporterMetabolite.smk` | Reporter metabolite identification |
