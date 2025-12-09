<h1>
  <img src="https://raw.githubusercontent.com/XinmengLiao/ENCORE/main/Logo/Blue.png" width="80">
  Welcome to ENCORE
</h1>

This is an end-to-end pipeline implemented using [Snakemake](https://snakemake.readthedocs.io/) to reconstruct enzyme-constrained community metabolic models and identify reporter metabolites from metagenomes.

<img width="1501" height="767" alt="image" src="https://github.com/user-attachments/assets/6179646a-712e-4f2b-82b9-87897ad3fe84" />

<br>
<br>

# Installation

## 📋 Requirements
- Linux
- Snakemake ≥ 3.5  
- Python ≥ 3.10  
- Conda (MetaGEM, MetaWRAP)
- GDTK-db
- Configuration file `ENCORE_config.yaml` (Should be in the same directory with `ENCORE.sh`)

## ⚙️ Environment and tools
#### MetaGEM and MetaWRAP

ENCORE adapted functions from metaGEM (includes metaWRAP). We recommend users install metaGEM via `conda`

```bash
conda create -n metagem -c bioconda metagem
conda create -n metawrap -c bioconda metawrap
```

#### CheckM

CheckM is used to evaluate the binning quality. Although it is embeded in the `metawrap` environment, extra data needed to be downloaded:

- https://data.ace.uq.edu.au/public/CheckM_databases
- [CheckM v1 reference data](https://zenodo.org/record/7401545#.Y44ymHbMJD8)

CheckM path should be added to the environment variable: `export CHECKM_DATA_PATH=/path/to/my_checkm_data`

#### GTDB-tk

GTDB-tk is used for taxonomic classifications of the metagenome assembled genomes. ENCORE analysis is based on [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/index.html)(≥2.4) R220. The database need ~ 100G storage space. GTDB-tk can be installed in the `metagem` conda environment.

```bash
conda activate metagem
conda install bioconda::gtdbtk=2.4

# GTDB-tk reference datasets need to be downloaded manually and its path should be added to the environment variable
export GTDBTK_DATA_PATH=/PATH/TO/THE/RELEASE/
```

GTDB-tk reference datasets could be donwloaded [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing).

More GTDB-tk databases and versions could be found [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#installing).

#### GECKO 3.0

[GECKO 3.0](https://github.com/SysBioChalmers/GECKO) is a MATLAB toolbox for reconstructing enzyme-constrained genome-scale metabolic models (ecModels). It extends genome-scale metabolic models (GEMs) by including enzymatic constraints, enabling more accurate simulations of metabolism.

Please follow the installations in [GECKO3 tutorial](https://github.com/SysBioChalmers/GECKO/wiki). Note: RAVEN Toolbox and Gurobi Optimizer are required. [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) to reconstruct GEMs and [Gurobi Optimizer](https://www.gurobi.com/downloads/) to solve flux balanced analysis (FBA) and flux variable analysis (FVA).

⚠️ **Note:**

1. Since GECKO 3.0 is applied via **MATLAB**, a valid MATLAB license is required.
  
2. Several functions of GECKO 3.0 need to be replaced by the following scripts before running GECKO 3.0.
  

```bash
gecko_dir=$(find / -type d -name "GECKO\ Tool" 2>/dev/null | head -n 1)
cp Installation/GECKO3_adapted_functions/fuzzyKcatMatching.m ${gecko_dir}/src/geckomat/gather_kcats/fuzzyKcatMatching.m
cp Installation/GECKO3_adapted_functions/getECfromDatabase_kegg.m ${gecko_dir}/src/geckomat/get_enzyme_data/getECfromDatabase_kegg.m
cp Installation/GECKO3_adapted_functions/getECfromDatabase_uniprot.m ${gecko_dir}/src/geckomat/get_enzyme_data/getECfromDatabase_uniprot.m
cp Installation/GECKO3_adapted_functions/loadConventionalGEM.m ${gecko_dir}/src/geckomat/utilities/loadConventionalGEM.m
cp Installation/GECKO3_adapted_functions/makeEcModel.m ${gecko_dir}/src/geckomat/change_model/makeEcModel.m
```

## 🚀 Running ENCORE
#### Clone the repository
```bash
git clone https://github.com/yourusername/ENCORE.git
cd ENCORE
```

#### Make the pipeline script executable
```bash
chmod +x ENCORE.sh
```
<br>

# Usage
```bash
bash ENCORE.sh [OPTIONS]
```

### ✨ Required arguments
- `-o, --output-dir` : Output directory for results  
- `-c, --cohort` : Cohort/sample prefix
- `-d, --data-folder` : Input data folder containing raw sequencing reads (default: `Toy_Dataset`)

### ✨ Optional arguments
- `-s, --scripts-dir` : Directory containing all scripts (default: `Scripts/`)  
- `-t, --threads` : Number of CPU cores to use (default: 1)
- `--dry-run` : Perform a dry run without executing workflows
- `-q, --quiet` : Run in quiet mode with minimal output
- `-h, --help` : Show help message  

### ✨ Available modules

#### Full pipeline
- `--metagenome` : Run all metagenome analysis modules (quality, trim, assembly, crossmap, maxbin, concoct, metabat, refinement, reassembly, abundance, taxonomy)
- `--metabolic` : Run all metabolic modeling modules (extraction, gem, ecgem)
- `--reporter` : Run reporter metabolite analysis (network, reporter)

#### Metagenome Quality Control & Assembly
- `--quality` : FastQC - Read quality check  
- `--trim` : fastp - Read trimming and quality filtering
- `--assembly` : MEGAHIT - Contig assembly from reads

#### Metagenome Coverage & Binning
- `--crossmap` : Prepare depth files for binning using BWA mapping and JGI scripts
- `--maxbin` : MaxBin2 - Binning algorithm for genome recovery
- `--concoct` : CONCOCT - Co-assembly clustering algorithm for binning
- `--metabat` : MetaBat2 - Metagenome binning tool

#### Genome & Abundance Analysis
- `--refinement` : metaWRAP - Refine bins from multiple binning methods
- `--reassembly` : metaWRAP - Reassemble refined bins with original reads
- `--extraction` : Extract genomic (DNA) and proteomic (protein) sequences from bins
- `--abundance` : Calculate bin abundance fractions from read mapping

#### Taxonomic Classification & Modeling
- `--taxonomy` : GTDB-Tk - Taxonomic classification of MAGs to species level
- `--gem` : CarveMe - Reconstruction of genome-scale metabolic models (GEMs)
- `--ecgem` : GECKO 3.0 - Reconstruction of enzyme-constrained GEMs (ecGEMs)

#### Community Metabolic Analysis
- `--network` : Construct microbial community ecGEM metabolic networks
- `--reporter` : Identify reporter metabolites in community models

### ✨ Examples

```bash
# Run the single module with 8 jobs at the same time
bash ENCORE.sh -o ./results -c test -d my_reads -t 8 --metagenome

# Run multiple modules
bash ENCORE.sh -o ./results -c test -d my_reads -t 4 --quality --trim

# Perform a dry run to check the workflow before execution
bash ENCORE.sh -o ./results -c test -d my_reads --dry-run --metagenome

# Run the complete analysis pipeline with all modules
bash ENCORE.sh -o ./results -c test -d my_reads -t 16 --metagenome --metabolic --reporter
```
