<h1>
  <img src="https://raw.githubusercontent.com/XinmengLiao/ENCORE/main/Logo/Blue.png" width="80">
  Welcome to ENCORE
</h1>

This is an end-to-end pipeline implemented using [Snakemake](https://snakemake.readthedocs.io/) to reconstruct enzyme-constrained community metabolic models and identify reporter metabolites from metagenomes.

<img width="3480" height="2009" alt="image" src="https://github.com/user-attachments/assets/6587de75-c7b4-4bcd-b067-5300aacb765b" />


<br>
<br>

# Installation

## 📋 Requirements
- Linux or macOS
- Snakemake ≥ 7 (all modules except `--fva` and `--reporter`, which run in plain bash)
- Python ≥ 3.10  
- Conda (metaGEM, metaWRAP)
- GTDB-Tk reference data (R220)
- MATLAB (with the `matlab` command on PATH) with GECKO 3, RAVEN Toolbox and Gurobi, and Docker (DLKcat)
- R with `dplyr`, `tidyr` and `Maaslin2` (`--reporter`)
- Configuration file `ENCORE_config.yaml` (in the same directory as `ENCORE.sh`)

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

Please follow the installations in [GECKO3 tutorial](https://github.com/SysBioChalmers/GECKO/wiki). Note: RAVEN Toolbox, Gurobi Optimizer, and DLKcat are required. [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) to reconstruct GEMs and [Gurobi Optimizer](https://www.gurobi.com/downloads/) to solve flux balanced analysis (FBA) and flux variable analysis (FVA). DLKcat to predict the enzymes' Kcat values. 

⚠️ **Note:**

1. Since GECKO 3.0 is applied via **MATLAB**, a valid MATLAB license is required.
  
2. ENCORE uses its own versions of several GECKO 3.0 functions (`Scripts/InstallationScripts/`: `fuzzyKcatMatching.m`, `getECfromDatabase_kegg.m`, `getECfromDatabase_uniprot.m`, `loadConventionalGEM.m`, `makeEcModel.m`), and `getStandardKcat_noUniprot.m` for bins whose proteome is not in UniProt. `--ecgem` puts this folder on top of the MATLAB path, so the GECKO installation does not need to be changed.

3. Set the paths of RAVEN, GECKO and Gurobi under `matlab:` in `ENCORE_config.yaml` (leave them empty if the toolboxes are already on the MATLAB path).

## 🚀 Running ENCORE
#### Clone the repository
```bash
git clone https://github.com/XinmengLiao/ENCORE.git
cd ENCORE
```

#### Edit the configuration file
`ENCORE_config.yaml` sets the paths of the conda environments (`envs:`), the GTDB-Tk reference data (`path: gtdbtk`), the MATLAB toolboxes (`matlab:`), the number of threads of each tool (`cores:`) and the tool parameters (`params:`). Paths of scripts and data in the config are relative to the folder of the config file.

<br>

# Usage
```bash
bash ENCORE.sh -o <output_dir> -c <name> -d <data_folder> [OPTIONS] <modules>
```

### ✨ Required arguments
- `-o, --output-dir` : Output directory for results  
- `-c, --cohort` : Name of the run
- `-d, --data-folder` : Folder with the raw reads, one sub-folder per sample: `<sample>/<sample>_R1.fastq.gz` and `<sample>/<sample>_R2.fastq.gz`

### ✨ Optional arguments
- `-t, --cores` : Maximum number of CPU cores used by Snakemake (default: 1). The threads of each tool are set under `cores:` in `ENCORE_config.yaml`
- `-s, --scripts-dir` : Directory containing all scripts (default: `Scripts/` of this repository)  
- `-n, --dry-run` : Show what would be run without executing
- `-k, --keep-incomplete` : Keep incomplete output files of failed jobs
- `-l, --list` : List all modules
- `-h, --help` : Show help message  

### ✨ Reporter options (`--reporter`)
- `--metadata` : Sample metadata (required), tab-separated, the first column is the sample ID
- `--fixed-effect` : Metadata column with the two groups compared (default: `Group`); MaAsLin2 uses this column only, without covariates or a paired design

### ✨ Available modules

#### Workflow groups
- `--preprocessing` : quality, trim, assembly
- `--downstream` : crossmap, maxbin, concoct, metabat, refinement, reassembly, abundance, extraction, taxonomy, gem, ecgem
- `--all` : `--preprocessing`, `--downstream` and `--fva` (from raw reads to the mandatory exchange reactions of the ecGEMs)

#### Metagenome Quality Control & Assembly
- `--quality` : FastQC - Read quality check  
- `--trim` : fastp - Read trimming and quality filtering
- `--assembly` : MEGAHIT - Contig assembly from reads

#### Metagenome Coverage & Binning
- `--crossmap` : Prepare depth files for binning using BWA mapping and JGI scripts
- `--maxbin` : MaxBin2 - Binning algorithm for genome recovery
- `--concoct` : CONCOCT - Co-assembly clustering algorithm for binning
- `--metabat` : MetaBAT2 - Metagenome binning tool

#### Genome & Abundance Analysis
- `--refinement` : metaWRAP - Refine bins from multiple binning methods
- `--reassembly` : metaWRAP - Reassemble refined bins with original reads
- `--abundance` : Calculate bin abundance from read mapping (bwa, samtools)
- `--extraction` : Extract genomic (DNA) and proteomic (protein) sequences from bins
- `--taxonomy` : GTDB-Tk - Taxonomic classification of MAGs to species level (after `--extraction`)

#### Metabolic Modeling
- `--gem` : CarveMe - Reconstruction of genome-scale metabolic models (GEMs)
- `--ecgem` : GECKO 3.0 - Reconstruction of enzyme-constrained GEMs (ecGEMs) and FVA of their exchange reactions (runs MATLAB in batch mode)

#### Community Metabolic Analysis
- `--fva` : Copy the FVA results to `<output>/fva/` and combine the mandatory exchange reactions of all samples and bins into `<output>/fva/fva_mandatory_exchange.txt`
- `--reporter` : Differential abundance analysis (MaAsLin2), microbe-metabolite network and reporter metabolites (`Scripts/RScripts/DAA_RM.R`), written to `<output>/reporter_metabolite/`

### ✨ Examples

```bash
# From raw reads to the mandatory exchange reactions of the ecGEMs
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset -t 16 --all

# Run multiple modules
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset -t 4 --quality --trim

# Perform a dry run to check the workflow before execution
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --dry-run --preprocessing

# DAA, microbe-metabolite network and reporter metabolites
bash ENCORE.sh -o ./output -c test -d ./Toy_Dataset --reporter \
    --metadata metadata.tsv --fixed-effect Group
```

### 🧪 Toy dataset
The toy dataset (three samples, ~7.5 GB of paired-end reads) is too large for GitHub. Download it [here](https://zenodo.org/records/23311537).

### 📖 Detailed tutorial: 
More detailed functions and usage illustrating with the toy datasets could be found in [ENCORE_Manual](https://github.com/XinmengLiao/ENCORE/blob/main/Scripts/Main_Functions/README.md). The runtime of each step on the toy dataset is summarized in [Example/README.md](https://github.com/XinmengLiao/ENCORE/blob/main/Example/README.md).

# License
ENCORE is released under the [MIT License](LICENSE).
