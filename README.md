# Welcome to ENCORE
This is an end-to-end pipeline implemented using [Snakemake](https://snakemake.readthedocs.io/) to reconstruct enzyme-constrained community metabolic models and identify reporter metabolites from metagenomes.

<img width="1494" height="761" alt="image" src="https://github.com/user-attachments/assets/54e71aa0-6e2b-4577-b8e3-d38ba06d6c8b" />


## Installation

#### Requirements
- Linux / macOS
- Snakemake ≥ 7.0  
- Python ≥ 3.10  
- Conda (MetaGEM, MetaWRAP)
- GDTK-db
- Configuration file `ENCORE_config.yaml` (Should be in the same directory with `ENCORE.sh`)

#### Clone the repository
```bash
git clone https://github.com/yourusername/ENCORE.git
cd ENCORE
```

#### Make the pipeline script executable
```bash
chmod +x ENCORE.sh
```


## Usage
```bash
bash ENCORE.sh [COMMAND] [cohort_prefix] [script_dir] [threads]
```

### Required arguments
- `-o, --output-directory` : Output directory for results  

### Optional arguments
- `-s, --scripts-dir` : Directory containing all scripts (default: `Scripts/`)  
- `-t, --thread` : Number of threads to use (default: 1)
- `-c, --cohort` : Cohort/sample prefix
- `--example`: Download example sequencing data for tests and showcases
- `-h, --help` : Show help message  

-----------------


### Example runs

```bash
# Run the full pipeline and 8 threads per job
bash ENCORE.sh -o results/ -t 8 all

# Run assemly with prefix 'cohort1' and 4 threads per job
bash ENCORE.sh -o results/ -cohort cohort1 -c 4 assembly

# Run CONCOCT binning with prefix 'cohort2' and 12 threads per job
bash ENCORE.sh -o results/ -cohort cohort1 -c 4 concoct
```
-----------------
### Available functions
#### Full pipeline
- `all` : Run the entire pipeline  

#### Metagenome
- `quality` : Read quality check (FastQC)  
- `trim` : Read trimming (fastp)  
- `assembly` : Contig assembly (MEGAHIT)  
- `crossmap` : Depth file preparation for binning  
- `maxbin` : Binning with MaxBin2  
- `concoct` : Binning with CONCOCT  
- `metabat` : Binning with MetaBat2  
- `refinement` : Bin refinement  
- `reassembly` : Bin reassembly  
- `abundance` : Abundance calculation  
- `taxonomy` : Taxonomic classification  
- `daa` : Differential Abundance Analysis

#### Genome-scale Metabolic Modelling
- `extraction` : Extract genomic & proteomic sequences  
- `gem` : GEM reconstruction  
- `ecgem` : ecGEM reconstruction  

#### Reporter Metabolite
- `network` : Construct microbial community ecGEMs network  
- `reporter` : Identify reporter metabolites  

