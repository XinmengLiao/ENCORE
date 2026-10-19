# ENCORE runtime on the toy dataset

This file only summarizes the specific steps, commands, and running time of the major part of ENCORE, from the raw metagenomic sequencing reads (quality control and assembly) to the ecGEM reconstruction (Step1 to Step11). Each step was run on the three toy samples (Sample1, Sample2 and Sample3), with the threads per sample listed in each table.

## Preparation
Download toy datasets

The toy dataset (paired-end reads of three samples, ~7.5 GB in total) is available on Zenodo: https://zenodo.org/records/23311537

ENCORE expects one sub-folder per sample in the data folder (`-d/--data-folder`), with the reads kept gzipped:

```
Toy_Dataset/
├── Sample1/
│   ├── Sample1_R1.fastq.gz
│   └── Sample1_R2.fastq.gz
├── Sample2/
│   ├── Sample2_R1.fastq.gz
│   └── Sample2_R2.fastq.gz
└── Sample3/
    ├── Sample3_R1.fastq.gz
    └── Sample3_R2.fastq.gz
```

```bash
# Download the six read files into Toy_Dataset/<sample>/
mkdir -p Toy_Dataset
for sample in Sample1 Sample2 Sample3; do
    mkdir -p Toy_Dataset/${sample}
    for read in R1 R2; do
        wget -c -O Toy_Dataset/${sample}/${sample}_${read}.fastq.gz \
            "https://zenodo.org/records/23311537/files/${sample}_${read}.fastq.gz?download=1"
    done
done
```

The reads do not need to be decompressed: ENCORE reads the `.fastq.gz` files directly. To look at the reads, use e.g. `zcat Toy_Dataset/Sample1/Sample1_R1.fastq.gz | head -8`. Then run ENCORE with `-d ./Toy_Dataset`.


## Step1 fastqc

| Sample | R1 file size (GB) | R2 file size (GB) | Threads per sample | FastQC runtime (min) |
|:---|---:|---:|---:|---:|
| Sample1 | 0.69 | 0.72 | 1 | ~2.9 |
| Sample2 | 1.07 | 1.13 | 1 | ~4.3 |
| Sample3 | 2.02 | 1.86 | 1 | ~6.9 |

## Step2 fastp

| Sample | Reads before filtering (million) | Reads after filtering (million) | Output R1 + R2 size (GB) | Threads per sample | fastp runtime (min) |
|:---|---:|---:|---:|---:|---:|
| Sample1 | 16.8 | 16.6 | 1.25 | 16 | ~0.6 |
| Sample2 | 26.7 | 26.5 | 1.95 | 16 | ~0.9 |
| Sample3 | 45.2 | 45.0 | 3.49 | 16 | ~1.6 |

## Step3 Assembly (MEGAHIT)

| Sample | Contigs file size (MB) | Threads per sample | MEGAHIT runtime (h) |
|:---|---:|---:|---:|
| Sample1 | 34.3 | 128 | ~2 |
| Sample2 | 49.5 | 128 | ~3 |
| Sample3 | 84.2 | 128 | ~4 |
## Step4 Cross-mapping (creating depth files for binning)

| Sample | Threads per sample | Cross-mapping runtime (min) |
|:---|---:|---:|
| Sample1 | 128 | ~21 |
| Sample2 | 128 | ~24 |
| Sample3 | 128 | ~30 |

## Step5 Binning

| Method | Sample | Threads per sample | Binning runtime (min) | Bins generated |
|:---|:---|---:|---:|---:|
| MaxBin2 | Sample1 | 16 | ~3 | 47 |
| MaxBin2 | Sample2 | 16 | ~4 | 67 |
| MaxBin2 | Sample3 | 16 | ~9 | 128 |
| MetaBAT2 | Sample1 | 16 | ~2 | 120 |
| MetaBAT2 | Sample2 | 16 | ~3 | 166 |
| MetaBAT2 | Sample3 | 16 | ~11 | 263 |
| CONCOCT | Sample1 | 16 | ~2 | 36 |
| CONCOCT | Sample2 | 16 | ~2 | 35 |
| CONCOCT | Sample3 | 16 | ~3 | 38 |

## Step6 Bin refinement and reassembly (metaWRAP)

| Sample | Input bins (MaxBin2 + MetaBAT2 + CONCOCT) | Refinement threads | Refinement runtime (min) | Refined bins (completeness ≥ 50%, contamination ≤ 10%) | Reassembly threads | Reassembly runtime (h) |
|:---|---:|---:|---:|---:|---:|---:|
| Sample1 | 203 | 16 | ~64 | 27 | 128 | ~5 |
| Sample2 | 268 | 16 | ~62 | 39 | 128 | ~8 |
| Sample3 | 429 | 16 | ~70 | 59 | 128 | ~13.5 |

Reassembly is run on the refined bins.

## Step7 Extraction of protein and gene sequences


This step could be finished in around 1min per sample.

## Step8 Taxonomic classification (GTDB-Tk)

| Sample | Bins classified | Threads per sample | Classification runtime (min) |
|:---|---:|---:|---:|
| Sample1 | 27 | 16 | ~3.1 |
| Sample2 | 39 | 16 | ~3.0 |
| Sample3 | 59 | 16 | ~2.9 |

## Step9 GEM reconstruction (CarveMe, gap-filled on M8)

| Sample | Protein bins | Threads per sample | CarveMe runtime (h) | GEMs generated |
|:---|---:|---:|---:|---:|
| Sample1 | 27 | 1 | ~1 | 27 |
| Sample2 | 39 | 1 | ~1.5 | 39 |
| Sample3 | 59 | 1 | ~3 | 59 |

## Step10 ecGEM reconstruction and mandatory exchange reaction identification (GECKO 3.0)
For each bin, ecGEM is reconstructed based on the corresponding GEM. Enzyme constraints are firstly added into the model. Flux balance analysis and flux variability analysis are further performed on the ecGEM to optimize the growth rate and predict the flux range for each exchange reaction under the optimal growth rate. 
The running time for each bin with 1 threads is around 2 mins with the pre-build SMILES inhouse database. The most time consuming process is querying SMILES from PubChem via API. 

## Step11 Abundance calculation

| Sample | Bins | Threads per sample | Abundance runtime (h) |
|:---|---:|---:|---:|
| Sample1 | 27 | 128 | ~1 |
| Sample2 | 39 | 128 | ~2 |
| Sample3 | 59 | 128 | ~4.5 |


## Downstream analysis
Running times are not summarized since downstream analyses depend on user-defined parameters and thresholds. These steps are generally fast and can be run using `bash ENCORE.sh --[function]` (e.g., `bash ENCORE.sh --fva` or `bash ENCORE.sh --reporter`).
