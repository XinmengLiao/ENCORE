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

Each sample was run as one job, and its bins were carved one after another. CarveMe uses a single thread (the CPU time Slurm recorded matches the wall-clock time, so about 1 core was busy). The `cores: carveme: 16` setting was not applied in this run. Slurm allocated a whole node (256 CPUs) to each job, but only about 1 core was used.

## Summary of total runtime

| Step | Threads per sample | Sample1 | Sample2 | Sample3 |
|:---|---:|---:|---:|---:|
| Step1 FastQC | 1 | ~2.9 min | ~4.3 min | ~6.9 min |
| Step2 fastp | 16 | ~0.6 min | ~0.9 min | ~1.6 min |
| Step3 Assembly (MEGAHIT) | 128 | ~2 h | ~3 h | ~4 h |
| Step4 Cross-mapping | 128 | ~21 min | ~24 min | ~30 min |
| Step5 Binning (3 binners, summed) | 16 | ~7 min | ~9 min | ~23 min |
| Step6 Bin refinement (16 threads) + reassembly (128 threads) | 16 / 128 | ~64 min + ~5 h | ~62 min + ~8 h | ~70 min + ~13.5 h |
| Step7 Extraction (protein + DNA) | 1 | ~1 min | ~1 min | ~1 min |
| Step8 Classification (GTDB-Tk) | 16 | ~3.1 min | ~3.0 min | ~2.9 min |
| Step9 GEM reconstruction (CarveMe) | 1 | ~1 h (27 GEMs) | ~1.5 h (39 GEMs) | ~3 h (59 GEMs) |
| **Total (approx.)** | | **~9.6 h** | **~14.2 h** | **~22.7 h** |

Note: binning time assumes the three binners run sequentially.
