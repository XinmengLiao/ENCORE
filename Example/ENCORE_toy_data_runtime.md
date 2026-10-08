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
| Sample1 | 34.3 | 256 | ~2 |
| Sample2 | 49.5 | 256 | ~3 |
| Sample3 | 84.2 | 256 | ~4 |
## Step4 Cross-mapping (depth files for binning)

| Sample | Threads per sample | Cross-mapping runtime (min) |
|:---|---:|---:|
| Sample1 | 256 | ~21 |
| Sample2 | 256 | ~24 |
| Sample3 | 256 | ~30 |
