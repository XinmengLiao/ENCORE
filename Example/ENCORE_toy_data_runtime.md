## Step1 fastqc

| Sample | R1 file size (GB) | R2 file size (GB) | Threads per sample | FastQC runtime (min) |
|:---|---:|---:|---:|---:|
| Sample1 | 0.69 | 0.72 | 1 | ~2.9 |
| Sample2 | 1.07 | 1.13 | 1 | ~4.3 |
| Sample3 | 2.02 | 1.86 | 1 | ~6.9 |

## Step2 fastp

| Sample | Reads before filtering (million) | Reads after filtering (million) | Output R1 + R2 size (GB) | fastp runtime (min) |
|:---|---:|---:|---:|---:|
| Sample1 | 16.8 | 16.6 | 1.25 | ~0.6 |
| Sample2 | 26.7 | 26.5 | 1.95 | ~0.9 |
| Sample3 | 45.2 | 45.0 | 3.49 | ~1.6 |
