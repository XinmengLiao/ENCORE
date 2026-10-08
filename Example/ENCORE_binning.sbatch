#!/bin/bash
#SBATCH -A naiss2026-4-670
#SBATCH -p main
#SBATCH -N 1
#SBATCH -c 128
#SBATCH -t 12:00:00
#SBATCH -J encore_binning
#SBATCH -o encore_binning_%j.out

# Driver job: launches maxbin, metabat and concoct side by side.
# With --slurm, ENCORE.sh submits every Snakemake job of each module via sbatch
# (resources from ENCORE_config.yaml), so this job itself only needs 1 core.
cd /cfs/klemming/projects/supr/naiss2024-22-1134/ENCORE/Example || exit 1

for mod in maxbin metabat concoct; do
    bash ENCORE.sh -o ./output -c test -d Toy_Dataset --$mod --slurm > encore_${mod}_slurm.log 2>&1 &
done

wait
