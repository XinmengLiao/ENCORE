#!/bin/bash
set -euo pipefail

# ================================
# Default settings
# ================================
OUTPUT_DIR=""
COHORT=""
SCRIPTS_DIR="Scripts"  
CORES=1

# ================================
# Snakemake wrapper functions
# ================================

run_snakemake () {
    local smk_file=$1
    echo ">>> Running $smk_file ..."
    snakemake -s "${SCRIPTS_DIR}/${smk_file}" --cores $CORES --config OUTPUT_DIR="$OUTPUT_DIR" COHORT="$COHORT"
    echo ">>> Finished $smk_file"
}

quality()       { run_snakemake "fastqc.smk"; }
trim()          { run_snakemake "fastp.smk"; }
assembly()      { run_snakemake "megahit.smk"; }
crossmap()      { run_snakemake "crossmapping.smk"; }
maxbin()        { run_snakemake "maxbin.smk"; }
concoct()       { run_snakemake "concoct.smk"; }
metabat()       { run_snakemake "metabat.smk"; }
refinement()    { run_snakemake "refinement.smk"; }
reassembly()    { run_snakemake "reassembly.smk"; }
abundance()     { run_snakemake "AbundanceCalculation.smk"; }
taxonomy()      { run_snakemake "taxonomy.smk"; }
daa()           { run_snakemake "daa.smk"; }

extraction()    { run_snakemake "extraction.smk"; }
gem()           { run_snakemake "carveme.smk"; }
ecgem()         { run_snakemake "GEMtoECGEM.smk"; }

network()       { run_snakemake "CoGEMNetwork.smk"; }
reporter()      { run_snakemake "ReporterMetabolite.smk"; }

all_pipeline()  { run_snakemake "WholePipeline.smk"; }

# ================================
# Usage function
# ================================
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    
    -a, --all           Run the whole ENCORE pipeline
    -o, --output_dir    Output directory
    --cohort            Cohort name prefix
    -t, --thread        Threads number for running jobs 
    -h, --help          Display this help message
    

    ---- Metagenome --------------------------

    --quality           Read quality check
    --trim              Read trimming
    --assembly          Contig assembly
    --crossmap          Prepare depth files for binning
    --maxbin            Binning with MaxBin2
    --concoct           Binning with CONCOCT
    --metabat           Binning with MetaBat2
    --refinement        Bin refinement
    --reassembly        Bin reassembly
    --abundance         Abundance Calculation
    --taxonomy          Bin classifications to species level
    -daa                Differential Abundance Analysis
    
    ------------------------------------------
    

    ---- Genome-scale Metabolic Modelling ----

    --extract           Extract genomic and proteomic sequences
    --gem               Reconstruction of GEMs
    --ecgem             Reconstruction of ecGEMs
    ------------------------------------------
    

    ---- Reporter Metabolite -----------------

    --network           Construct microbial community ecGEMs network
    --reporter          Identify reporter metabolites
    ------------------------------------------


REQUIRED ARGUMENTS:
    - input-sample: Sample ID
    - output-directory: Output directory path
    - vcf: User uploaded vcf file
EOF
}

# ================================
# Parse arguments
# ================================
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-directory)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--cohort)
            COHORT="$2"
            shift 2
            ;;
        -s|--scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        -j|--cores)
            CORES="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            COMMAND="$1"
            shift
            ;;
    esac
done

# ================================
# Validate required arguments
# ================================
if [[ -z "${COMMAND:-}" ]]; then
    echo "Error: No command specified"
    usage
    exit 1
fi

# ================================
# Run the chosen command
# ================================
case "$COMMAND" in
    quality|trim|assembly|crossmap|maxbin|concoct|metabat|refinement|reassembly|abundance|taxonomy|daa|extraction|gem|ecgem|network|reporter|all)
        $COMMAND
        ;;
    *)
        echo "Error: Unknown command $COMMAND"
        usage
        exit 1
        ;;
esac
