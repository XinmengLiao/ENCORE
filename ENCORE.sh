#!/bin/bash
set -euo pipefail

# ================================
# Default settings
# ================================
OUTPUT_DIR=""
COHORT=""
DATA_FOLDER=""
SCRIPTS_DIR="Scripts"  
CORES=1
DRY_RUN=false
QUIET=false
KEEP_INCOMPLETE=false
declare -a MODULES=()

# ================================
# Module definitions (using simple arrays for zsh compatibility)
# ================================
PREPROCESSING_MODULES="quality trim assembly"
DOWNSTREAM_MODULES="crossmap maxbin concoct metabat refinement reassembly abundance taxonomy extraction gem ecgem"
ALL_MODULES="quality trim assembly crossmap maxbin concoct metabat refinement reassembly abundance taxonomy extraction gem ecgem network reporter"

# Function to get snakemake file for a module
get_smk_file() {
    case $1 in
        quality) echo "fastqc.smk" ;;
        trim) echo "fastp.smk" ;;
        assembly) echo "megahit.smk" ;;
        preprocessing) echo "preprocessing.smk" ;;
        crossmap) echo "crossmapping.smk" ;;
        maxbin) echo "maxbin.smk" ;;
        concoct) echo "concoct.smk" ;;
        metabat) echo "metabat.smk" ;;
        refinement) echo "refinement.smk" ;;
        reassembly) echo "reassembly.smk" ;;
        abundance) echo "AbundanceCalculation.smk" ;;
        taxonomy) echo "taxonomy.smk" ;;
        extraction) echo "extraction.smk" ;;
        gem) echo "carveme.smk" ;;
        ecgem) echo "GEMtoECGEM.smk" ;;
        downstream) echo "downstream.smk" ;;
        network) echo "CoGEMNetwork.smk" ;;
        reporter) echo "ReporterMetabolite.smk" ;;
        *) return 1 ;;
    esac
}

# Function to get description for a module
get_module_desc() {
    case $1 in
        quality) echo "Read quality check" ;;
        trim) echo "Read trimming" ;;
        assembly) echo "Contig assembly" ;;
        preprocessing) echo "Complete preprocessing (quality, trim, assembly)" ;;
        crossmap) echo "Prepare depth files for binning" ;;
        maxbin) echo "Binning with MaxBin2" ;;
        concoct) echo "Binning with CONCOCT" ;;
        metabat) echo "Binning with MetaBat2" ;;
        refinement) echo "Bin refinement" ;;
        reassembly) echo "Bin reassembly" ;;
        abundance) echo "Abundance Calculation" ;;
        taxonomy) echo "Bin classifications to species level" ;;
        extraction) echo "Extract genomic and proteomic sequences" ;;
        gem) echo "Reconstruction of GEMs" ;;
        ecgem) echo "Reconstruction of ecGEMs" ;;
        downstream) echo "Complete downstream analysis (binning to ecGEM)" ;;
        network) echo "Construct microbial community ecGEMs network" ;;
        reporter) echo "Identify reporter metabolites" ;;
        *) return 1 ;;
    esac
}

# ================================
# Helper functions
# ================================
log_info() {
    echo -e "[INFO] $*"

}

log_success() {
    echo -e "[✓] $*"
}

log_warning() {
    echo -e "[WARNING] $*"
}

log_error() {
    echo -e "[ERROR] $*"
}

validate_env() {
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "Output directory is required (-o/--output-dir)"
        return 1
    fi
    
    if [[ -z "$COHORT" ]]; then
        log_error "Cohort name is required (-c/--cohort)"
        return 1
    fi
    
    if [[ -z "$DATA_FOLDER" ]]; then
        log_error "Data folder is required (-d/--data-folder)"
        return 1
    fi
    
    if ! command -v snakemake &> /dev/null; then
        log_error "Snakemake is not installed or not in PATH"
        return 1
    fi
    
    if [[ ! -d "$SCRIPTS_DIR/Main_Functions" ]]; then
        log_error "Scripts directory not found: $SCRIPTS_DIR/Main_Functions"
        return 1
    fi
    
    return 0
}

validate_module() {
    local module=$1
    # Check if module is in any of the module lists
    if [[ "$ALL_MODULES" == *"$module"* ]]; then
        return 0
    fi
    return 1
}

list_modules() {
    log_info "Available modules:"
    echo ""
    echo -e "All available modules:"
    for mod in $ALL_MODULES; do
        printf "  %-15s %s\n" "$mod" "$(get_module_desc $mod)"
    done
}

# ================================
# Snakemake wrapper functions
# ================================

run_snakemake() {
    local module=$1
    local smk_file
    smk_file=$(get_smk_file "$module") || {
        log_error "Unknown module: $module"
        return 1
    }
    
    if [[ ! -f "${SCRIPTS_DIR}/Main_Functions/${smk_file}" ]]; then
        log_error "Snakemake file not found: ${SCRIPTS_DIR}/Main_Functions/${smk_file}"
        return 1
    fi
    
    local desc
    desc=$(get_module_desc "$module")
    log_info "Running: $desc ($module)"
    
    local snakemake_opts="--cores $CORES --config OUTPUT_DIR=\"$OUTPUT_DIR\" COHORT=\"$COHORT\" DATA_FOLDER=\"$DATA_FOLDER\""
    
    if [[ "$DRY_RUN" == true ]]; then
        snakemake_opts="$snakemake_opts --dry-run"
    fi
    
    if [[ "$KEEP_INCOMPLETE" == true ]]; then
        snakemake_opts="$snakemake_opts --keep-incomplete"
    fi
    
    eval "snakemake -s \"${SCRIPTS_DIR}/Main_Functions/${smk_file}\" $snakemake_opts"
    
    if [[ $? -eq 0 ]]; then
        log_success "Completed: $module"
        return 0
    else
        log_error "Failed: $module"
        return 1
    fi
}

run_modules() {
    local failed_modules=()
    
    log_info "Running ${#MODULES[@]} module(s)"
    echo ""
    
    for module in "${MODULES[@]}"; do
        if ! validate_module "$module"; then
            log_error "Unknown module: $module"
            failed_modules+=("$module")
            continue
        fi
        
        if ! run_snakemake "$module"; then
            failed_modules+=("$module")
        fi
        echo ""
    done
    
    # Summary
    echo -e "=== Pipeline Summary ==="
    local completed=$((${#MODULES[@]} - ${#failed_modules[@]}))
    log_success "Completed: $completed/${#MODULES[@]} module(s)"
    
    if [[ ${#failed_modules[@]} -gt 0 ]]; then
        log_error "Failed modules: ${failed_modules[*]}"
        return 1
    fi
    
    return 0
}


quality()       { run_snakemake "quality"; }
trim()          { run_snakemake "trim"; }
assembly()      { run_snakemake "assembly"; }
preprocessing() { run_snakemake "preprocessing"; }
crossmap()      { run_snakemake "crossmap"; }
maxbin()        { run_snakemake "maxbin"; }
concoct()       { run_snakemake "concoct"; }
metabat()       { run_snakemake "metabat"; }
refinement()    { run_snakemake "refinement"; }
reassembly()    { run_snakemake "reassembly"; }
abundance()     { run_snakemake "abundance"; }
taxonomy()      { run_snakemake "taxonomy"; }

extract()       { run_snakemake "extract"; }
gem()           { run_snakemake "gem"; }
ecgem()         { run_snakemake "ecgem"; }

downstream()    { run_snakemake "downstream"; }
network()       { run_snakemake "network"; }
reporter()      { run_snakemake "reporter"; }

# ================================
# Usage function
# ================================
usage() {
    cat << EOF
╔═══════════════════════════════════════════════════════════════╗
║          ENCORE - Metagenomic Pipeline                          ║
║   Enhanced metagenome analysis with metabolic modelling         ║
╚═══════════════════════════════════════════════════════════════╝

Usage:
  $0 [OPTIONS] <modules>

OPTIONS:
  -o, --output-dir <dir>      Output directory (required)
  -c, --cohort <name>         Cohort/sample prefix (required)
  -d, --data-folder <name>    Data folder name (default: Toy_Dataset)
  -t, --cores <num>           Number of CPU cores (default: 1)
  -s, --scripts-dir <dir>     Scripts directory (default: Scripts)
  
  -n, --dry-run               Perform a dry run without executing
  -k, --keep-incomplete       Keep incomplete output files (don't delete intermediate files)
  -l, --list                  List all available modules
  -h, --help                  Show this help message

EXAMPLES:
  # Run a single module
  bash $0 -o ./output -c sample1 -d my_data --quality

  # Run multiple modules in sequence
  bash $0 -o ./output -c sample1 -d my_data --quality --trim --assembly

  # Run an entire workflow
  bash $0 -o ./output -c sample1 -d my_data --metagenome --metabolic

  # Dry run to see what would be executed
  bash $0 -o ./output -c sample1 --dry-run --assembly

  # Keep intermediate files if job fails
  bash $0 -o ./output -c sample1 -d my_data --assembly --keep-incomplete

  # List all available modules
  bash $0 --list

WORKFLOW GROUPS:
  --preprocessing Run preprocessing modules (quality, trim, assembly)
  --downstream    Run downstream analysis modules (crossmap through ecgem)
  --all           Run complete pipeline (preprocessing → crossmap → downstream → network → reporter)

INDIVIDUAL MODULES (Metagenome):
  --quality       Read quality check
  --trim          Read trimming
  --assembly      Contig assembly
  --preprocessing Complete preprocessing (quality, trim, assembly)
  --crossmap      Prepare depth files for binning
  --maxbin        Binning with MaxBin2
  --concoct       Binning with CONCOCT
  --metabat       Binning with MetaBat2
  --refinement    Bin refinement
  --reassembly    Bin reassembly
  --abundance     Abundance Calculation
  --taxonomy      Bin classifications to species level

INDIVIDUAL MODULES (Metabolic):
  --extract       Extract genomic and proteomic sequences
  --gem           Reconstruction of GEMs
  --ecgem         Reconstruction of ecGEMs
  --downstream    Complete downstream analysis (binning to ecGEM)

INDIVIDUAL MODULES (Reporter):
  --network       Construct microbial community ecGEMs network
  --reporter      Identify reporter metabolites

EOF
}

# ================================
# Parse arguments
# ================================
if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--cohort)
            COHORT="$2"
            shift 2
            ;;
        -d|--data-folder)
            DATA_FOLDER="$2"
            shift 2
            ;;
        -s|--scripts-dir)
            SCRIPTS_DIR="$2"
            shift 2
            ;;
        -t|--cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -k|--keep-incomplete)
            KEEP_INCOMPLETE=true
            shift
            ;;
        -l|--list)
            list_modules
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --preprocessing)
            for mod in $PREPROCESSING_MODULES; do
                MODULES+=("$mod")
            done
            shift
            ;;
        --downstream)
            for mod in $DOWNSTREAM_MODULES; do
                MODULES+=("$mod")
            done
            shift
            ;;
        --reporter)
            for mod in $REPORTER_MODULES; do
                MODULES+=("$mod")
            done
            shift
            ;;
        --all)
            for mod in $ALL_MODULES; do
                MODULES+=("$mod")
            done
            shift
            ;;
        --quality|--trim|--assembly|--preprocessing|--crossmap|--maxbin|--concoct|--metabat|--refinement|--reassembly|--abundance|--taxonomy|--extract|--gem|--ecgem|--downstream|--network|--reporter)
            # Remove leading dashes
            module="${1#--}"
            # Handle alias: --extract -> extraction
            if [[ "$module" == "extract" ]]; then
                module="extraction"
            fi
            MODULES+=("$module")
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# ================================
# Validate and execute
# ================================
if ! validate_env; then
    exit 1
fi

if [[ ${#MODULES[@]} -eq 0 ]]; then
    log_error "No modules specified"
    echo ""
    usage
    exit 1
fi

log_info "ENCORE Pipeline - Configuration"
log_info "Output directory: $OUTPUT_DIR"
log_info "Cohort: $COHORT"
log_info "Data folder: $DATA_FOLDER"
log_info "CPU cores: $CORES"
[[ "$DRY_RUN" == true ]] && log_warning "DRY RUN MODE - No modules will be executed"
[[ "$KEEP_INCOMPLETE" == true ]] && log_warning "KEEP INCOMPLETE MODE - Intermediate files will be preserved"
echo ""

if ! run_modules; then
    exit 1
fi
