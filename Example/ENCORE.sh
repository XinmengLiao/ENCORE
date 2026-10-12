#!/bin/bash
set -euo pipefail

# ================================
# Default settings
# ================================
OUTPUT_DIR=""
COHORT=""
DATA_FOLDER=""
METADATA_FILE=""
FIXED_EFFECT="Group"
RANDOM_EFFECT="NA"
CHECKM_FILE="NA"
SCRIPT_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Scripts/ is next to this file (repository root) or one level up (Example/)
if [[ -d "${SCRIPT_HOME}/Scripts" ]]; then
    SCRIPTS_DIR="${SCRIPT_HOME}/Scripts"
else
    SCRIPTS_DIR="${SCRIPT_HOME}/../Scripts"
fi
CONFIG_FILE="${SCRIPT_HOME}/ENCORE_config.yaml"
CORES=1
DRY_RUN=false
KEEP_INCOMPLETE=false
declare -a MODULES=()

# ================================
# Module definitions (using simple arrays for zsh compatibility)
# ================================
# extraction must run before taxonomy (GTDB-Tk reads dna_bins/)
PREPROCESSING_MODULES="quality trim assembly"
DOWNSTREAM_MODULES="crossmap maxbin concoct metabat refinement reassembly abundance extraction taxonomy gem ecgem"
# --all runs up to fva; reporter needs the sample metadata (--metadata) and is run on its own
ALL_MODULES="$PREPROCESSING_MODULES $DOWNSTREAM_MODULES fva reporter"

# Function to get snakemake file for a module
get_smk_file() {
    case $1 in
        quality) echo "fastqc.smk" ;;
        trim) echo "fastp.smk" ;;
        assembly) echo "megahit.smk" ;;
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
        *) return 1 ;;
    esac
}

# Function to get description for a module
get_module_desc() {
    case $1 in
        quality) echo "Read quality check" ;;
        trim) echo "Read trimming" ;;
        assembly) echo "Contig assembly" ;;
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
        fva) echo "Collect FVA results and extract mandatory exchange reactions" ;;
        reporter) echo "DAA, microbe-metabolite network and reporter metabolites" ;;
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
    
    # fva and reporter run in plain bash (see fva_data_manage, run_reporter_bash);
    # every other module needs Snakemake
    if needs_snakemake && ! command -v snakemake &> /dev/null; then
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

# Read a key from a top-level block (e.g. 'folder:' or 'cores:') of the config file
read_config_value() {
    awk -v blk="$1" -v key="$2" '
        $0 ~ "^"blk":" {in_blk=1; next}
        /^[^[:space:]#]/ {in_blk=0}
        in_blk && $1 == key":" {sub(/[[:space:]]*#.*/, "", $0); sub(/^[^:]*:[[:space:]]*/, "", $0); print; exit}
    ' "$CONFIG_FILE"
}

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
    
    local snakemake_opts="--cores $CORES --configfile \"$CONFIG_FILE\" --config OUTPUT_DIR=\"$OUTPUT_DIR\" COHORT=\"$COHORT\" DATA_FOLDER=\"$DATA_FOLDER\" CONFIG_DIR=\"$(dirname "$CONFIG_FILE")\""
    
    if [[ "$DRY_RUN" == true ]]; then
        snakemake_opts="$snakemake_opts --dry-run"
    fi
    
    if [[ "$KEEP_INCOMPLETE" == true ]]; then
        snakemake_opts="$snakemake_opts --keep-incomplete"
    fi

    # Rules without a 'threads:' directive default to 1 thread,
    # but their shell commands use cores:<key> from the config. Set it explicitly.
    case "$module" in
        crossmap) snakemake_opts="$snakemake_opts --set-threads crossMapSeries=$(read_config_value cores crossMap)" ;;
        refinement) snakemake_opts="$snakemake_opts --set-threads binRefine=$(read_config_value cores refine)" ;;
        reassembly) snakemake_opts="$snakemake_opts --set-threads binReassemble=$(read_config_value cores reassemble)" ;;
        gem) snakemake_opts="$snakemake_opts --set-threads carveme=$(read_config_value cores carveme)" ;;
        abundance) snakemake_opts="$snakemake_opts --set-threads abundance=$(read_config_value cores abundance)" ;;
        taxonomy) snakemake_opts="$snakemake_opts --set-threads GTDBtk=$(read_config_value cores gtdbtk)" ;;
    esac

    # Modules listed under 'module_envs:' in the config run inside that conda env,
    # stacked on top of the current one (conda activate --stack).
    local stack_key stack_env stacked=false
    stack_key=$(read_config_value module_envs "$module")
    if [[ -n "$stack_key" ]]; then
        stack_env=$(read_config_value envs "$stack_key")
        if [[ -z "$stack_env" || ! -d "$stack_env" ]]; then
            log_error "Env '$stack_key' for module $module not found in 'envs:' of $CONFIG_FILE"
            return 1
        fi
        set +u
        eval "$(conda shell.bash hook)"
        conda activate --stack "$stack_env" || { set -u; log_error "conda activate --stack $stack_env failed"; return 1; }
        set -u
        stacked=true
        log_info "Activated (stacked) conda env: $stack_env"
    fi

    local rc=0
    eval "snakemake -s \"${SCRIPTS_DIR}/Main_Functions/${smk_file}\" $snakemake_opts" || rc=$?

    if [[ "$stacked" == true ]]; then
        set +u; conda deactivate; set -u
    fi

    if [[ $rc -eq 0 ]]; then
        log_success "Completed: $module"
        return 0
    else
        log_error "Failed: $module"
        return 1
    fi
}

# Does any selected module still go through Snakemake?
needs_snakemake() {
    local m
    [[ ${#MODULES[@]} -eq 0 ]] && return 1
    for m in "${MODULES[@]}"; do
        [[ "$m" != "fva" && "$m" != "reporter" ]] && return 0
    done
    return 1
}

# Config value as a path: strip quotes, expand ~, resolve relative to the config file folder
config_path() {
    local val
    val=$(read_config_value "$1" "$2")
    val="${val%\"}"; val="${val#\"}"; val="${val%\'}"; val="${val#\'}"
    val="${val/#\~/$HOME}"
    [[ -z "$val" ]] && return 0
    [[ "$val" != /* ]] && val="$(dirname "$CONFIG_FILE")/$val"
    echo "$val"
}

# ================================
# FVA results of the ecGEMs (plain bash, no Snakemake)
# meta-gecko.m writes the FVA of the exchange reactions of each bin to
# <ecGEMs>/<sample>/<bin>/output/<sample>.<binID>FVA.txt (exchangeRxns, equations, LB_exc, UB_exc).
# Copy them to <output>/<fva>/ and combine the mandatory exchange reactions of all samples and bins
# (flux range does not include 0: UB_exc < 0 Uptake, LB_exc > 0 Secretion)
# into <output>/<fva>/fva_mandatory_exchange.txt with Sample and Bin columns.
# ================================
fva_data_manage() {
    log_info "Running: $(get_module_desc fva) (fva) in bash, without Snakemake"

    local ecgem_root fva_root combined
    ecgem_root="$OUTPUT_DIR/$(read_config_value folder ecGEMs)"
    fva_root="$OUTPUT_DIR/$(read_config_value folder fva)"
    combined="$fva_root/fva_mandatory_exchange.txt"
    [[ -n "$(read_config_value folder fva)" ]] || { log_error "folder: fva missing in $CONFIG_FILE"; return 1; }
    [[ -d "$ecgem_root" ]] || { log_error "ecGEM folder not found: $ecgem_root (run --ecgem first)"; return 1; }

    local sample_dir bin_dir sample bin fva_file
    local files=() samples=() bins=()
    for sample_dir in "$ecgem_root"/*/; do
        [[ -d "$sample_dir" ]] || continue
        sample=$(basename "$sample_dir")
        for bin_dir in "$sample_dir"*/; do
            [[ -d "${bin_dir}output" ]] || continue
            bin=$(basename "$bin_dir")
            for fva_file in "${bin_dir}output/"*FVA.txt; do
                [[ -f "$fva_file" ]] || continue
                files+=("$fva_file"); samples+=("$sample"); bins+=("$bin")
            done
        done
    done
    [[ ${#files[@]} -gt 0 ]] || { log_error "No FVA files (*FVA.txt) found in $ecgem_root/<sample>/<bin>/output/"; return 1; }
    log_info "Found ${#files[@]} FVA file(s) in $(printf '%s\n' "${samples[@]}" | sort -u | wc -l | tr -d ' ') sample(s)"

    if [[ "$DRY_RUN" == true ]]; then
        log_info "[dry-run] would copy them to $fva_root and write $combined"
        return 0
    fi

    mkdir -p "$fva_root"
    printf 'exchangeRxns\tequations\tLB_exc\tUB_exc\tDirection\tSample\tBin\n' > "$combined"
    local i
    for i in "${!files[@]}"; do
        cp "${files[$i]}" "$fva_root/" || { log_error "Copy failed: ${files[$i]}"; return 1; }
        # Columns are found by name in the header; bin is the folder name (with dots, e.g. bin.1.s)
        awk -F '\t' -v OFS='\t' -v sample="${samples[$i]}" -v bin="${bins[$i]}" '
            NR == 1 { for (c = 1; c <= NF; c++) col[$c] = c; next }
            {
                lb = $(col["LB_exc"]) + 0; ub = $(col["UB_exc"]) + 0
                if (ub < 0) dir = "Uptake"; else if (lb > 0) dir = "Secretion"; else next
                print $(col["exchangeRxns"]), $(col["equations"]), $(col["LB_exc"]), $(col["UB_exc"]), dir, sample, bin
            }' "${files[$i]}" >> "$combined"
    done

    log_success "${#files[@]} FVA file(s) copied to $fva_root"
    log_success "$(($(wc -l < "$combined") - 1)) mandatory exchange reaction(s) written to $combined"
    return 0
}

# ================================
# Reporter metabolites (plain bash, no Snakemake)
# Runs Scripts/RScripts/DAA_RM.R: DAA of the species with MaAsLin2, microbe-metabolite network
# from the mandatory exchange reactions (--fva), and directional reporter metabolites.
# Inputs: <abundance>/<sample>/*.map, <classification>/<sample>/*_gtdbtk_summary.tsv,
# <fva>/fva_mandatory_exchange.txt and the --metadata file. Output: <output>/<reporter>/
# ================================
run_reporter_bash() {
    log_info "Running: $(get_module_desc reporter) (reporter) in bash, without Snakemake"

    local abund_root taxo_root fva_file out_dir daa_rm
    abund_root="$OUTPUT_DIR/$(read_config_value folder abundance)"
    taxo_root="$OUTPUT_DIR/$(read_config_value folder classification)"
    fva_file="$OUTPUT_DIR/$(read_config_value folder fva)/fva_mandatory_exchange.txt"
    out_dir="$OUTPUT_DIR/$(read_config_value folder reporter)"
    daa_rm=$(config_path Rscripts DAA_RM)

    [[ -f "$daa_rm" ]] || { log_error "Rscripts: DAA_RM not found in $CONFIG_FILE: $daa_rm"; return 1; }
    [[ -f "$METADATA_FILE" ]] || { log_error "Metadata file not found: $METADATA_FILE (--metadata)"; return 1; }
    [[ "$CHECKM_FILE" == "NA" || -f "$CHECKM_FILE" ]] || { log_error "CheckM file not found: $CHECKM_FILE (--checkm)"; return 1; }
    [[ -d "$abund_root" ]] || { log_error "Abundance folder not found: $abund_root (run --abundance first)"; return 1; }
    [[ -d "$taxo_root" ]] || { log_error "Taxonomy folder not found: $taxo_root (run --taxonomy first)"; return 1; }
    [[ -f "$fva_file" ]] || { log_error "FVA file not found: $fva_file (run --fva first)"; return 1; }
    command -v Rscript &> /dev/null || { log_error "Rscript is not installed or not in PATH"; return 1; }

    log_info "Metadata: $METADATA_FILE (fixed effect: $FIXED_EFFECT, random effect: $RANDOM_EFFECT, CheckM: $CHECKM_FILE)"
    if [[ "$DRY_RUN" == true ]]; then
        log_info "[dry-run] would run $daa_rm and write the results to $out_dir"
        return 0
    fi

    mkdir -p "$out_dir/taxonomy"
    cp "$taxo_root"/*/*_gtdbtk_summary.tsv "$out_dir/taxonomy/" || { log_error "No *_gtdbtk_summary.tsv in $taxo_root/<sample>/"; return 1; }

    # Mapped reads and bin length of each bin from the flagstat .map files of --abundance;
    # bin names as in the ecGEMs (bin.1.orig -> bin.1.o, bin.1.strict -> bin.1.s, bin.1.permissive -> bin.1.p)
    local summary="$out_dir/abundance.summary.txt" map_file sample bin
    printf 'Sample\tBin\tMapped_reads\tPrimary_mapped\tBin_length\n' > "$summary"
    for map_file in "$abund_root"/*/*.map; do
        [[ -f "$map_file" ]] || continue
        sample=$(basename "$(dirname "$map_file")")
        bin=$(basename "$map_file" .map | sed -E 's/\.orig$/.o/; s/\.strict$/.s/; s/\.permissive$/.p/')
        awk -v OFS='\t' -v sample="$sample" -v bin="$bin" '
            / mapped \(/ && !/primary/ { mapped = $1 }
            /primary mapped \(/ { primary = $1 }
            /^Bin Length =/ { len = $4 }
            END { print sample, bin, mapped, primary, len }' "$map_file" >> "$summary"
    done
    [[ $(wc -l < "$summary") -gt 1 ]] || { log_error "No .map files in $abund_root/<sample>/"; return 1; }

    log_info "Now running DAA, microbe-metabolite network and reporter metabolites ..."
    Rscript "$daa_rm" "$summary" "$out_dir/taxonomy" "$METADATA_FILE" "$fva_file" "$out_dir" \
        "$FIXED_EFFECT" "$RANDOM_EFFECT" "$CHECKM_FILE" 2>&1 | tee "$out_dir/DAA_RM.Rlog"
    [[ ${PIPESTATUS[0]} -eq 0 ]] || { log_error "DAA_RM.R failed, see $out_dir/DAA_RM.Rlog"; return 1; }

    log_success "Reporter metabolites written to $out_dir"
    return 0
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
        
        if [[ "$module" == "fva" ]]; then
            fva_data_manage || failed_modules+=("$module")
        elif [[ "$module" == "reporter" ]]; then
            run_reporter_bash || failed_modules+=("$module")
        elif ! run_snakemake "$module"; then
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
crossmap()      { run_snakemake "crossmap"; }
maxbin()        { run_snakemake "maxbin"; }
concoct()       { run_snakemake "concoct"; }
metabat()       { run_snakemake "metabat"; }
refinement()    { run_snakemake "refinement"; }
reassembly()    { run_snakemake "reassembly"; }
abundance()     { run_snakemake "abundance"; }
taxonomy()      { run_snakemake "taxonomy"; }

extraction()    { run_snakemake "extraction"; }
gem()           { run_snakemake "gem"; }
ecgem()         { run_snakemake "ecgem"; }
fva()           { fva_data_manage; }

reporter()      { run_reporter_bash; }

# ================================
# Usage function
# ================================
usage() {
    cat << EOF
ENCORE - enzyme-constrained community metabolic models from metagenomes

Usage:
  bash $0 -o <output_dir> -c <name> -d <data_folder> [OPTIONS] <modules>

REQUIRED:
  -o, --output-dir <dir>      Output directory
  -c, --cohort <name>         Name of the run
  -d, --data-folder <dir>     Folder with the raw reads, one sub-folder per sample:
                              <sample>/<sample>_R1.fastq.gz and <sample>/<sample>_R2.fastq.gz

OPTIONS:
  -t, --cores <num>           Maximum number of CPU cores used by Snakemake (default: 1);
                              threads per tool are set under 'cores:' in ENCORE_config.yaml
  -s, --scripts-dir <dir>     Scripts directory (default: Scripts/ of this repository)
  -n, --dry-run               Show what would be run without executing
  -k, --keep-incomplete       Keep incomplete output files of failed jobs
  -l, --list                  List all modules
  -h, --help                  Show this help message

REPORTER OPTIONS (--reporter):
  --metadata <file>           Sample metadata (required): tab-separated, first column = sample ID
  --fixed-effect <col>        Metadata column with the two groups compared (default: Group);
                              extra covariates comma separated, e.g. Group,Age
  --random-effect <col>       Metadata column for paired / repeated samples (default: none)
  --checkm <file>             Bin quality (Sample, bin, completeness, contamination); keeps bins with
                              completeness > 70 and contamination < 10 (default: no filtering)

WORKFLOW GROUPS:
  --preprocessing   quality, trim, assembly
  --downstream      crossmap, maxbin, concoct, metabat, refinement, reassembly, abundance,
                    extraction, taxonomy, gem, ecgem
  --all             --preprocessing, --downstream and --fva (raw reads to the mandatory exchange reactions)

MODULES (Metagenome):
  --quality         FastQC - read quality check
  --trim            fastp - read trimming
  --assembly        MEGAHIT - contig assembly
  --crossmap        Prepare depth files for binning
  --maxbin          Binning with MaxBin2
  --concoct         Binning with CONCOCT
  --metabat         Binning with MetaBAT2
  --refinement      metaWRAP - bin refinement
  --reassembly      metaWRAP - bin reassembly
  --abundance       Abundance of each bin (bwa, samtools)
  --extraction      Extract DNA and protein sequences of the bins
  --taxonomy        GTDB-Tk - taxonomic classification of the bins (after --extraction)

MODULES (Metabolic):
  --gem             CarveMe - reconstruction of GEMs
  --ecgem           GECKO 3 - reconstruction of ecGEMs and FVA of their exchange reactions
                    (runs MATLAB in batch mode: 'matlab' must be on PATH, see 'matlab:' in the config)
  --fva             Copy the FVA results of the ecGEMs to <output>/fva and combine the
                    mandatory exchange reactions into <output>/fva/fva_mandatory_exchange.txt
  --reporter        DAA (MaAsLin2), microbe-metabolite network and reporter metabolites
                    (needs --metadata and the output of --abundance, --taxonomy and --fva)

EXAMPLES:
  # Raw reads to the mandatory exchange reactions of the ecGEMs
  bash $0 -o ./output -c test -d ./Toy_Dataset -t 16 --all

  # Single modules
  bash $0 -o ./output -c test -d ./Toy_Dataset --quality --trim

  # Dry run
  bash $0 -o ./output -c test -d ./Toy_Dataset --dry-run --assembly

  # DAA, microbe-metabolite network and reporter metabolites
  bash $0 -o ./output -c test -d ./Toy_Dataset --reporter \\
      --metadata metadata.tsv --fixed-effect Group --random-effect subject
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
        --metadata)
            METADATA_FILE="$2"
            shift 2
            ;;
        --fixed-effect)
            FIXED_EFFECT="$2"
            shift 2
            ;;
        --random-effect)
            RANDOM_EFFECT="$2"
            shift 2
            ;;
        --checkm)
            CHECKM_FILE="$2"
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
        --all)
            for mod in $PREPROCESSING_MODULES $DOWNSTREAM_MODULES fva; do
                MODULES+=("$mod")
            done
            shift
            ;;
        --quality|--trim|--assembly|--crossmap|--maxbin|--concoct|--metabat|--refinement|--reassembly|--abundance|--taxonomy|--extraction|--gem|--ecgem|--fva|--reporter)
            # Remove leading dashes
            module="${1#--}"
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

# Put the conda envs listed in the config (envs: metagem / metawrap) on PATH
for env_key in metagem metawrap; do
    env_dir=$(awk -v k="$env_key" '$1==k":" {print $2}' "$CONFIG_FILE" | tr -d '[]')
    if [[ -n "$env_dir" && -d "$env_dir/bin" ]]; then PATH="$PATH:$env_dir/bin"; fi
done
export PATH

# GTDB-Tk reference data (path: gtdbtk in the config)
# (a [path/to/...] placeholder is ignored)
if [[ -z "${GTDBTK_DATA_PATH:-}" ]]; then
    GTDBTK_DATA_PATH=$(read_config_value path gtdbtk)
    [[ "$GTDBTK_DATA_PATH" == \[* ]] && GTDBTK_DATA_PATH=""
fi
export GTDBTK_DATA_PATH

if ! validate_env; then
    exit 1
fi

# The .smk files resolve a relative OUTPUT_DIR against the Scripts/ folder, so use absolute paths
# (snakemake runs from the Snakefile directory context).
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
SCRIPTS_DIR="$(cd "$SCRIPTS_DIR" && pwd)"

# The Snakefiles read raw data from <OUTPUT_DIR>/<DATA_FOLDER>/<sample>/<sample>_R{1,2}.fastq.gz.
# To keep the data where it is (no copy / symlink in the output directory), pass DATA_FOLDER
# as a path relative to OUTPUT_DIR (e.g. ../Toy_Dataset).
# (perl instead of 'realpath --relative-to', which macOS realpath does not support)
if [[ -d "$DATA_FOLDER" ]]; then
    data_abs="$(cd "$DATA_FOLDER" && pwd)"
elif [[ -d "$OUTPUT_DIR/$DATA_FOLDER" ]]; then
    data_abs="$(cd "$OUTPUT_DIR/$DATA_FOLDER" && pwd)"
else
    log_error "Data folder not found: $DATA_FOLDER"
    exit 1
fi
DATA_FOLDER="$(perl -MFile::Spec -e 'print File::Spec->abs2rel($ARGV[0], $ARGV[1])' "$data_abs" "$OUTPUT_DIR")"

# Check if reporter module is in MODULES and the metadata file is provided
if [[ " ${MODULES[@]} " =~ " reporter " ]] && [[ -z "$METADATA_FILE" ]]; then
    log_error "--reporter requires --metadata parameter"
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
