#!/usr/bin/env bash
#
# check_dependencies.sh - Audit system for pipeline dependencies
#
# Description:
#   This script checks for all tools, packages, and reference files required
#   by the Plasmidsaurus RNA-seq pipeline. Run this before attempting to
#   execute the pipeline to identify missing dependencies.
#
# Usage:
#   ./check_dependencies.sh [-v] [-o output_file]
#
# Arguments:
#   -v, --verbose       Show detailed version information
#   -o, --output        Save report to file (default: print to stdout)
#   -h, --help          Display this help message
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.0.0

set -uo pipefail

# ==============================================================================
# CONSTANTS
# ==============================================================================
readonly SCRIPT_NAME="$(basename "$0")"

# Color codes for output (disabled if not a terminal)
if [[ -t 1 ]]; then
    readonly RED='\033[0;31m'
    readonly GREEN='\033[0;32m'
    readonly YELLOW='\033[0;33m'
    readonly BLUE='\033[0;34m'
    readonly NC='\033[0m' # No Color
else
    readonly RED=''
    readonly GREEN=''
    readonly YELLOW=''
    readonly BLUE=''
    readonly NC=''
fi

# Counters
PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

print_header() {
    echo ""
    echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
}

print_subheader() {
    echo ""
    echo -e "${BLUE}--- $1 ---${NC}"
}

pass() {
    echo -e "  [${GREEN}✓${NC}] $1"
    ((PASS_COUNT++))
}

fail() {
    echo -e "  [${RED}✗${NC}] $1"
    ((FAIL_COUNT++))
}

warn() {
    echo -e "  [${YELLOW}!${NC}] $1"
    ((WARN_COUNT++))
}

info() {
    echo -e "  [${BLUE}i${NC}] $1"
}

usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} [-v] [-o output_file]

Check system for pipeline dependencies.

Options:
  -v, --verbose    Show detailed version information
  -o, --output     Save report to file
  -h, --help       Display this help message

EOF
    exit 0
}

# ==============================================================================
# CHECK FUNCTIONS
# ==============================================================================

#######################################
# Check if a command exists and optionally get version
#
# Arguments:
#   $1 - Command name
#   $2 - Version flag (optional, e.g., "--version")
#   $3 - Minimum version (optional, for display only)
#######################################
check_command() {
    local cmd="$1"
    local version_flag="${2:---version}"
    local min_version="${3:-}"

    if command -v "${cmd}" &> /dev/null; then
        local version_info
        version_info=$("${cmd}" ${version_flag} 2>&1 | head -3 | tr '\n' ' ' | cut -c1-60)

        if [[ -n "${min_version}" ]]; then
            pass "${cmd} (need >= ${min_version}): ${version_info}"
        else
            pass "${cmd}: ${version_info}"
        fi

        # Also show the path
        if [[ "${VERBOSE}" == true ]]; then
            info "  Path: $(which "${cmd}")"
        fi
        return 0
    else
        if [[ -n "${min_version}" ]]; then
            fail "${cmd} (need >= ${min_version}): NOT FOUND"
        else
            fail "${cmd}: NOT FOUND"
        fi
        return 1
    fi
}

#######################################
# Check for R package
#
# Arguments:
#   $1 - Package name
#   $2 - Minimum version (optional)
#######################################
check_r_package() {
    local package="$1"
    local min_version="${2:-}"

    if ! command -v Rscript &> /dev/null; then
        fail "R package ${package}: R not installed"
        return 1
    fi

    local result
    result=$(Rscript -e "if(requireNamespace('${package}', quietly=TRUE)) { cat('INSTALLED:', as.character(packageVersion('${package}'))) } else { cat('MISSING') }" 2>/dev/null)

    if [[ "${result}" == MISSING ]]; then
        fail "R::${package}${min_version:+ (>= ${min_version})}: NOT INSTALLED"
        return 1
    else
        local version="${result#INSTALLED: }"
        pass "R::${package}${min_version:+ (>= ${min_version})}: v${version}"
        return 0
    fi
}

#######################################
# Check for Python package
#
# Arguments:
#   $1 - Package name
#   $2 - Minimum version (optional)
#######################################
check_python_package() {
    local package="$1"
    local min_version="${2:-}"

    if ! command -v python3 &> /dev/null; then
        fail "Python package ${package}: Python3 not installed"
        return 1
    fi

    local result
    result=$(python3 -c "
try:
    import ${package}
    try:
        print('INSTALLED:', ${package}.__version__)
    except:
        print('INSTALLED: unknown version')
except ImportError:
    print('MISSING')
" 2>/dev/null)

    if [[ "${result}" == MISSING ]]; then
        fail "Python::${package}${min_version:+ (>= ${min_version})}: NOT INSTALLED"
        return 1
    else
        local version="${result#INSTALLED: }"
        pass "Python::${package}${min_version:+ (>= ${min_version})}: v${version}"
        return 0
    fi
}

#######################################
# Search for reference files in common locations
#
# Arguments:
#   $1 - Description
#   $2 - File pattern (glob)
#   $3... - Directories to search
#######################################
search_reference() {
    local description="$1"
    local pattern="$2"
    shift 2
    local dirs=("$@")

    local found_files=()

    for dir in "${dirs[@]}"; do
        if [[ -d "${dir}" ]]; then
            while IFS= read -r -d '' file; do
                found_files+=("${file}")
            done < <(find "${dir}" -maxdepth 3 -name "${pattern}" -type f -print0 2>/dev/null | head -z -5)
        fi
    done

    if [[ ${#found_files[@]} -gt 0 ]]; then
        pass "${description}: Found ${#found_files[@]} candidate(s)"
        for f in "${found_files[@]:0:3}"; do
            info "  ${f}"
        done
        if [[ ${#found_files[@]} -gt 3 ]]; then
            info "  ... and $((${#found_files[@]} - 3)) more"
        fi
        return 0
    else
        warn "${description}: Not found in common locations"
        return 1
    fi
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

VERBOSE=false
OUTPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# ==============================================================================
# MAIN CHECKS
# ==============================================================================

# Redirect output if file specified
if [[ -n "${OUTPUT_FILE}" ]]; then
    exec > >(tee "${OUTPUT_FILE}")
fi

echo ""
echo -e "${BLUE}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║     PLASMIDSAURUS RNA-SEQ PIPELINE DEPENDENCY CHECK          ║${NC}"
echo -e "${BLUE}║     $(date '+%Y-%m-%d %H:%M:%S')                                      ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════════════════════════════╝${NC}"

# ------------------------------------------------------------------------------
# System Information
# ------------------------------------------------------------------------------
print_header "SYSTEM INFORMATION"

info "Hostname: $(hostname)"
info "User: $(whoami)"
info "Working directory: $(pwd)"
info "Shell: ${SHELL}"

if [[ -f /etc/os-release ]]; then
    source /etc/os-release
    info "OS: ${PRETTY_NAME:-${NAME} ${VERSION}}"
fi

info "Kernel: $(uname -r)"
info "Architecture: $(uname -m)"

# Check available memory
if command -v free &> /dev/null; then
    total_mem=$(free -g | awk '/^Mem:/{print $2}')
    info "Total RAM: ${total_mem} GB"
fi

# Check available cores
n_cores=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "unknown")
info "CPU cores: ${n_cores}"

# ------------------------------------------------------------------------------
# Core Command-Line Tools
# ------------------------------------------------------------------------------
print_header "CORE COMMAND-LINE TOOLS"

print_subheader "Read Processing"
check_command "fastp" "--version" "0.24.0"

print_subheader "Alignment"
check_command "STAR" "--version" "2.7.11"

print_subheader "BAM Processing"
check_command "samtools" "--version" "1.22"

print_subheader "UMI Deduplication"
# UMICollapse can be installed multiple ways
if command -v umicollapse &> /dev/null; then
    pass "umicollapse: $(umicollapse --version 2>&1 | head -1 || echo 'installed')"
elif command -v java &> /dev/null; then
    warn "umicollapse command not found, but Java is available"
    info "  Check if UMICOLLAPSE_JAR environment variable is set: ${UMICOLLAPSE_JAR:-NOT SET}"
    # Search for JAR file
    for jar_path in \
        "${HOME}/tools/umicollapse"/*.jar \
        "/opt/umicollapse"/*.jar \
        "/usr/local/share/umicollapse"/*.jar; do
        if [[ -f "${jar_path}" ]]; then
            info "  Found JAR: ${jar_path}"
        fi
    done
else
    fail "umicollapse: NOT FOUND (and no Java available)"
fi

print_subheader "Quality Control"
check_command "fastqc" "--version" || true
check_command "multiqc" "--version" "1.32"

# RSeQC tools (Python-based)
for tool in infer_experiment.py read_distribution.py bam_stat.py geneBody_coverage.py; do
    check_command "${tool}" "--version" || true
done

check_command "qualimap" "--help" "2.3" || true

print_subheader "Quantification"
check_command "featureCounts" "-v" "2.1"

# ------------------------------------------------------------------------------
# R Environment
# ------------------------------------------------------------------------------
print_header "R ENVIRONMENT"

if check_command "R" "--version"; then
    check_command "Rscript" "--version"

    print_subheader "Required R Packages"
    check_r_package "edgeR" "4.0"
    check_r_package "ggplot2"
    check_r_package "pheatmap"
    check_r_package "RColorBrewer"
    check_r_package "ggrepel"
    check_r_package "optparse"

    print_subheader "Optional R Packages"
    check_r_package "DESeq2" || true
    check_r_package "limma" || true
fi

# ------------------------------------------------------------------------------
# Python Environment
# ------------------------------------------------------------------------------
print_header "PYTHON ENVIRONMENT"

if check_command "python3" "--version"; then

    print_subheader "Required Python Packages"
    check_python_package "gseapy" "0.12"
    check_python_package "pandas"
    check_python_package "matplotlib"
    check_python_package "numpy"

    print_subheader "Optional Python Packages"
    check_python_package "scipy" || true
    check_python_package "seaborn" || true
fi

# ------------------------------------------------------------------------------
# Reference File Locations
# ------------------------------------------------------------------------------
print_header "REFERENCE FILE SEARCH"

info "Searching common locations for reference files..."
info "(This may take a moment)"

# Common reference genome locations
ref_dirs=(
    "/data/references"
    "/data/ref"
    "/data/genomes"
    "/references"
    "/ref"
    "/home/*/references"
    "/home/*/data/references"
    "${HOME}/references"
    "${HOME}/data/references"
    "/opt/references"
    "/usr/local/share/references"
    "/scratch/references"
    "/work/references"
)

# Also check for any environment variables pointing to references
if [[ -n "${GENOME_DIR:-}" ]]; then
    ref_dirs+=("${GENOME_DIR}")
    info "GENOME_DIR env variable: ${GENOME_DIR}"
fi

if [[ -n "${STAR_INDEX:-}" ]]; then
    info "STAR_INDEX env variable: ${STAR_INDEX}"
fi

print_subheader "Genome FASTA Files"
search_reference "Human genome (hg38/GRCh38)" "*.fa" "${ref_dirs[@]}" /data/**/hg38* /data/**/GRCh38* 2>/dev/null || true
search_reference "Mouse genome (mm10/GRCm38)" "*.fa" "${ref_dirs[@]}" /data/**/mm10* /data/**/GRCm38* 2>/dev/null || true

print_subheader "GTF Annotation Files"
search_reference "GTF annotations" "*.gtf" "${ref_dirs[@]}" /data/**/annotation* 2>/dev/null || true

print_subheader "STAR Index Directories"
# Look for STAR index (contains SA file)
for dir in "${ref_dirs[@]}" /data/**/star* /data/**/STAR*; do
    if [[ -d "${dir}" ]] && [[ -f "${dir}/SA" ]]; then
        pass "STAR index found: ${dir}"
    fi
done 2>/dev/null || warn "No STAR indices found in common locations"

print_subheader "BED Gene Model Files"
search_reference "BED gene models" "*.bed" "${ref_dirs[@]}" 2>/dev/null || true

# ------------------------------------------------------------------------------
# Environment Variables
# ------------------------------------------------------------------------------
print_header "ENVIRONMENT VARIABLES"

env_vars=(
    "PATH"
    "CONDA_PREFIX"
    "CONDA_DEFAULT_ENV"
    "VIRTUAL_ENV"
    "R_HOME"
    "R_LIBS"
    "R_LIBS_USER"
    "PYTHONPATH"
    "STAR_INDEX"
    "GENOME_DIR"
    "UMICOLLAPSE_JAR"
)

for var in "${env_vars[@]}"; do
    value="${!var:-}"
    if [[ -n "${value}" ]]; then
        # Truncate long paths
        if [[ ${#value} -gt 60 ]]; then
            info "${var}: ${value:0:57}..."
        else
            info "${var}: ${value}"
        fi
    fi
done

# Check if running in conda environment
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    info "Running in conda environment: ${CONDA_DEFAULT_ENV:-unknown}"
fi

# ------------------------------------------------------------------------------
# Module System (HPC)
# ------------------------------------------------------------------------------
print_header "HPC MODULE SYSTEM"

if command -v module &> /dev/null; then
    pass "Module system available"

    info "Currently loaded modules:"
    module list 2>&1 | head -20 | while read -r line; do
        info "  ${line}"
    done

    print_subheader "Available Relevant Modules (sampling)"
    for pattern in star samtools fastp R python; do
        avail=$(module avail "${pattern}" 2>&1 | head -5)
        if [[ -n "${avail}" ]]; then
            info "${pattern}: modules available"
        fi
    done
else
    info "No module system detected (not an HPC environment or modules not configured)"
fi

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
print_header "SUMMARY"

echo ""
echo -e "  ${GREEN}Passed:${NC}  ${PASS_COUNT}"
echo -e "  ${RED}Failed:${NC}  ${FAIL_COUNT}"
echo -e "  ${YELLOW}Warnings:${NC} ${WARN_COUNT}"
echo ""

if [[ ${FAIL_COUNT} -eq 0 ]]; then
    echo -e "${GREEN}All required dependencies found!${NC}"
    echo "The pipeline should be ready to run (after configuring reference paths)."
    exit 0
elif [[ ${FAIL_COUNT} -lt 5 ]]; then
    echo -e "${YELLOW}Some dependencies missing.${NC}"
    echo "Install missing tools before running the pipeline."
    exit 1
else
    echo -e "${RED}Multiple dependencies missing.${NC}"
    echo "Significant setup required before running the pipeline."
    exit 2
fi
