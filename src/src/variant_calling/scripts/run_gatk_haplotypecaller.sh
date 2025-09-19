#!/bin/bash

# ==============================================================================
# SCRIPT: run_gatk_haplotypecaller.sh
# AUTHOR: Genomic Data Analysis Pipeline Team
# VERSION: 1.0
# DATE: 2025-09-19
# DESCRIPTION: Executes GATK HaplotypeCaller for variant calling from BAM files
# USAGE: ./run_gatk_haplotypecaller.sh -b <BAM> -r <REF> -o <OUT> [-t <THREADS>] [-c <CONFIG>]
# ==============================================================================

# Exit on any error
set -euo pipefail

# Default values
THREADS=8
CONFIG_FILE=""
LOG_DIR="./logs"
LOG_FILE="$LOG_DIR/gatk_haplotypecaller_$(date +%Y%m%d_%H%M%S).log"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

# Function to log messages with timestamp
log_message() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
}

# Function to log info messages
log_info() {
    log_message "INFO" "$1"
}

# Function to log warning messages
log_warn() {
    log_message "WARN" "$1"
}

# Function to log error messages
log_error() {
    log_message "ERROR" "$1"
}

# ==============================================================================
# ERROR HANDLING
# ==============================================================================

# Function to handle script termination
cleanup() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log_error "Script terminated with exit code: $exit_code"
        log_error "Check log file for details: $LOG_FILE"
    else
        log_info "Script completed successfully"
    fi
    exit $exit_code
}

# Trap to handle unexpected exits
trap cleanup EXIT

# ==============================================================================
# USAGE FUNCTION
# ==============================================================================

usage() {
    cat << EOF
USAGE: $0 -b <BAM_FILE> -r <REFERENCE_GENOME> -o <OUTPUT_VCF> [-t <THREADS>] [-c <CONFIG_FILE>]

DESCRIPTION:
    Executes GATK HaplotypeCaller for variant calling from aligned BAM files.
    Performs quality control checks and generates comprehensive logging.

REQUIRED ARGUMENTS:
    -b, --bam       Input BAM file (sorted and indexed)
    -r, --reference Reference genome FASTA file (with .fai index)
    -o, --output    Output VCF file path

OPTIONAL ARGUMENTS:
    -t, --threads   Number of threads to use (default: 8)
    -c, --config    Configuration file with additional GATK parameters
    -h, --help      Show this help message

EXAMPLE:
    $0 -b sample.sorted.bam -r hg38.fa -o variants.vcf -t 16
    $0 --bam sample.bam --reference genome.fa --output out.vcf --config gatk.conf

NOTE:
    - BAM file must be sorted and indexed (.bai file required)
    - Reference genome must have corresponding .fai index file
    - Output directory must be writable
    - GATK must be available in PATH or GATK_PATH environment variable

EOF
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

# Initialize variables
BAM_FILE=""
REF_GENOME=""
OUT_VCF=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILE="$2"
            shift 2
            ;;
        -r|--reference)
            REF_GENOME="$2"
            shift 2
            ;;
        -o|--output)
            OUT_VCF="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_info "Starting GATK HaplotypeCaller pipeline"
log_info "Log file: $LOG_FILE"

# Check required arguments
if [[ -z "$BAM_FILE" || -z "$REF_GENOME" || -z "$OUT_VCF" ]]; then
    log_error "Missing required arguments. BAM, reference, and output files are mandatory."
    usage
    exit 1
fi

# Validate input files exist
if [[ ! -f "$BAM_FILE" ]]; then
    log_error "BAM file not found: $BAM_FILE"
    exit 1
fi

if [[ ! -f "$REF_GENOME" ]]; then
    log_error "Reference genome file not found: $REF_GENOME"
    exit 1
fi

# Check for BAM index
if [[ ! -f "${BAM_FILE}.bai" && ! -f "${BAM_FILE%.*}.bai" ]]; then
    log_error "BAM index file not found. Please index your BAM file with: samtools index $BAM_FILE"
    exit 1
fi

# Check for reference index
if [[ ! -f "${REF_GENOME}.fai" ]]; then
    log_error "Reference index file not found. Please index with: samtools faidx $REF_GENOME"
    exit 1
fi

# Validate threads parameter
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    log_error "Invalid threads value: $THREADS. Must be a positive integer."
    exit 1
fi

# Check if GATK is available
if ! command -v gatk &> /dev/null && [[ -z "${GATK_PATH:-}" ]]; then
    log_error "GATK not found in PATH. Please install GATK or set GATK_PATH environment variable."
    exit 1
fi

# Set GATK command
GATK_CMD="gatk"
if [[ -n "${GATK_PATH:-}" ]]; then
    GATK_CMD="$GATK_PATH/gatk"
fi

# Create output directory if it doesn't exist
OUT_DIR=$(dirname "$OUT_VCF")
mkdir -p "$OUT_DIR"

log_info "Input validation completed successfully"
log_info "BAM file: $BAM_FILE"
log_info "Reference: $REF_GENOME"
log_info "Output VCF: $OUT_VCF"
log_info "Threads: $THREADS"
[[ -n "$CONFIG_FILE" ]] && log_info "Config file: $CONFIG_FILE"

# ==============================================================================
# CONFIGURATION FILE PROCESSING
# ==============================================================================

EXTRA_ARGS=""
if [[ -n "$CONFIG_FILE" ]]; then
    if [[ -f "$CONFIG_FILE" ]]; then
        log_info "Loading configuration from: $CONFIG_FILE"
        # Read additional arguments from config file (one per line, ignore comments)
        while IFS= read -r line || [[ -n "$line" ]]; do
            # Skip empty lines and comments
            [[ "$line" =~ ^[[:space:]]*$ ]] && continue
            [[ "$line" =~ ^[[:space:]]*# ]] && continue
            EXTRA_ARGS="$EXTRA_ARGS $line"
        done < "$CONFIG_FILE"
        log_info "Additional arguments from config: $EXTRA_ARGS"
    else
        log_warn "Configuration file not found: $CONFIG_FILE. Proceeding without additional parameters."
    fi
fi

# ==============================================================================
# GATK HAPLOTYPECALLER EXECUTION
# ==============================================================================

log_info "Starting GATK HaplotypeCaller execution"

# Build GATK command with comprehensive parameters
GATK_COMMAND="$GATK_CMD HaplotypeCaller \\
    --input $BAM_FILE \\
    --reference $REF_GENOME \\
    --output $OUT_VCF \\
    --native-pair-hmm-threads $THREADS \\
    --verbosity INFO \\
    --create-output-variant-index true \\
    --annotation AlleleFraction \\
    --annotation Coverage \\
    --annotation FisherStrand \\
    --annotation MappingQualityRankSumTest \\
    --annotation QualByDepth \\
    --annotation RMSMappingQuality \\
    --annotation ReadPosRankSumTest \\
    --annotation StrandOddsRatio$EXTRA_ARGS"

log_info "Executing GATK command:"
log_info "$GATK_COMMAND"

# Execute GATK HaplotypeCaller with error handling
start_time=$(date +%s)

if eval "$GATK_COMMAND" 2>&1 | tee -a "$LOG_FILE"; then
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    log_info "GATK HaplotypeCaller completed successfully in ${duration} seconds"
    
    # Verify output file was created and has content
    if [[ -f "$OUT_VCF" && -s "$OUT_VCF" ]]; then
        variant_count=$(grep -v "^#" "$OUT_VCF" | wc -l)
        log_info "Output VCF created: $OUT_VCF"
        log_info "Number of variants called: $variant_count"
        
        # Log file sizes for reference
        bam_size=$(du -h "$BAM_FILE" | cut -f1)
        vcf_size=$(du -h "$OUT_VCF" | cut -f1)
        log_info "Input BAM size: $bam_size"
        log_info "Output VCF size: $vcf_size"
    else
        log_error "Output VCF file was not created or is empty: $OUT_VCF"
        exit 1
    fi
else
    log_error "GATK HaplotypeCaller failed. Check log file for details: $LOG_FILE"
    exit 1
fi

# ==============================================================================
# POST-PROCESSING AND SUMMARY
# ==============================================================================

log_info "Generating execution summary"

# Create summary file
SUMMARY_FILE="${OUT_VCF%.vcf}_summary.txt"
cat > "$SUMMARY_FILE" << EOF
GATK HaplotypeCaller Execution Summary
======================================
Date: $(date)
Script: $0
BAM File: $BAM_FILE
Reference: $REF_GENOME
Output VCF: $OUT_VCF
Threads Used: $THREADS
Execution Time: ${duration} seconds
Variants Called: $variant_count
Log File: $LOG_FILE

Command Executed:
$GATK_COMMAND
EOF

log_info "Summary written to: $SUMMARY_FILE"
log_info "Pipeline completed successfully. Output files:"
log_info "  - VCF: $OUT_VCF"
log_info "  - Index: ${OUT_VCF}.idx"
log_info "  - Summary: $SUMMARY_FILE"
log_info "  - Log: $LOG_FILE"

exit 0
