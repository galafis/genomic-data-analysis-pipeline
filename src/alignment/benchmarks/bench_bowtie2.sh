#!/bin/bash

#=============================================================================
# Bowtie2 Alignment Benchmarking Script
# 
# Description: Professional benchmarking script for Bowtie2 aligner
# Author: Genomic Data Analysis Pipeline
# Created: $(date '+%Y-%m-%d')
# Version: 1.0
#=============================================================================

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Script configuration
SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
RESULTS_DIR="${SCRIPT_DIR}/results"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

# Default parameters
DEFAULT_THREADS=8
DEFAULT_MEMORY="8G"
DEFAULT_PRESET="--sensitive"
VERBOSE=false

#=============================================================================
# FUNCTIONS
#=============================================================================

# Display usage information
show_usage() {
    cat << EOF
Usage: $SCRIPT_NAME [OPTIONS] -r REFERENCE -1 READ1 -2 READ2 -o OUTPUT

Bowtie2 Alignment Benchmarking Script

Required Arguments:
  -r, --reference FILE    Reference genome FASTA file
  -1, --read1 FILE        First paired-end FASTQ file
  -2, --read2 FILE        Second paired-end FASTQ file
  -o, --output FILE       Output BAM file

Optional Arguments:
  -t, --threads INT       Number of threads (default: $DEFAULT_THREADS)
  -m, --memory STR        Memory limit (default: $DEFAULT_MEMORY)
  -p, --preset STR        Bowtie2 preset (default: $DEFAULT_PRESET)
  -v, --verbose           Enable verbose output
  -h, --help              Show this help message

Presets:
  --very-fast-local       Very fast local alignment
  --fast-local            Fast local alignment
  --sensitive-local       Sensitive local alignment
  --very-sensitive-local  Very sensitive local alignment
  --very-fast             Very fast end-to-end alignment
  --fast                  Fast end-to-end alignment
  --sensitive             Sensitive end-to-end alignment (default)
  --very-sensitive        Very sensitive end-to-end alignment

Example:
  $SCRIPT_NAME -r genome.fa -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o aligned.bam
  $SCRIPT_NAME -r genome.fa -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o aligned.bam -t 16 -p --very-sensitive

Output Files:
  - OUTPUT.bam              Aligned reads in BAM format
  - OUTPUT.bam.bai          BAM index file
  - OUTPUT_stats.txt        Alignment statistics
  - OUTPUT_benchmark.log    Detailed benchmark log
  - OUTPUT_time_memory.log  Time and memory usage log

EOF
}

# Logging function
log_message() {
    local level="$1"
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
}

# Error handling function
error_exit() {
    log_message "ERROR" "$1"
    exit 1
}

# Check if required tools are available
check_dependencies() {
    local tools=("bowtie2" "bowtie2-build" "samtools" "/usr/bin/time")
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error_exit "Required tool not found: $tool"
        fi
    done
    
    log_message "INFO" "All dependencies verified"
}

# Validate input files
validate_inputs() {
    local files=("$REFERENCE" "$READ1" "$READ2")
    
    for file in "${files[@]}"; do
        if [[ ! -f "$file" ]]; then
            error_exit "Input file not found: $file"
        fi
        if [[ ! -r "$file" ]]; then
            error_exit "Input file not readable: $file"
        fi
    done
    
    # Check if reference index exists
    if [[ ! -f "${REFERENCE}.1.bt2" ]]; then
        log_message "WARN" "Bowtie2 index not found. Will build index."
        BUILD_INDEX=true
    else
        log_message "INFO" "Bowtie2 index found"
        BUILD_INDEX=false
    fi
    
    log_message "INFO" "Input validation completed"
}

# Build Bowtie2 index if needed
build_index() {
    if [[ "$BUILD_INDEX" == "true" ]]; then
        log_message "INFO" "Building Bowtie2 index for $REFERENCE"
        
        /usr/bin/time -v bowtie2-build \
            "$REFERENCE" \
            "$REFERENCE" \
            --threads "$THREADS" \
            2>&1 | tee -a "$TIME_LOG"
        
        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            error_exit "Failed to build Bowtie2 index"
        fi
        
        log_message "INFO" "Bowtie2 index building completed"
    fi
}

# Perform alignment with benchmarking
run_alignment() {
    log_message "INFO" "Starting Bowtie2 alignment"
    log_message "INFO" "Reference: $REFERENCE"
    log_message "INFO" "Read 1: $READ1"
    log_message "INFO" "Read 2: $READ2"
    log_message "INFO" "Output: $OUTPUT"
    log_message "INFO" "Threads: $THREADS"
    log_message "INFO" "Preset: $PRESET"
    
    # Create output directory if needed
    mkdir -p "$(dirname "$OUTPUT")"
    
    # Run Bowtie2 alignment with time and memory monitoring
    /usr/bin/time -v bowtie2 \
        $PRESET \
        -x "$REFERENCE" \
        -1 "$READ1" \
        -2 "$READ2" \
        -p "$THREADS" \
        --rg-id "$SAMPLE_ID" \
        --rg "SM:$SAMPLE_ID" \
        --rg "PL:ILLUMINA" \
        --rg "LB:$SAMPLE_ID" \
        2> "${OUTPUT%.bam}_bowtie2.log" | \
    samtools view -@ "$((THREADS-1))" -bS - | \
    samtools sort -@ "$((THREADS-1))" -m "$MEMORY" -o "$OUTPUT" - \
    2>&1 | tee -a "$TIME_LOG"
    
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        error_exit "Bowtie2 alignment failed"
    fi
    
    log_message "INFO" "Alignment completed successfully"
}

# Generate alignment statistics
generate_stats() {
    log_message "INFO" "Generating alignment statistics"
    
    # Index BAM file
    samtools index "$OUTPUT"
    
    # Generate comprehensive statistics
    {
        echo "=== Bowtie2 Alignment Statistics ==="
        echo "Date: $(date)"
        echo "Reference: $REFERENCE"
        echo "Sample: $SAMPLE_ID"
        echo "Output: $OUTPUT"
        echo ""
        
        echo "=== Bowtie2 Alignment Summary ==="
        grep -E "(reads processed|aligned|paired|singleton|discordant)" "${OUTPUT%.bam}_bowtie2.log" || true
        echo ""
        
        echo "=== SAMtools Statistics ==="
        samtools stats "$OUTPUT" | head -20
        echo ""
        
        echo "=== SAMtools Flagstat ==="
        samtools flagstat "$OUTPUT"
        echo ""
        
        echo "=== SAMtools Index Stats ==="
        samtools idxstats "$OUTPUT" | head -10
        echo ""
        
        echo "=== File Sizes ==="
        ls -lh "$OUTPUT"* | awk '{print $5 "\t" $9}'
        
    } > "$STATS_FILE"
    
    log_message "INFO" "Statistics generated: $STATS_FILE"
}

# Display benchmark summary
show_summary() {
    log_message "INFO" "Benchmark Summary"
    log_message "INFO" "=================="
    
    # Extract key metrics
    local total_reads=$(grep -o "[0-9]\+ reads processed" "${OUTPUT%.bam}_bowtie2.log" | cut -d' ' -f1 || echo "N/A")
    local aligned_reads=$(grep -o "[0-9]\+ ([0-9.]\+%) aligned" "${OUTPUT%.bam}_bowtie2.log" | head -1 | cut -d' ' -f1 || echo "N/A")
    local runtime=$(grep "Elapsed (wall clock) time" "$TIME_LOG" | tail -1 | cut -d' ' -f8 || echo "N/A")
    local max_memory=$(grep "Maximum resident set size" "$TIME_LOG" | tail -1 | awk '{print $6/1024 "MB"}' || echo "N/A")
    
    log_message "INFO" "Total reads processed: $total_reads"
    log_message "INFO" "Aligned reads: $aligned_reads"
    log_message "INFO" "Runtime: $runtime"
    log_message "INFO" "Peak memory: $max_memory"
    log_message "INFO" "Output BAM: $OUTPUT"
    log_message "INFO" "Statistics: $STATS_FILE"
    log_message "INFO" "Logs: $LOG_FILE"
}

#=============================================================================
# ARGUMENT PARSING
#=============================================================================

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -1|--read1)
            READ1="$2"
            shift 2
            ;;
        -2|--read2)
            READ2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -p|--preset)
            PRESET="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            show_usage
            exit 1
            ;;
    esac
done

#=============================================================================
# MAIN EXECUTION
#=============================================================================

# Check required arguments
if [[ -z "${REFERENCE:-}" || -z "${READ1:-}" || -z "${READ2:-}" || -z "${OUTPUT:-}" ]]; then
    echo "Error: Missing required arguments" >&2
    show_usage
    exit 1
fi

# Set defaults
THREADS="${THREADS:-$DEFAULT_THREADS}"
MEMORY="${MEMORY:-$DEFAULT_MEMORY}"
PRESET="${PRESET:-$DEFAULT_PRESET}"

# Generate sample ID from output filename
SAMPLE_ID="$(basename "${OUTPUT%.bam}")"

# Setup output files
mkdir -p "$LOG_DIR" "$RESULTS_DIR"
LOG_FILE="$LOG_DIR/${SAMPLE_ID}_${TIMESTAMP}_benchmark.log"
TIME_LOG="$LOG_DIR/${SAMPLE_ID}_${TIMESTAMP}_time_memory.log"
STATS_FILE="$RESULTS_DIR/${SAMPLE_ID}_${TIMESTAMP}_stats.txt"

# Start logging
log_message "INFO" "Starting Bowtie2 benchmark for $SAMPLE_ID"
log_message "INFO" "Script: $SCRIPT_NAME"
log_message "INFO" "Command: $0 $*"

# Main execution flow
check_dependencies
validate_inputs
build_index
run_alignment
generate_stats
show_summary

log_message "INFO" "Bowtie2 benchmark completed successfully"
exit 0
