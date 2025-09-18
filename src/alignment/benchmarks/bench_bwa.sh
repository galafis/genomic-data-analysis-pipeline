#!/bin/bash

# bench_bwa.sh - BWA-MEM Alignment Benchmarking Script
# Description: Benchmarks BWA-MEM performance with configurable parameters
# Author: Genomic Data Analysis Pipeline
# Version: 1.0
# Date: $(date '+%Y-%m-%d')

set -euo pipefail

# Default configuration
DEFAULT_THREADS=4
DEFAULT_OUTPUT_DIR="./benchmark_results"
TIME_LOG="benchmark_time.log"
METRICS_LOG="alignment_metrics.log"

# Function to display usage
show_usage() {
    cat << EOF
Usage: $0 -r REFERENCE -1 FASTQ1 [-2 FASTQ2] -o OUTPUT_PREFIX [OPTIONS]

Required arguments:
  -r, --reference    Reference genome FASTA file
  -1, --fastq1      Input FASTQ file (read 1)
  -o, --output      Output BAM file prefix

Optional arguments:
  -2, --fastq2      Input FASTQ file (read 2, for paired-end)
  -t, --threads     Number of threads (default: $DEFAULT_THREADS)
  -d, --outdir      Output directory (default: $DEFAULT_OUTPUT_DIR)
  -h, --help        Show this help message

Example:
  # Single-end reads
  $0 -r genome.fa -1 sample.fastq -o sample_aligned
  
  # Paired-end reads with 8 threads
  $0 -r genome.fa -1 sample_R1.fastq -2 sample_R2.fastq -o sample_aligned -t 8
EOF
}

# Function to log messages with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$TIME_LOG"
}

# Function to check if required tools are available
check_dependencies() {
    local tools=("bwa" "samtools")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "Error: $tool is not installed or not in PATH" >&2
            exit 1
        fi
    done
}

# Function to validate input files
validate_inputs() {
    if [[ ! -f "$REFERENCE" ]]; then
        echo "Error: Reference file '$REFERENCE' not found" >&2
        exit 1
    fi
    
    if [[ ! -f "$FASTQ1" ]]; then
        echo "Error: FASTQ1 file '$FASTQ1' not found" >&2
        exit 1
    fi
    
    if [[ -n "$FASTQ2" && ! -f "$FASTQ2" ]]; then
        echo "Error: FASTQ2 file '$FASTQ2' not found" >&2
        exit 1
    fi
}

# Function to run BWA-MEM alignment with benchmarking
run_bwa_benchmark() {
    local start_time=$(date +%s)
    local output_bam="${OUTPUT_DIR}/${OUTPUT_PREFIX}.bam"
    local sorted_bam="${OUTPUT_DIR}/${OUTPUT_PREFIX}_sorted.bam"
    
    log_message "Starting BWA-MEM alignment benchmark"
    log_message "Reference: $REFERENCE"
    log_message "FASTQ1: $FASTQ1"
    [[ -n "$FASTQ2" ]] && log_message "FASTQ2: $FASTQ2"
    log_message "Threads: $THREADS"
    log_message "Output: $output_bam"
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Run BWA-MEM alignment
    log_message "Running BWA-MEM alignment..."
    if [[ -n "$FASTQ2" ]]; then
        # Paired-end alignment
        /usr/bin/time -v bwa mem -t "$THREADS" "$REFERENCE" "$FASTQ1" "$FASTQ2" 2>> "$TIME_LOG" | \
        samtools view -@ "$THREADS" -bS - > "$output_bam" 2>> "$TIME_LOG"
    else
        # Single-end alignment
        /usr/bin/time -v bwa mem -t "$THREADS" "$REFERENCE" "$FASTQ1" 2>> "$TIME_LOG" | \
        samtools view -@ "$THREADS" -bS - > "$output_bam" 2>> "$TIME_LOG"
    fi
    
    # Sort BAM file
    log_message "Sorting BAM file..."
    /usr/bin/time -v samtools sort -@ "$THREADS" -o "$sorted_bam" "$output_bam" 2>> "$TIME_LOG"
    
    # Index sorted BAM
    log_message "Indexing BAM file..."
    samtools index "$sorted_bam"
    
    # Generate alignment statistics
    log_message "Generating alignment metrics..."
    samtools flagstat "$sorted_bam" > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_flagstat.txt"
    samtools stats "$sorted_bam" > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_stats.txt"
    
    # Calculate execution time
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_message "BWA-MEM alignment completed successfully"
    log_message "Total execution time: ${duration} seconds"
    
    # Display summary
    echo ""
    echo "=== BENCHMARK SUMMARY ==="
    echo "Reference: $REFERENCE"
    echo "Input reads: $FASTQ1$([ -n "$FASTQ2" ] && echo ", $FASTQ2")"
    echo "Output BAM: $sorted_bam"
    echo "Threads used: $THREADS"
    echo "Execution time: ${duration}s"
    echo "Results directory: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"*
    echo ""
    echo "Quick alignment stats:"
    head -n 5 "${OUTPUT_DIR}/${OUTPUT_PREFIX}_flagstat.txt"
}

# Parse command line arguments
REFERENCE=""
FASTQ1=""
FASTQ2=""
OUTPUT_PREFIX=""
THREADS="$DEFAULT_THREADS"
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -1|--fastq1)
            FASTQ1="$2"
            shift 2
            ;;
        -2|--fastq2)
            FASTQ2="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -d|--outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1" >&2
            show_usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$REFERENCE" || -z "$FASTQ1" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: Missing required arguments" >&2
    show_usage
    exit 1
fi

# Main execution
check_dependencies
validate_inputs
run_bwa_benchmark

# Example simulated output:
# [2024-09-18 15:30:01] Starting BWA-MEM alignment benchmark
# [2024-09-18 15:30:01] Reference: /data/reference/hg38.fa
# [2024-09-18 15:30:01] FASTQ1: /data/samples/sample001_R1.fastq.gz
# [2024-09-18 15:30:01] FASTQ2: /data/samples/sample001_R2.fastq.gz
# [2024-09-18 15:30:01] Threads: 8
# [2024-09-18 15:30:01] Output: ./benchmark_results/sample001.bam
# [2024-09-18 15:30:01] Running BWA-MEM alignment...
# [2024-09-18 15:32:45] Sorting BAM file...
# [2024-09-18 15:33:12] Indexing BAM file...
# [2024-09-18 15:33:15] Generating alignment metrics...
# [2024-09-18 15:33:18] BWA-MEM alignment completed successfully
# [2024-09-18 15:33:18] Total execution time: 197 seconds
#
# === BENCHMARK SUMMARY ===
# Reference: /data/reference/hg38.fa
# Input reads: /data/samples/sample001_R1.fastq.gz, /data/samples/sample001_R2.fastq.gz
# Output BAM: ./benchmark_results/sample001_sorted.bam
# Threads used: 8
# Execution time: 197s
# Results directory: ./benchmark_results
#
# Generated files:
# -rw-r--r-- 1 user user 1.2G Sep 18 15:33 sample001.bam
# -rw-r--r-- 1 user user 892M Sep 18 15:33 sample001_sorted.bam
# -rw-r--r-- 1 user user  256 Sep 18 15:33 sample001_sorted.bam.bai
# -rw-r--r-- 1 user user 1.1K Sep 18 15:33 sample001_flagstat.txt
# -rw-r--r-- 1 user user  15K Sep 18 15:33 sample001_stats.txt
#
# Quick alignment stats:
# 2000000 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 18459 + 0 supplementary
# 0 + 0 duplicates
# 1945231 + 0 mapped (97.26% : N/A)
