#!/bin/bash

# bench_minimap2.sh - Minimap2 Alignment Benchmarking Script
# Part of the Genomic Data Analysis Pipeline
# Author: Genomic Data Analysis Pipeline Team
# Date: $(date +"%Y-%m-%d")
# Version: 1.0

# Description:
# Professional benchmarking script for Minimap2 long-read alignment tool
# Measures performance metrics including runtime, memory usage, and alignment statistics

set -euo pipefail

# Default configuration
DEFAULT_PRESET="map-ont"  # Oxford Nanopore preset
DEFAULT_THREADS=8
DEFAULT_OUTPUT_DIR="./benchmark_results"
LOG_LEVEL="INFO"

# Function to display usage information
usage() {
    cat << EOF
Usage: $0 -i <input_fastq> -r <reference> -o <output_bam> [OPTIONS]

Required arguments:
  -i, --input         Input FASTQ file (long reads)
  -r, --reference     Reference genome file (FASTA)
  -o, --output        Output BAM file

Optional arguments:
  -p, --preset        Minimap2 preset (default: ${DEFAULT_PRESET})
                      Options: map-ont, map-pb, map-hifi, asm5, asm10, asm20
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -d, --output-dir    Output directory for benchmark results (default: ${DEFAULT_OUTPUT_DIR})
  -l, --log-level     Logging level: DEBUG, INFO, WARN, ERROR (default: ${LOG_LEVEL})
  -h, --help          Display this help message

Example usage:
  $0 -i reads.fastq -r genome.fa -o aligned.bam -p map-ont -t 16
  $0 -i nanopore_reads.fastq.gz -r hg38.fa -o output.bam --preset map-ont --threads 32

Preset descriptions:
  map-ont    : Oxford Nanopore reads
  map-pb     : PacBio CLR reads
  map-hifi   : PacBio HiFi/CCS reads
  asm5       : Long assembly contigs (5% divergence)
  asm10      : Long assembly contigs (10% divergence)
  asm20      : Long assembly contigs (20% divergence)
EOF
}

# Logging function
log() {
    local level=$1
    shift
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] $*" >&2
}

# Function to check if required tools are installed
check_dependencies() {
    local missing_tools=()
    
    for tool in minimap2 samtools; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log "ERROR" "Missing required tools: ${missing_tools[*]}"
        log "INFO" "Please install missing dependencies and try again"
        exit 1
    fi
}

# Function to validate input files
validate_inputs() {
    local input_fastq="$1"
    local reference="$2"
    
    if [[ ! -f "$input_fastq" ]]; then
        log "ERROR" "Input FASTQ file not found: $input_fastq"
        exit 1
    fi
    
    if [[ ! -f "$reference" ]]; then
        log "ERROR" "Reference file not found: $reference"
        exit 1
    fi
    
    log "INFO" "Input validation completed successfully"
}

# Function to run Minimap2 with benchmarking
run_minimap2_benchmark() {
    local input_fastq="$1"
    local reference="$2"
    local output_bam="$3"
    local preset="$4"
    local threads="$5"
    local output_dir="$6"
    
    local benchmark_file="${output_dir}/minimap2_benchmark.txt"
    local log_file="${output_dir}/minimap2.log"
    local temp_sam="${output_dir}/temp_alignment.sam"
    
    log "INFO" "Starting Minimap2 alignment benchmark"
    log "INFO" "Input: $input_fastq"
    log "INFO" "Reference: $reference"
    log "INFO" "Output: $output_bam"
    log "INFO" "Preset: $preset"
    log "INFO" "Threads: $threads"
    
    # Run Minimap2 with time measurement
    log "INFO" "Running Minimap2 alignment..."
    /usr/bin/time -v minimap2 \
        -ax "$preset" \
        -t "$threads" \
        --secondary=no \
        "$reference" \
        "$input_fastq" \
        > "$temp_sam" \
        2> >(tee "$log_file" | grep -E "(real|user|sys|Maximum resident|Minor|Major)" > "$benchmark_file")
    
    # Convert SAM to BAM and sort
    log "INFO" "Converting to BAM and sorting..."
    /usr/bin/time -v samtools view -@ "$threads" -bS "$temp_sam" | \
        samtools sort -@ "$threads" -o "$output_bam" - \
        2>> "$benchmark_file"
    
    # Index BAM file
    log "INFO" "Indexing BAM file..."
    samtools index -@ "$threads" "$output_bam"
    
    # Clean up temporary files
    rm -f "$temp_sam"
    
    log "INFO" "Minimap2 alignment completed successfully"
}

# Function to calculate alignment statistics
calculate_alignment_stats() {
    local output_bam="$1"
    local output_dir="$2"
    local stats_file="${output_dir}/alignment_stats.txt"
    
    log "INFO" "Calculating alignment statistics..."
    
    {
        echo "=== ALIGNMENT STATISTICS ==="
        echo "Generated on: $(date)"
        echo ""
        
        echo "Basic Statistics:"
        samtools flagstat "$output_bam"
        echo ""
        
        echo "Detailed Statistics:"
        samtools stats "$output_bam" | head -n 50
        echo ""
        
        echo "Insert Size Statistics:"
        samtools stats "$output_bam" | grep "^IS" | head -n 20
        echo ""
        
        echo "Coverage Statistics:"
        samtools depth "$output_bam" | awk '{sum+=$3; count++} END {if(count>0) print "Average depth: " sum/count; print "Total positions: " count}'
        
    } > "$stats_file"
    
    log "INFO" "Alignment statistics saved to: $stats_file"
}

# Function to generate benchmark summary
generate_summary() {
    local output_dir="$1"
    local input_fastq="$2"
    local reference="$3"
    local output_bam="$4"
    local preset="$5"
    local threads="$6"
    
    local summary_file="${output_dir}/benchmark_summary.txt"
    local benchmark_file="${output_dir}/minimap2_benchmark.txt"
    
    log "INFO" "Generating benchmark summary..."
    
    {
        echo "====================================================="
        echo "           MINIMAP2 BENCHMARK SUMMARY"
        echo "====================================================="
        echo "Date: $(date)"
        echo "Host: $(hostname)"
        echo "User: $(whoami)"
        echo ""
        echo "CONFIGURATION:"
        echo "  Input FASTQ: $input_fastq"
        echo "  Reference: $reference"
        echo "  Output BAM: $output_bam"
        echo "  Preset: $preset"
        echo "  Threads: $threads"
        echo ""
        echo "INPUT FILE SIZES:"
        echo "  FASTQ: $(du -h "$input_fastq" | cut -f1)"
        echo "  Reference: $(du -h "$reference" | cut -f1)"
        echo "  Output BAM: $(du -h "$output_bam" | cut -f1)"
        echo ""
        echo "SYSTEM INFORMATION:"
        echo "  CPU: $(nproc) cores available"
        echo "  Memory: $(free -h | awk '/^Mem:/ {print $2}') total"
        echo "  OS: $(uname -s -r)"
        echo ""
        echo "PERFORMANCE METRICS:"
        if [[ -f "$benchmark_file" ]]; then
            cat "$benchmark_file"
        else
            echo "  Benchmark data not available"
        fi
        echo ""
        echo "MINIMAP2 VERSION:"
        minimap2 --version 2>/dev/null || echo "  Version information not available"
        echo ""
        echo "SAMTOOLS VERSION:"
        samtools --version | head -n 1
        echo ""
        echo "====================================================="
        
    } > "$summary_file"
    
    log "INFO" "Benchmark summary saved to: $summary_file"
}

# Main function
main() {
    local input_fastq=""
    local reference=""
    local output_bam=""
    local preset="$DEFAULT_PRESET"
    local threads="$DEFAULT_THREADS"
    local output_dir="$DEFAULT_OUTPUT_DIR"
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                input_fastq="$2"
                shift 2
                ;;
            -r|--reference)
                reference="$2"
                shift 2
                ;;
            -o|--output)
                output_bam="$2"
                shift 2
                ;;
            -p|--preset)
                preset="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -d|--output-dir)
                output_dir="$2"
                shift 2
                ;;
            -l|--log-level)
                LOG_LEVEL="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                log "ERROR" "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done
    
    # Validate required arguments
    if [[ -z "$input_fastq" || -z "$reference" || -z "$output_bam" ]]; then
        log "ERROR" "Missing required arguments"
        usage
        exit 1
    fi
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Check dependencies
    check_dependencies
    
    # Validate inputs
    validate_inputs "$input_fastq" "$reference"
    
    # Run benchmark
    run_minimap2_benchmark "$input_fastq" "$reference" "$output_bam" "$preset" "$threads" "$output_dir"
    
    # Calculate statistics
    calculate_alignment_stats "$output_bam" "$output_dir"
    
    # Generate summary
    generate_summary "$output_dir" "$input_fastq" "$reference" "$output_bam" "$preset" "$threads"
    
    log "INFO" "Minimap2 benchmark completed successfully!"
    log "INFO" "Results available in: $output_dir"
    log "INFO" "Summary file: ${output_dir}/benchmark_summary.txt"
}

# Execute main function if script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
