#!/bin/bash
# bench_star.sh - STAR Alignment Benchmarking Script
# 
# This script benchmarks STAR RNA-seq aligner performance with various parameters
# and generates detailed performance metrics for genomic data analysis pipeline.
#
# Usage: ./bench_star.sh <fastq_r1> <fastq_r2> <genome_dir> <gtf_file> <output_dir> <threads>
# Example: ./bench_star.sh sample_R1.fastq sample_R2.fastq /path/to/star_index annotations.gtf results 8
#
# Requirements:
# - STAR aligner installed and in PATH
# - samtools for BAM processing and metrics
# - Sufficient disk space for temporary files
# - Pre-built STAR genome index

set -euo pipefail

# Check command line arguments
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <fastq_r1> <fastq_r2> <genome_dir> <gtf_file> <output_dir> <threads>"
    echo "Example: $0 sample_R1.fastq sample_R2.fastq /path/to/star_index annotations.gtf results 8"
    exit 1
fi

# Input parameters
FASTQ_R1="$1"
FASTQ_R2="$2"
GENOME_DIR="$3"
GTF_FILE="$4"
OUTPUT_DIR="$5"
THREADS="$6"

# Validate input files
if [ ! -f "$FASTQ_R1" ]; then
    echo "Error: FASTQ R1 file not found: $FASTQ_R1"
    exit 1
fi

if [ ! -f "$FASTQ_R2" ]; then
    echo "Error: FASTQ R2 file not found: $FASTQ_R2"
    exit 1
fi

if [ ! -d "$GENOME_DIR" ]; then
    echo "Error: STAR genome directory not found: $GENOME_DIR"
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF annotation file not found: $GTF_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Setup logging
LOG_FILE="$OUTPUT_DIR/star_benchmark.log"
METRICS_FILE="$OUTPUT_DIR/star_metrics.txt"
TIME_FILE="$OUTPUT_DIR/star_timing.txt"

# Function for structured logging
log_message() {
    local level="$1"
    local message="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" | tee -a "$LOG_FILE"
}

# Function to collect system metrics
collect_system_metrics() {
    log_message "INFO" "Collecting system metrics..."
    {
        echo "=== SYSTEM INFORMATION ==="
        echo "Hostname: $(hostname)"
        echo "Date: $(date)"
        echo "CPU Info: $(grep 'model name' /proc/cpuinfo | head -n1 | cut -d':' -f2 | xargs)"
        echo "CPU Cores: $(nproc)"
        echo "Memory: $(free -h | grep '^Mem:' | awk '{print $2}')"
        echo "Disk Space: $(df -h . | tail -n1 | awk '{print $4}' | xargs) available"
        echo ""
    } >> "$METRICS_FILE"
}

# Function to collect STAR alignment metrics
collect_star_metrics() {
    local star_log="$1"
    log_message "INFO" "Collecting STAR alignment metrics..."
    {
        echo "=== STAR ALIGNMENT METRICS ==="
        if [ -f "$star_log" ]; then
            grep -E "Number of input reads|Uniquely mapped reads|% of reads mapped to multiple loci|% of reads unmapped" "$star_log" || true
        fi
        echo ""
    } >> "$METRICS_FILE"
}

# Function to collect BAM statistics
collect_bam_stats() {
    local bam_file="$1"
    log_message "INFO" "Collecting BAM statistics..."
    {
        echo "=== BAM STATISTICS ==="
        samtools flagstat "$bam_file" 2>/dev/null || echo "Error collecting BAM stats"
        echo ""
        echo "=== BAM INDEX STATS ==="
        samtools idxstats "$bam_file" 2>/dev/null | head -n5 || echo "Error collecting index stats"
        echo ""
    } >> "$METRICS_FILE"
}

# Main benchmarking function
run_star_benchmark() {
    log_message "INFO" "Starting STAR alignment benchmark"
    log_message "INFO" "Parameters: R1=$FASTQ_R1, R2=$FASTQ_R2, GenomeDir=$GENOME_DIR, GTF=$GTF_FILE, Threads=$THREADS"
    
    # Collect initial system metrics
    collect_system_metrics
    
    # Setup output files
    local output_prefix="$OUTPUT_DIR/star_aligned"
    local bam_file="${output_prefix}Aligned.sortedByCoord.out.bam"
    local log_file="${output_prefix}Log.final.out"
    
    # Run STAR alignment with timing
    log_message "INFO" "Running STAR alignment..."
    /usr/bin/time -v -o "$TIME_FILE" \
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --sjdbGTFfile "$GTF_FILE" \
        --readFilesIn "$FASTQ_R1" "$FASTQ_R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$output_prefix" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --quantMode GeneCounts \
        --outSAMattributes All \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        2>&1 | tee -a "$LOG_FILE"
    
    if [ $? -eq 0 ]; then
        log_message "INFO" "STAR alignment completed successfully"
        
        # Index BAM file for statistics
        log_message "INFO" "Indexing BAM file..."
        samtools index "$bam_file" 2>&1 | tee -a "$LOG_FILE"
        
        # Collect alignment metrics
        collect_star_metrics "$log_file"
        collect_bam_stats "$bam_file"
        
        # Append timing information
        {
            echo "=== TIMING INFORMATION ==="
            cat "$TIME_FILE"
            echo ""
        } >> "$METRICS_FILE"
        
        log_message "INFO" "Benchmark completed successfully"
        log_message "INFO" "Results saved to: $OUTPUT_DIR"
        log_message "INFO" "BAM file: $bam_file"
        log_message "INFO" "Metrics: $METRICS_FILE"
        log_message "INFO" "Timing: $TIME_FILE"
        
    else
        log_message "ERROR" "STAR alignment failed"
        exit 1
    fi
}

# Cleanup function
cleanup() {
    log_message "INFO" "Cleaning up temporary files..."
    # Remove STAR temporary directory if exists
    [ -d "${OUTPUT_DIR}/star_aligned_STARtmp" ] && rm -rf "${OUTPUT_DIR}/star_aligned_STARtmp"
}

# Set trap for cleanup
trap cleanup EXIT

# Execute benchmark
log_message "INFO" "Starting STAR benchmarking script"
run_star_benchmark
log_message "INFO" "STAR benchmarking completed"

# Example usage information
cat << 'EOF'

=== EXAMPLE USAGE ===

# Basic usage with paired-end FASTQ files
./bench_star.sh sample_R1.fastq.gz sample_R2.fastq.gz /path/to/star_index annotations.gtf results 8

# With different thread counts for performance comparison
for threads in 4 8 16; do
    ./bench_star.sh sample_R1.fastq.gz sample_R2.fastq.gz /path/to/star_index annotations.gtf "results_${threads}t" $threads
done

# Benchmark multiple samples
for sample in sample1 sample2 sample3; do
    ./bench_star.sh "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" /path/to/star_index annotations.gtf "results_${sample}" 8
done

EOF
