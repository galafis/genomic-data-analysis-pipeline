#!/usr/bin/env nextflow

/*
 * Genomic Data Analysis Pipeline
 * DNA-seq Workflow
 * 
 * This workflow processes DNA-seq data from raw reads to variant calling and annotation.
 * 
 * Author: Gabriel Demetrios Lafis
 * Date: June 2025
 */

// Default parameters
params.reads = "data/samples/*/fastq/*.fastq.gz"
params.genome = "data/reference/genome.fa"
params.outdir = "results/dna_seq"
params.help = false

// Show help message
if (params.help) {
    log.info"""
    =========================================
    Genomic Data Analysis Pipeline - DNA-seq
    =========================================
    
    Usage:
    
    nextflow run dna_seq.nf [options]
    
    Options:
      --reads         Path to input reads (default: $params.reads)
      --genome        Path to reference genome (default: $params.genome)
      --outdir        Output directory (default: $params.outdir)
      --help          Show this help message and exit
    
    """.stripIndent()
    exit 0
}

// Log parameters
log.info"""
=========================================
Genomic Data Analysis Pipeline - DNA-seq
=========================================

Parameters:
  reads: ${params.reads}
  genome: ${params.genome}
  outdir: ${params.outdir}

=========================================
"""

// Create channels
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

// Reference genome file
genome_file = file(params.genome)

// Process 1: Quality Control with FastQC
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from read_pairs_ch
    
    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -q ${reads}
    """
}

// Process 2: Adapter Trimming with Trimmomatic
process trimmomatic {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from read_pairs_ch
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz") into trimmed_reads_ch
    
    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \
        ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Alignment with BWA-MEM
process bwa_mem {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from trimmed_reads_ch
    path genome from genome_file
    
    output:
    tuple val(sample_id), path("${sample_id}.bam") into aligned_bam_ch
    
    script:
    """
    # Index the genome if index doesn't exist
    if [ ! -f ${genome}.bwt ]; then
        bwa index ${genome}
    fi
    
    # Align reads
    bwa mem -t ${task.cpus} ${genome} ${reads[0]} ${reads[1]} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    
    # Index BAM file
    samtools index ${sample_id}.bam
    """
}

// Process 4: Mark Duplicates with Picard
process mark_duplicates {
    tag "$sample_id"
    publishDir "${params.outdir}/dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam) from aligned_bam_ch
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bai") into dedup_bam_ch
    path "${sample_id}.metrics.txt"
    
    script:
    """
    picard MarkDuplicates \
        INPUT=${bam} \
        OUTPUT=${sample_id}.dedup.bam \
        METRICS_FILE=${sample_id}.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        ASSUME_SORTED=true
    """
}

// Process 5: Base Recalibration with GATK
process base_recalibration {
    tag "$sample_id"
    publishDir "${params.outdir}/bqsr", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai) from dedup_bam_ch
    path genome from genome_file
    
    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bai") into recal_bam_ch
    
    script:
    """
    # Create sequence dictionary if it doesn't exist
    if [ ! -f ${genome.baseName}.dict ]; then
        picard CreateSequenceDictionary R=${genome} O=${genome.baseName}.dict
    fi
    
    # Index the genome if index doesn't exist
    if [ ! -f ${genome}.fai ]; then
        samtools faidx ${genome}
    fi
    
    # Base recalibration
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${genome} \
        --known-sites known_variants.vcf.gz \
        -O ${sample_id}.recal.table
    
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${genome} \
        --bqsr-recal-file ${sample_id}.recal.table \
        -O ${sample_id}.recal.bam
    """
}

// Process 6: Variant Calling with GATK HaplotypeCaller
process haplotype_caller {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai) from recal_bam_ch
    path genome from genome_file
    
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi") into gvcf_ch
    
    script:
    """
    gatk HaplotypeCaller \
        -I ${bam} \
        -R ${genome} \
        -ERC GVCF \
        -O ${sample_id}.g.vcf.gz
    """
}

// Process 7: Joint Genotyping
process joint_genotyping {
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    path genome from genome_file
    path gvcfs from gvcf_ch.map { it[1] }.collect()
    path indices from gvcf_ch.map { it[2] }.collect()
    
    output:
    path "cohort.vcf.gz" into vcf_ch
    path "cohort.vcf.gz.tbi"
    
    script:
    """
    # Create sample map file
    for gvcf in ${gvcfs}; do
        sample=\$(basename \$gvcf .g.vcf.gz)
        echo "\$sample \$gvcf" >> sample_map.txt
    done
    
    # Import GVCFs to GenomicsDB
    gatk GenomicsDBImport \
        --sample-name-map sample_map.txt \
        --genomicsdb-workspace-path genomicsdb \
        -L intervals.list
    
    # Joint genotyping
    gatk GenotypeGVCFs \
        -R ${genome} \
        -V gendb://genomicsdb \
        -O cohort.vcf.gz
    """
}

// Process 8: Variant Filtering
process variant_filtering {
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    path vcf from vcf_ch
    path genome from genome_file
    
    output:
    path "cohort.filtered.vcf.gz" into filtered_vcf_ch
    path "cohort.filtered.vcf.gz.tbi"
    
    script:
    """
    # SNP filtering
    gatk SelectVariants \
        -R ${genome} \
        -V ${vcf} \
        --select-type-to-include SNP \
        -O raw_snps.vcf.gz
    
    gatk VariantFiltration \
        -R ${genome} \
        -V raw_snps.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "SNP_filter" \
        -O filtered_snps.vcf.gz
    
    # Indel filtering
    gatk SelectVariants \
        -R ${genome} \
        -V ${vcf} \
        --select-type-to-include INDEL \
        -O raw_indels.vcf.gz
    
    gatk VariantFiltration \
        -R ${genome} \
        -V raw_indels.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "INDEL_filter" \
        -O filtered_indels.vcf.gz
    
    # Merge filtered variants
    gatk MergeVcfs \
        -I filtered_snps.vcf.gz \
        -I filtered_indels.vcf.gz \
        -O cohort.filtered.vcf.gz
    """
}

// Process 9: Variant Annotation with SnpEff
process variant_annotation {
    publishDir "${params.outdir}/annotated", mode: 'copy'
    
    input:
    path vcf from filtered_vcf_ch
    
    output:
    path "cohort.annotated.vcf.gz" into annotated_vcf_ch
    path "cohort.annotated.vcf.gz.tbi"
    path "snpeff_summary.html"
    
    script:
    """
    snpEff -Xmx8g -v GRCh38.99 \
        -stats snpeff_summary.html \
        ${vcf} > cohort.annotated.vcf
    
    bgzip cohort.annotated.vcf
    tabix -p vcf cohort.annotated.vcf.gz
    """
}

// Process 10: Generate Report
process generate_report {
    publishDir "${params.outdir}/report", mode: 'copy'
    
    input:
    path vcf from annotated_vcf_ch
    
    output:
    path "variant_report.html"
    
    script:
    """
    # Generate variant statistics
    bcftools stats ${vcf} > variant_stats.txt
    
    # Create HTML report
    echo "<html><head><title>Variant Calling Report</title></head><body>" > variant_report.html
    echo "<h1>Variant Calling Report</h1>" >> variant_report.html
    echo "<h2>Variant Statistics</h2>" >> variant_report.html
    echo "<pre>" >> variant_report.html
    cat variant_stats.txt >> variant_report.html
    echo "</pre>" >> variant_report.html
    echo "</body></html>" >> variant_report.html
    """
}

// Workflow completion message
workflow.onComplete {
    log.info"""
    =========================================
    Genomic Data Analysis Pipeline - DNA-seq
    =========================================
    
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    
    Results are available in: ${params.outdir}
    =========================================
    """
}

