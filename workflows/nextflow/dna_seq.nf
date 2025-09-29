#!/usr/bin/env nextflow

/*
 * DNA-seq Analysis Pipeline
 * Pipeline completo para análise de sequenciamento de DNA
 * 
 * Autor: Gabriel Demetrios Lafis
 * Versão: 1.0.0
 */

nextflow.enable.dsl=2

// Parâmetros padrão
params.reads = "data/samples/*_{R1,R2}.fastq.gz"
params.genome = "data/reference/genome.fa"
params.outdir = "results/dna_seq"
params.threads = 8
params.min_qual = 20
params.min_len = 50

// Validação de parâmetros
if (!params.reads) {
    error "Por favor, forneça os arquivos de entrada com --reads"
}

if (!params.genome) {
    error "Por favor, forneça o genoma de referência com --genome"
}

// Log de início
log.info """
========================================
 Pipeline DNA-seq Analysis v1.0.0
========================================
 Reads       : ${params.reads}
 Genome      : ${params.genome}
 Outdir      : ${params.outdir}
 Threads     : ${params.threads}
 Min Quality : ${params.min_qual}
 Min Length  : ${params.min_len}
========================================
""".stripIndent()

// Processos

process FASTQC {
    tag "QC on $sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc -t ${params.threads} ${reads}
    """
}

process TRIMMING {
    tag "Trimming $sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz")
    path "${sample_id}_fastp.{json,html}"
    
    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --cut_front --cut_tail \
        --cut_mean_quality ${params.min_qual} \
        --length_required ${params.min_len} \
        --thread ${params.threads} \
        --json ${sample_id}_fastp.json \
        --html ${sample_id}_fastp.html
    """
}

process BWA_INDEX {
    tag "Indexing genome"
    
    input:
    path genome
    
    output:
    tuple path(genome), path("${genome}.*")
    
    script:
    """
    bwa index ${genome}
    """
}

process BWA_MEM {
    tag "Aligning $sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    tuple path(genome), path(index)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem -t ${params.threads} ${genome} ${reads[0]} ${reads[1]} | \
    samtools sort -@ ${params.threads} -o ${sample_id}.bam -
    samtools index ${sample_id}.bam
    """
}

process MARK_DUPLICATES {
    tag "MarkDup $sample_id"
    publishDir "${params.outdir}/dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam")
    path "${sample_id}.dedup_metrics.txt"
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.dedup.bam \
        -M ${sample_id}.dedup_metrics.txt
    samtools index ${sample_id}.dedup.bam
    """
}

process VARIANT_CALLING {
    tag "Calling variants $sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    tuple path(genome), path(index)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")
    
    script:
    """
    gatk HaplotypeCaller \
        -R ${genome} \
        -I ${bam} \
        -O ${sample_id}.vcf.gz
    """
}

process VARIANT_FILTERING {
    tag "Filtering variants $sample_id"
    publishDir "${params.outdir}/variants_filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz")
    
    script:
    """
    gatk VariantFiltration \
        -V ${vcf} \
        -O ${sample_id}.filtered.vcf.gz \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40"
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc . -n multiqc_report.html
    """
}

// Workflow principal

workflow {
    // Criar canal de reads
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
    
    // Controle de qualidade inicial
    FASTQC(reads_ch)
    
    // Trimagem
    TRIMMING(reads_ch)
    
    // Indexação do genoma
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    BWA_INDEX(genome_ch)
    
    // Alinhamento
    BWA_MEM(TRIMMING.out[0], BWA_INDEX.out)
    
    // Remoção de duplicatas
    MARK_DUPLICATES(BWA_MEM.out)
    
    // Chamada de variantes
    VARIANT_CALLING(MARK_DUPLICATES.out[0], BWA_INDEX.out)
    
    // Filtragem de variantes
    VARIANT_FILTERING(VARIANT_CALLING.out)
    
    // Relatório consolidado
    qc_files = FASTQC.out
        .mix(TRIMMING.out[1])
        .mix(MARK_DUPLICATES.out[1])
        .collect()
    
    MULTIQC(qc_files)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Duration: $workflow.duration"
}
