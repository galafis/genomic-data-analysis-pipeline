#!/usr/bin/env nextflow

/*
 * Genomic Data Analysis Pipeline
 * RNA-seq Workflow
 * 
 * This workflow processes RNA-seq data from raw reads to differential expression analysis.
 * 
 * Author: Gabriel Demetrios Lafis
 * Date: June 2025
 */

// Default parameters
params.reads = "data/samples/*/fastq/*.fastq.gz"
params.genome = "data/reference/genome.fa"
params.annotation = "data/reference/genes.gtf"
params.outdir = "results/rna_seq"
params.help = false

// Show help message
if (params.help) {
    log.info"""
    =========================================
    Genomic Data Analysis Pipeline - RNA-seq
    =========================================
    
    Usage:
    
    nextflow run rna_seq.nf [options]
    
    Options:
      --reads         Path to input reads (default: $params.reads)
      --genome        Path to reference genome (default: $params.genome)
      --annotation    Path to genome annotation (default: $params.annotation)
      --outdir        Output directory (default: $params.outdir)
      --help          Show this help message and exit
    
    """.stripIndent()
    exit 0
}

// Log parameters
log.info"""
=========================================
Genomic Data Analysis Pipeline - RNA-seq
=========================================

Parameters:
  reads: ${params.reads}
  genome: ${params.genome}
  annotation: ${params.annotation}
  outdir: ${params.outdir}

=========================================
"""

// Create channels
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs_for_fastqc_ch }

// Reference genome and annotation files
genome_file = file(params.genome)
annotation_file = file(params.annotation)

// Process 1: Quality Control with FastQC
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from read_pairs_for_fastqc_ch
    
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

// Process 3: Alignment with STAR
process star_align {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from trimmed_reads_ch
    path genome from genome_file
    path annotation from annotation_file
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai") into aligned_bam_ch
    path "${sample_id}.Log.final.out" into star_logs_ch
    
    script:
    """
    # Create STAR index if it doesn't exist
    if [ ! -d star_index ]; then
        mkdir -p star_index
        STAR --runMode genomeGenerate \
            --runThreadN ${task.cpus} \
            --genomeDir star_index \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${annotation} \
            --sjdbOverhang 100
    fi
    
    # Align reads
    STAR --runMode alignReads \
        --runThreadN ${task.cpus} \
        --genomeDir star_index \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}. \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard \
        --quantMode GeneCounts
    
    # Index BAM file
    samtools index ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}

// Process 4: Quantification with featureCounts
process feature_counts {
    tag "$sample_id"
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai) from aligned_bam_ch
    path annotation from annotation_file
    
    output:
    path "${sample_id}.counts.txt" into feature_counts_ch
    
    script:
    """
    featureCounts -T ${task.cpus} \
        -p -s 2 \
        -a ${annotation} \
        -o ${sample_id}.counts.txt \
        ${bam}
    """
}

// Process 5: Merge Count Tables
process merge_counts {
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    path counts from feature_counts_ch.collect()
    
    output:
    path "merged_counts.txt" into merged_counts_ch
    
    script:
    """
    # Extract sample names
    for file in ${counts}; do
        sample=\$(basename \$file .counts.txt)
        echo \$sample >> sample_names.txt
    done
    
    # Merge count tables
    Rscript -e '
    library(data.table)
    library(dplyr)
    
    # Read sample names
    samples <- readLines("sample_names.txt")
    
    # Initialize merged counts
    merged <- NULL
    
    # Process each count file
    for (i in 1:length(samples)) {
        sample <- samples[i]
        file <- paste0(sample, ".counts.txt")
        
        # Skip header lines
        counts <- fread(file, skip=2)
        
        # Keep only gene ID and counts
        counts <- counts[, c(1, 7)]
        colnames(counts)[2] <- sample
        
        # Merge with existing data
        if (is.null(merged)) {
            merged <- counts
        } else {
            merged <- merge(merged, counts, by="Geneid")
        }
    }
    
    # Write merged counts
    write.table(merged, "merged_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)
    '
    """
}

// Process 6: Differential Expression Analysis with DESeq2
process deseq2_analysis {
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    
    input:
    path counts from merged_counts_ch
    
    output:
    path "deseq2_results.csv" into deseq2_results_ch
    path "normalized_counts.csv"
    path "pca_plot.pdf"
    path "heatmap.pdf"
    path "volcano_plot.pdf"
    
    script:
    """
    Rscript -e '
    library(DESeq2)
    library(pheatmap)
    library(ggplot2)
    library(EnhancedVolcano)
    
    # Read count data
    countData <- read.table("${counts}", header=TRUE, row.names=1)
    
    # Create sample metadata
    # This is a simplified example - in real analysis, you would have a proper metadata file
    samples <- colnames(countData)
    condition <- factor(ifelse(grepl("control", samples), "control", "treatment"))
    colData <- data.frame(row.names=samples, condition=condition)
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
    
    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds, contrast=c("condition", "treatment", "control"))
    res <- res[order(res$padj),]
    
    # Write results
    write.csv(as.data.frame(res), "deseq2_results.csv")
    
    # Normalized counts
    normalized_counts <- counts(dds, normalized=TRUE)
    write.csv(normalized_counts, "normalized_counts.csv")
    
    # PCA plot
    vsd <- vst(dds, blind=FALSE)
    pdf("pca_plot.pdf", width=8, height=6)
    plotPCA(vsd, intgroup="condition")
    dev.off()
    
    # Heatmap of top genes
    top_genes <- head(rownames(res[order(res$padj),]), 50)
    mat <- assay(vsd)[top_genes,]
    mat <- mat - rowMeans(mat)
    pdf("heatmap.pdf", width=8, height=10)
    pheatmap(mat, annotation_col=colData)
    dev.off()
    
    # Volcano plot
    pdf("volcano_plot.pdf", width=10, height=8)
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    pointSize = 3.0,
                    labSize = 3.0)
    dev.off()
    '
    """
}

// Process 7: Functional Enrichment Analysis
process functional_enrichment {
    publishDir "${params.outdir}/enrichment", mode: 'copy'
    
    input:
    path deseq2_results from deseq2_results_ch
    
    output:
    path "go_enrichment.csv"
    path "kegg_enrichment.csv"
    path "go_dotplot.pdf"
    path "kegg_dotplot.pdf"
    
    script:
    """
    Rscript -e '
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
    
    # Read DESeq2 results
    res <- read.csv("${deseq2_results}", row.names=1)
    
    # Get significant genes
    sig_genes <- rownames(res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,])
    
    # Convert gene IDs to Entrez IDs
    # This assumes your gene IDs are Ensembl IDs - adjust as needed
    gene_ids <- sig_genes
    entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_ids, keytype="ENSEMBL", column="ENTREZID")
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    # GO enrichment analysis
    go_results <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
    
    # KEGG pathway enrichment analysis
    kegg_results <- enrichKEGG(gene = entrez_ids,
                              organism = "hsa",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)
    
    # Save results
    write.csv(as.data.frame(go_results), "go_enrichment.csv")
    write.csv(as.data.frame(kegg_results), "kegg_enrichment.csv")
    
    # Create dotplots
    pdf("go_dotplot.pdf", width=10, height=8)
    print(dotplot(go_results, showCategory=15))
    dev.off()
    
    pdf("kegg_dotplot.pdf", width=10, height=8)
    print(dotplot(kegg_results, showCategory=15))
    dev.off()
    '
    """
}

// Process 8: Generate Report
process generate_report {
    publishDir "${params.outdir}/report", mode: 'copy'
    
    input:
    path star_logs from star_logs_ch.collect()
    
    output:
    path "rna_seq_report.html"
    
    script:
    """
    # Generate alignment statistics
    echo "<html><head><title>RNA-seq Analysis Report</title></head><body>" > rna_seq_report.html
    echo "<h1>RNA-seq Analysis Report</h1>" >> rna_seq_report.html
    
    echo "<h2>Alignment Statistics</h2>" >> rna_seq_report.html
    echo "<table border='1'>" >> rna_seq_report.html
    echo "<tr><th>Sample</th><th>Total Reads</th><th>Uniquely Mapped</th><th>% Uniquely Mapped</th></tr>" >> rna_seq_report.html
    
    for log in ${star_logs}; do
        sample=\$(basename \$log .Log.final.out)
        total=\$(grep "Number of input reads" \$log | cut -f2)
        unique=\$(grep "Uniquely mapped reads number" \$log | cut -f2)
        percent=\$(grep "Uniquely mapped reads %" \$log | cut -f2)
        
        echo "<tr><td>\$sample</td><td>\$total</td><td>\$unique</td><td>\$percent</td></tr>" >> rna_seq_report.html
    done
    
    echo "</table>" >> rna_seq_report.html
    
    echo "<h2>Differential Expression Analysis</h2>" >> rna_seq_report.html
    echo "<p>Differential expression analysis was performed using DESeq2. Results are available in the deseq2 directory.</p>" >> rna_seq_report.html
    
    echo "<h2>Functional Enrichment Analysis</h2>" >> rna_seq_report.html
    echo "<p>Functional enrichment analysis was performed using clusterProfiler. Results are available in the enrichment directory.</p>" >> rna_seq_report.html
    
    echo "</body></html>" >> rna_seq_report.html
    """
}

// Workflow completion message
workflow.onComplete {
    log.info"""
    =========================================
    Genomic Data Analysis Pipeline - RNA-seq
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

