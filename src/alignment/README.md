# 🧬 Alinhamento de Dados NGS (QC Pós-Alinhamento, Indexação, SAM/BAM)

## 📋 Visão Geral

Este módulo é responsável pelo alinhamento de sequências NGS contra genomas de referência, incluindo controle de qualidade pós-alinhamento, indexação eficiente, processamento de arquivos SAM/BAM, e integração com módulos de pré-processamento e chamada de variantes. Garante alinhamentos de alta qualidade e reprodutibilidade para análises downstream.

## 🎯 Funções do Módulo

• **Alinhamento de Reads**: BWA-MEM2, Bowtie2, STAR, minimap2 para diferentes tipos de dados NGS.
• **QC Pós-Alinhamento**: Avaliação de qualidade, taxa de mapeamento, distribuição MAPQ, cobertura.
• **Processamento SAM/BAM**: Conversão, ordenação, indexação, remoção de duplicatas (samtools/Picard).
• **Indexação de Referência**: Construção de índices otimizados para diferentes aligners.
• **Métricas de Qualidade**: Relatórios detalhados de performance e estatísticas de alinhamento.

## 📁 Estrutura de Diretórios

```
src/alignment/
├── aligners/
│   ├── bwa_mem2_wrapper.py
│   ├── bowtie2_wrapper.py
│   ├── star_wrapper.py
│   └── minimap2_wrapper.py
├── post_alignment_qc/
│   ├── alignment_metrics.py
│   ├── coverage_analysis.py
│   └── mapping_stats.py
├── sam_bam_processing/
│   ├── sam_bam_converter.py
│   ├── duplicate_remover.py
│   └── bam_indexer.py
└── reference_indexing/
    ├── bwa_indexer.py
    ├── bowtie2_indexer.py
    └── star_indexer.py
```

## 🔧 Ferramentas Recomendadas e Pipelines

• **BWA-MEM2**: Alinhamento rápido para WGS/WES com genomas grandes.
• **Bowtie2**: Alinhamento eficiente para reads curtos e dados de baixa complexidade.
• **STAR**: Alinhamento de RNA-seq com splice junction detection.
• **minimap2**: Long-read alignment para PacBio/Oxford Nanopore.
• **samtools**: Manipulação e processamento de arquivos SAM/BAM.
• **Picard**: Ferramentas Java para QC e processamento de dados NGS.

### Exemplos de Execução (scripts)

• **Indexação de referência (BWA-MEM2)**:
```bash
# Indexar genoma de referência
bwa-mem2 index -p reference_index reference.fasta
```

• **Alinhamento paired-end (BWA-MEM2)**:
```bash
bwa-mem2 mem -t 16 -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA" \
  reference_index \
  trimmed/sample_R1_paired.fq.gz trimmed/sample_R2_paired.fq.gz \
  | samtools sort -@ 8 -o alignment/sample.sorted.bam -
samtools index alignment/sample.sorted.bam
```

• **Alinhamento RNA-seq (STAR)**:
```bash
# Indexação
STAR --runMode genomeGenerate --genomeDir star_index/ \
  --genomeFastaFiles reference.fasta --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 99 --runThreadN 16

# Alinhamento
STAR --runMode alignReads --genomeDir star_index/ \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/sample_ --runThreadN 16
```

• **Long-read alignment (minimap2)**:
```bash
minimap2 -ax map-ont -t 16 reference.fasta long_reads.fastq.gz \
  | samtools sort -@ 8 -o alignment/sample_lr.sorted.bam -
samtools index alignment/sample_lr.sorted.bam
```

• **Remoção de duplicatas (Picard)**:
```bash
picard MarkDuplicates \
  INPUT=sample.sorted.bam \
  OUTPUT=sample.dedup.bam \
  METRICS_FILE=qc/sample.dedup_metrics.txt \
  REMOVE_DUPLICATES=true
samtools index sample.dedup.bam
```

• **QC pós-alinhamento (samtools + Picard)**:
```bash
# Estatísticas básicas
samtools flagstat sample.dedup.bam > qc/sample.flagstat
samtools stats sample.dedup.bam > qc/sample.stats

# Métricas de alinhamento
picard CollectAlignmentSummaryMetrics \
  INPUT=sample.dedup.bam \
  OUTPUT=qc/sample.alignment_metrics.txt \
  REFERENCE_SEQUENCE=reference.fasta

# Cobertura
samtools depth sample.dedup.bam > qc/sample.depth
```

## ⚙️ Orientações de Uso e Configuração

• **Instalação via conda (recomendado)**:
```bash
conda install -c bioconda bwa-mem2 bowtie2 star minimap2 samtools picard
```

• **Parâmetros sugeridos**:
  - **Threads**: 16-32 (ajuste conforme recursos disponíveis)
  - **Qualidade mínima MAPQ**: ≥20 para análises downstream
  - **Read Groups**: sempre adicionar informações de amostra/lane
  - **Remoção de duplicatas**: recomendado para WGS/WES

• **Layout de entradas/saídas**:
  - **Input**: `data/processed/trimmed/*_paired.fq.gz` (do pré-processamento)
  - **Output**: `data/aligned/*.sorted.bam`, `data/aligned/*.bai`, `qc/alignment/`

• **Validação pós-alinhamento**:
  - Taxa de mapeamento >85% para WGS/WES
  - Taxa de duplicação <20% (ideal <10%)
  - Distribuição uniforme de cobertura
  - MAPQ score médio >30

## 🔗 Integração com Preprocessing e Workflows

• **Pipeline preprocessing → alignment**:
  - Utilizar arquivos trimados paired como input
  - Preservar sample IDs e metadados de sequenciamento
  - Aplicar read groups consistentes

```bash
# Exemplo de pipeline integrado
nextflow run workflows/align.nf \
  --input "data/processed/trimmed/*_{R1,R2}_paired.fq.gz" \
  --reference references/genome.fasta \
  --outdir data/aligned/
```

• **Integração com variant calling**:
  - BAMs indexados como input para GATK/FreeBayes
  - Aplicar BQSR (Base Quality Score Recalibration) se necessário
  - Manter arquivos intermediários para troubleshooting

• **Workflows automatizados**:
```bash
# Snakemake
snakemake -j 32 --configfile config/alignment_config.yaml alignment_all

# Nextflow
nextflow run workflows/full_pipeline.nf \
  --input "data/raw/*_{R1,R2}.fastq.gz" --outdir results/
```

## 🧪 QC Pós-Alinhamento e Métricas-Chave

• **Métricas essenciais**:
  - Taxa de mapeamento (mapped reads %)
  - Distribuição MAPQ (>20, >30, >40)
  - Taxa de duplicação PCR
  - Insert size distribution (paired-end)
  - Cobertura média e uniformidade

• **Relatórios automatizados**:
```bash
# MultiQC para agregação
multiqc -o qc/alignment_multiqc/ qc/alignment/

# Relatórios customizados
python scripts/alignment_report.py \
  --bam_dir data/aligned/ --output qc/custom_report.html
```

• **Thresholds de qualidade**:
  - Taxa de mapeamento: >85% (WGS/WES), >70% (RNA-seq)
  - MAPQ médio: >25
  - Taxa de duplicação: <20%
  - Cobertura mínima: 10x (WES), 30x (WGS)

## ✅ Melhores Práticas

• **Read Groups obrigatórios**: sempre incluir @RG com ID, SM, PL, LB.
• **Versionamento**: fixar versões de aligners e parâmetros em lockfiles.
• **Backup de índices**: manter índices pré-construídos para referências padrão.
• **Paralelização**: otimizar threads por amostra vs. número de amostras.
• **Monitoramento**: logs detalhados de tempo/memória por etapa.
• **Testes de regressão**: validar com datasets conhecidos (tests/test_alignment.py).
• **Armazenamento**: compressão adequada de BAMs intermediários.

## 🔌 Configuração (exemplo de YAML)

```yaml
alignment:
  threads: 16
  min_mapq: 20
  remove_duplicates: true
  
  reference:
    genome_fasta: references/genome.fasta
    
  aligners:
    dna_default: bwa-mem2
    rna_default: star
    long_read: minimap2
    
  outputs:
    aligned_dir: data/aligned
    qc_dir: qc/alignment
    
  read_groups:
    platform: ILLUMINA
    add_rg: true
    
  post_processing:
    mark_duplicates: true
    create_index: true
    collect_metrics: true
```

## 🔗 Documentação Relacionada

• **Main Pipeline README**: ../README.md
• **Preprocessing README**: ../preprocessing/README.md
• **Variant Calling README**: ../variant_calling/README.md
• **Annotation README**: ../annotation/README.md
• **Workflows README**: ../workflows/README.md

## 🖼️ Padrão Visual do Projeto

• **Cabeçalhos com emojis temáticos**, listas claras, blocos de código com syntax highlighting, paleta consistente com demais READMEs.
• **Estrutura**: Visão Geral → Funções → Ferramentas/Exemplos → Configuração → Integração → QC → Melhores Práticas → Links.

## 📞 Suporte

• **Documentação**: Consulte `docs/alignment/` e issues do repositório.
• **Contato**: Equipe de desenvolvimento do projeto.
• **Troubleshooting**: Verifique logs em `logs/alignment/` e métricas QC.
