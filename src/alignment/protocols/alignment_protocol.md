# Protocolo de Alinhamento Genômico

## Descrição do Propósito

Este protocolo estabelece os padrões e procedimentos para alinhamento de sequências genômicas, garantindo consistência, reprodutibilidade e qualidade em todo o pipeline de análise de dados genômicos. O protocolo abrange desde a preparação dos dados até a validação final dos resultados de alinhamento.

## Requisitos

### Requisitos de Sistema
- **Sistema Operacional**: Linux/Unix (recomendado Ubuntu 18.04+)
- **Memória RAM**: Mínimo 16 GB (recomendado 32 GB+)
- **Espaço em disco**: 100 GB livres por amostra
- **CPU**: Mínimo 8 cores (recomendado 16+ cores)

### Dependências de Software
```bash
# Ferramentas de alinhamento
bwa >= 0.7.17
bowtie2 >= 2.4.0
star >= 2.7.0
minimap2 >= 2.17

# Processamento e análise
samtools >= 1.10
picard >= 2.23.0
bedtools >= 2.30.0

# Controle de qualidade
fastqc >= 0.11.9
qualimap >= 2.2.1

# Python e bibliotecas
python >= 3.8
pysam >= 0.16.0
biopython >= 1.78
```

### Dados de Entrada
- **Reads de sequenciamento**: Arquivos FASTQ (single-end ou paired-end)
- **Genoma de referência**: Arquivo FASTA indexado
- **Anotações** (opcional): Arquivos GTF/GFF para RNA-seq
- **Metadados**: Informações sobre amostras e condições experimentais

## Passo a Passo Padrão

### 1. Preparação e Validação dos Dados

#### 1.1 Controle de Qualidade Inicial
```bash
# Executar FastQC nos dados brutos
fastqc -t 8 -o qc_reports/ sample_R1.fastq sample_R2.fastq

# Análise de qualidade com MultiQC
multiqc qc_reports/ -o multiqc_report/
```

#### 1.2 Filtragem e Trimming (se necessário)
```bash
# Remoção de adaptadores e bases de baixa qualidade
trimmomatic PE -threads 8 \
  sample_R1.fastq sample_R2.fastq \
  sample_R1_trimmed.fastq sample_R1_unpaired.fastq \
  sample_R2_trimmed.fastq sample_R2_unpaired.fastq \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:15 MINLEN:36
```

### 2. Indexação do Genoma de Referência

#### 2.1 Indexação BWA
```bash
# Criar índices BWA para o genoma de referência
bwa index reference_genome.fasta
```

#### 2.2 Indexação SAMtools
```bash
# Indexar FASTA para SAMtools
samtools faidx reference_genome.fasta
```

#### 2.3 Dicionário de Sequência (Picard)
```bash
# Criar dicionário de sequências
picard CreateSequenceDictionary \
  R=reference_genome.fasta \
  O=reference_genome.dict
```

### 3. Alinhamento

#### 3.1 Alinhamento com BWA-MEM (DNA-seq)
```bash
# Alinhamento paired-end
bwa mem -t 16 -M \
  -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" \
  reference_genome.fasta \
  sample_R1_trimmed.fastq sample_R2_trimmed.fastq \
  > sample_aligned.sam
```

#### 3.2 Alinhamento com STAR (RNA-seq)
```bash
# Indexação do genoma para STAR
STAR --runMode genomeGenerate \
     --genomeDir star_index/ \
     --genomeFastaFiles reference_genome.fasta \
     --sjdbGTFfile annotations.gtf \
     --sjdbOverhang 100

# Alinhamento
STAR --genomeDir star_index/ \
     --readFilesIn sample_R1_trimmed.fastq sample_R2_trimmed.fastq \
     --outFileNamePrefix sample_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --quantMode GeneCounts
```

### 4. Processamento Pós-Alinhamento

#### 4.1 Conversão SAM para BAM e Ordenação
```bash
# Conversão e ordenação
samtools view -@ 8 -bS sample_aligned.sam | \
samtools sort -@ 8 -o sample_sorted.bam

# Indexação do BAM
samtools index sample_sorted.bam
```

#### 4.2 Marcação e Remoção de Duplicatas
```bash
# Marcar duplicatas com Picard
picard MarkDuplicates \
  I=sample_sorted.bam \
  O=sample_dedup.bam \
  M=sample_duplicates_metrics.txt \
  REMOVE_DUPLICATES=true

# Reindexar após remoção de duplicatas
samtools index sample_dedup.bam
```

### 5. Controle de Qualidade (QC)

#### 5.1 Estatísticas Básicas de Alinhamento
```bash
# Estatísticas gerais
samtools stats sample_dedup.bam > sample_stats.txt

# Estatísticas de flag
samtools flagstat sample_dedup.bam > sample_flagstat.txt
```

#### 5.2 Análise de Qualidade com Qualimap
```bash
# QC abrangente do BAM
qualimap bamqc \
  -bam sample_dedup.bam \
  -outdir qualimap_results/ \
  -outformat HTML \
  --java-mem-size=8G
```

#### 5.3 Análise de Cobertura
```bash
# Cobertura por posição
bedtools genomecov -ibam sample_dedup.bam -d > coverage_per_base.txt

# Histograma de cobertura
bedtools genomecov -ibam sample_dedup.bam > coverage_histogram.txt
```

### 6. Validação e Filtros Finais

#### 6.1 Filtros de Qualidade
```bash
# Filtrar por qualidade de mapeamento
samtools view -b -q 20 sample_dedup.bam > sample_filtered.bam

# Remover alinhamentos múltiplos
samtools view -b -F 256 sample_filtered.bam > sample_primary.bam
```

#### 6.2 Validação com Picard
```bash
# Validar integridade do arquivo BAM
picard ValidateSamFile \
  I=sample_primary.bam \
  MODE=SUMMARY
```

## Parâmetros Principais

### BWA-MEM
- **-t**: Número de threads (recomendado: 8-16)
- **-M**: Marcar reads mais curtos como secundários
- **-R**: String de grupo de reads (obrigatório)
- **-k**: Comprimento mínimo de seed (padrão: 19)
- **-w**: Largura da banda para alinhamento gapped (padrão: 100)

### STAR (RNA-seq)
- **--runThreadN**: Número de threads
- **--sjdbOverhang**: ReadLength-1 (padrão: 100)
- **--outFilterMismatchNmax**: Máximo de mismatches (padrão: 10)
- **--alignIntronMax**: Tamanho máximo de intron (padrão: 1000000)

### Filtros de Qualidade
- **MAPQ >= 20**: Qualidade mínima de mapeamento
- **Remoção de reads duplicados**: Recomendado para a maioria dos casos
- **Filtros de insert size**: Para dados paired-end (50-1000 bp típico)

### Parâmetros de Trimming
- **LEADING/TRAILING**: Remoção de bases de baixa qualidade (Q < 3)
- **SLIDINGWINDOW**: Janela deslizante 4:15 (Q médio >= 15)
- **MINLEN**: Comprimento mínimo após trimming (36 bp)

## Notas de Boas Práticas

### Nomenclatura e Organização
- **Convenção de nomes**: `{projeto}_{amostra}_{tipo}_{data}.{extensão}`
- **Estrutura de diretórios**: Separar raw data, processed data e results
- **Documentação**: Manter log detalhado de todos os parâmetros utilizados

### Controle de Versão
- Registrar versões de todas as ferramentas utilizadas
- Manter backups dos dados intermediários críticos
- Documentar modificações no protocolo com data e justificativa

### Validação de Resultados
- **Taxa de mapeamento**: > 85% para dados de boa qualidade
- **Duplicatas**: < 20% para dados WGS, < 50% para dados de enriquecimento
- **Distribuição de insert size**: Verificar se está dentro do esperado
- **Cobertura uniforme**: Avaliar viés de GC e regiões de baixa complexidade

### Performance e Recursos
- **Paralelização**: Usar threads de acordo com recursos disponíveis
- **Memória**: Monitorar uso de RAM, especialmente com STAR
- **I/O**: Usar discos SSD para dados temporários quando possível
- **Batch processing**: Processar múltiplas amostras em paralelo

### Troubleshooting Comum
- **Baixa taxa de mapeamento**: Verificar qualidade dos reads, contaminação
- **Alto uso de memória**: Reduzir threads, usar ferramentas alternativas
- **Erros de indexação**: Verificar integridade dos arquivos de referência
- **Reads não mapeados**: Analisar composição das sequências não alinhadas

## Links Internos

### Scripts e Ferramentas
- [BWA Wrapper](../tools/bwa_aligner.py) - Script automatizado para alinhamento BWA
- [STAR Pipeline](../tools/star_aligner.py) - Pipeline completo para RNA-seq
- [Quality Control](../tools/qualimap_qc.py) - Controle de qualidade automatizado
- [Batch Processor](../pipelines/batch_alignment.py) - Processamento em lote

### Utilitários
- [SAM/BAM Utils](../utils/sam_processing.py) - Utilitários para manipulação de arquivos
- [Coverage Analysis](../utils/coverage_analysis.py) - Análise de cobertura
- [QC Metrics](../utils/qc_metrics.py) - Coleta de métricas de qualidade

### Benchmarks e Comparações
- [Aligner Benchmark](../benchmarks/compare_aligners.py) - Comparação entre ferramentas
- [Performance Tests](../benchmarks/performance_test.py) - Testes de performance
- [Accuracy Assessment](../benchmarks/accuracy_benchmark.py) - Avaliação de precisão

### Documentação Adicional
- [README Principal](../README.md) - Visão geral do módulo de alinhamento
- [Configuração](../docs/setup_guide.md) - Guia de instalação e configuração
- [Troubleshooting](../docs/troubleshooting.md) - Soluções para problemas comuns
- [Best Practices](../docs/best_practices.md) - Práticas recomendadas avançadas

## Histórico de Revisões

| Versão | Data | Autor | Modificações |
|--------|------|-------|-------------|
| 1.0 | 2025-09-18 | Pipeline Team | Versão inicial do protocolo |
| | | | |

## Contato e Suporte

Para dúvidas, sugestões ou reportar problemas:
- Abrir issue no repositório do projeto
- Email: bioinformatics-support@projeto.org
- Documentação: https://projeto.github.io/genomic-pipeline/

---

**Nota**: Este protocolo deve ser revisado periodicamente e atualizado conforme novas versões de ferramentas e melhores práticas sejam estabelecidas. Consulte sempre a versão mais recente antes de iniciar novos projetos.
