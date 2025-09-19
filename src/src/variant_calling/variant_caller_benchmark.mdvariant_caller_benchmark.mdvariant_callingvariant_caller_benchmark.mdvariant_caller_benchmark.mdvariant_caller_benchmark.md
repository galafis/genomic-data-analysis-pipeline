# Variant Caller Benchmarking

## Introdução

O benchmarking de variant callers é essencial para avaliar a performance e precisão de diferentes ferramentas de detecção de variantes genômicas. Este documento apresenta uma metodologia padronizada para comparar ferramentas populares como GATK HaplotypeCaller, FreeBayes e bcftools mpileup/call.

## Ferramentas Avaliadas

### GATK HaplotypeCaller
- **Algoritmo**: Assembly-based local realignment
- **Pontos fortes**: Alta precisão em regiões complexas, suporte robusto para indels
- **Limitações**: Maior consumo de recursos computacionais

### FreeBayes
- **Algoritmo**: Bayesian genetic variant detector
- **Pontos fortes**: Rápido, boa performance em populações
- **Limitações**: Pode gerar mais falsos positivos em regiões repetitivas

### bcftools mpileup/call
- **Algoritmo**: Statistical approach based on read depth
- **Pontos fortes**: Muito rápido, baixo uso de memória
- **Limitações**: Menor sensibilidade para indels complexos

## Critérios de Avaliação

### Métricas de Qualidade

- **Sensibilidade (Recall)**: TP / (TP + FN)
- **Precisão (Precision)**: TP / (TP + FP)
- **F1-Score**: 2 × (Precisão × Sensibilidade) / (Precisão + Sensibilidade)
- **Especificidade**: TN / (TN + FP)

### Métricas de Performance

- **Tempo de execução** (wall-clock time)
- **Uso de memória RAM** (pico)
- **Uso de CPU** (média)
- **I/O disk** (operações de leitura/escrita)

## Metodologia Recomendada

### 1. Preparação do Dataset

```bash
# Utilizar datasets de referência (ex: Genome in a Bottle)
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/

# Preparar BAM files alinhados com BWA-MEM
bwa mem -t 8 -R '@RG\tID:sample\tSM:NA12878' \
  reference.fa sample_R1.fastq sample_R2.fastq | \
  samtools sort -@ 8 -o sample.bam
```

### 2. Execução dos Variant Callers

#### GATK HaplotypeCaller
```bash
#!/bin/bash
time gatk --java-options "-Xmx8g" HaplotypeCaller \
  -R reference.fa \
  -I sample.bam \
  -O variants_gatk.vcf.gz \
  --native-pair-hmm-threads 8 \
  --standard-min-confidence-threshold-for-calling 10.0
```

#### FreeBayes
```bash
#!/bin/bash
time freebayes -f reference.fa \
  -b sample.bam \
  --min-base-quality 20 \
  --min-mapping-quality 20 \
  --ploidy 2 > variants_freebayes.vcf
```

#### bcftools
```bash
#!/bin/bash
time bcftools mpileup -Ou -f reference.fa sample.bam | \
  bcftools call -mv -Oz -o variants_bcftools.vcf.gz
```

### 3. Normalização e Filtragem

```bash
# Normalizar variantes
for vcf in *.vcf.gz; do
  bcftools norm -f reference.fa -m-both $vcf | \
  bcftools norm -d both -o ${vcf%.vcf.gz}_norm.vcf.gz -O z
done

# Aplicar filtros de qualidade
bcftools filter -i 'QUAL>=30 && DP>=10' \
  variants_norm.vcf.gz -o variants_filtered.vcf.gz
```

### 4. Comparação com Truth Set

```bash
# Usando hap.py (GIAB)
python hap.py truth.vcf.gz query.vcf.gz \
  -f confident_regions.bed \
  -r reference.fa \
  -o benchmark_results

# Usando bcftools isec
bcftools isec -p comparison_dir truth.vcf.gz query.vcf.gz
```

## Exemplo de Resultados

### Tabela Comparativa - Dataset NA12878 (30x coverage)

| Ferramenta | SNPs (TP) | SNPs (FP) | SNPs (FN) | Sensibilidade | Precisão | F1-Score | Tempo (min) | RAM (GB) |
|------------|-----------|-----------|-----------|---------------|----------|----------|-------------|----------|
| GATK-HC    | 2,847,392 | 1,234     | 2,891     | 99.90%        | 99.96%   | 99.93%   | 127         | 8.2      |
| FreeBayes  | 2,845,821 | 2,187     | 4,462     | 99.84%        | 99.92%   | 99.88%   | 89          | 4.1      |
| bcftools   | 2,841,459 | 856       | 8,824     | 99.69%        | 99.97%   | 99.83%   | 34          | 1.8      |

### Métricas de Indels

| Ferramenta | Indels (TP) | Indels (FP) | Indels (FN) | Sensibilidade | Precisão | F1-Score |
|------------|-------------|-------------|-------------|---------------|----------|-----------|
| GATK-HC    | 347,823     | 892         | 1,245       | 99.64%        | 99.74%   | 99.69%   |
| FreeBayes  | 346,891     | 1,634       | 2,177       | 99.38%        | 99.53%   | 99.45%   |
| bcftools   | 342,156     | 567         | 6,912       | 98.02%        | 99.83%   | 98.92%   |

## Monitoramento de Performance

### Script de Benchmark Automático

```bash
#!/bin/bash
# benchmark_variant_callers.sh

REF="reference.fa"
BAM="sample.bam"
OUTDIR="benchmark_results"

mkdir -p $OUTDIR

# Função para medir recursos
measure_resources() {
    local cmd="$1"
    local output="$2"
    
    /usr/bin/time -v -o "${output}_resources.txt" \
    bash -c "$cmd"
}

# GATK
echo "Running GATK HaplotypeCaller..."
measure_resources \
    "gatk HaplotypeCaller -R $REF -I $BAM -O $OUTDIR/gatk.vcf.gz" \
    "$OUTDIR/gatk"

# FreeBayes
echo "Running FreeBayes..."
measure_resources \
    "freebayes -f $REF $BAM > $OUTDIR/freebayes.vcf" \
    "$OUTDIR/freebayes"

# bcftools
echo "Running bcftools..."
measure_resources \
    "bcftools mpileup -Ou -f $REF $BAM | bcftools call -mv -Oz -o $OUTDIR/bcftools.vcf.gz" \
    "$OUTDIR/bcftools"
```

## Boas Práticas

### Configuração do Ambiente

1. **Recursos computacionais**: Mínimo 8 cores, 16GB RAM
2. **Storage**: SSD recomendado para I/O intensivo
3. **Reproducibilidade**: Fixar versões das ferramentas e parâmetros

### Validação Estatística

```bash
# Calcular intervalo de confiança para métricas
R --vanilla << 'EOF'
results <- read.csv("benchmark_results.csv")
binom.test(results$TP, results$TP + results$FN, conf.level=0.95)
EOF
```

### Controle de Qualidade

- Verificar distribuição de qualidades (bcftools stats)
- Analisar viés de strand (bcftools +strand-bias)
- Avaliar cobertura uniforme (mosdepth)

## Interpretação dos Resultados

### Critérios de Seleção

1. **Para aplicações clínicas**: Priorizar sensibilidade (minimizar FN)
2. **Para descoberta populacional**: Balance entre sensibilidade e precisão
3. **Para análises em larga escala**: Considerar tempo/custo computacional

### Limitações do Benchmarking

- Truth sets podem conter viés
- Regiões difíceis (repetitivas, GC extremo) sub-representadas
- Variabilidade entre diferentes populações

## Ferramentas Auxiliares

### Análise de Concordância

```bash
# Venn diagram de variantes
bedtools intersect -a gatk.vcf -b freebayes.vcf -u | wc -l

# Análise detalhada com RTG Tools
rtg vcfeval -b truth.vcf.gz -c query.vcf.gz \
  -t reference.sdf -o eval_output
```

### Visualização

```bash
# Plot de métricas com R
R --vanilla << 'EOF'
library(ggplot2)
data <- read.csv("metrics.csv")
ggplot(data, aes(x=Tool, y=F1Score)) + 
  geom_bar(stat="identity") + 
  theme_minimal()
EOF
```

## Referências

1. Zook, J.M. et al. "Extensive sequencing of seven human genomes to characterize benchmark reference materials." *Sci Data* 3, 160025 (2016)

2. Krusche, P. et al. "Best practices for benchmarking germline small-variant calls in human genomes." *Nat Biotechnol* 37, 555–560 (2019)

3. McKenna, A. et al. "The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data." *Genome Research* 20:1297-1303 (2010)

4. Garrison, E. & Marth, G. "Haplotype-based variant detection from short-read sequencing." *arXiv preprint* arXiv:1207.3907 (2012)

5. Li, H. et al. "The Sequence Alignment/Map format and SAMtools." *Bioinformatics* 25, 2078-2079 (2009)

6. Cleary, J.G. et al. "Comparing Variant Call Files for Performance Benchmarking of Next-Generation Sequencing Variant Calling Pipelines." *bioRxiv* (2015)

## Anexos

### Configurações Recomendadas

```yaml
# config.yaml
gatk:
  java_opts: "-Xmx8g"
  min_confidence: 10.0
  native_threads: 8

freebayes:
  min_base_quality: 20
  min_mapping_quality: 20
  ploidy: 2

bcftools:
  mpileup_opts: "-q 20 -Q 20"
  call_opts: "-mv"
```

### Scripts de Automação

Ver diretório `scripts/` para implementações completas dos pipelines de benchmarking.
