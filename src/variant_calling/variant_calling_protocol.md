# Protocolo de Chamada de Variantes

## Introdução

Este protocolo descreve o processo padrão para chamada de variantes genômicas a partir de dados de sequenciamento de alta throughput. O pipeline suporta múltiplos algoritmos de chamada de variantes e inclui etapas de pré-processamento, chamada, filtragem e anotação.

## Requisitos

### Hardware Recomendado
- **CPU**: Mínimo 8 cores, recomendado 16+ cores
- **RAM**: Mínimo 32GB, recomendado 64GB+
- **Armazenamento**: SSD com pelo menos 500GB livres
- **Rede**: Conexão estável para download de bancos de dados

### Software Necessário
- **bcftools** (>=1.15)
- **samtools** (>=1.15)
- **FreeBayes** (>=1.3.6)
- **GATK** (>=4.3.0)
- **VEP** ou **SnpEff** (anotação)
- **Python** (>=3.8) com pandas, numpy
- **R** (>=4.0) com VariantAnnotation, vcfR

### Bancos de Dados
- **dbSNP** (versão mais recente)
- **gnomAD** (frequências populacionais)
- **ClinVar** (significância clínica)
- **RefSeq/Ensembl** (anotação funcional)

## Preparação dos Dados

### 1. Verificação dos Arquivos de Entrada
```bash
# Verificar integridade dos BAMs
samtools quickcheck *.bam

# Verificar indexação
for bam in *.bam; do
  if [ ! -f "${bam}.bai" ]; then
    samtools index "$bam"
  fi
done
```

### 2. Validação da Qualidade
```bash
# Estatísticas básicas
samtools flagstat input.bam > flagstat_report.txt
samtools stats input.bam > stats_report.txt

# Verificar cobertura média
samtools depth input.bam | awk '{sum+=$3} END {print "Average depth:", sum/NR}'
```

## Chamada de Variantes

### Método 1: bcftools

#### Parâmetros Recomendados
- **-q**: Qualidade mínima de base (20)
- **-Q**: Qualidade mínima de mapeamento (20)
- **-d**: Profundidade máxima (250)
- **-a**: Anotações adicionais (AD,DP,SP)

#### Comandos de Exemplo
```bash
# Chamada básica
bcftools mpileup -f reference.fa -q 20 -Q 20 -d 250 \
  -a AD,DP,SP input.bam | \
bcftools call -mv -Oz -o variants_bcftools.vcf.gz

# Indexar resultado
bcftools index variants_bcftools.vcf.gz

# Normalizar e split de multi-alélicos
bcftools norm -f reference.fa -m -both -Oz \
  -o variants_bcftools_norm.vcf.gz variants_bcftools.vcf.gz
```

### Método 2: FreeBayes

#### Parâmetros Recomendados
- **--min-base-quality**: 20
- **--min-mapping-quality**: 20
- **--min-coverage**: 5
- **--min-alternate-fraction**: 0.2

#### Comandos de Exemplo
```bash
# Chamada com FreeBayes
freebayes -f reference.fa \
  --min-base-quality 20 \
  --min-mapping-quality 20 \
  --min-coverage 5 \
  --min-alternate-fraction 0.2 \
  input.bam > variants_freebayes.vcf

# Compressão e indexação
bgzip variants_freebayes.vcf
tabix -p vcf variants_freebayes.vcf.gz
```

### Método 3: GATK HaplotypeCaller

#### Parâmetros Recomendados
- **--minimum-mapping-quality**: 20
- **--base-quality-score-threshold**: 18
- **--standard-min-confidence-threshold-for-calling**: 30

#### Comandos de Exemplo
```bash
# Chamada com GATK
gatk HaplotypeCaller \
  -R reference.fa \
  -I input.bam \
  --minimum-mapping-quality 20 \
  --base-quality-score-threshold 18 \
  --standard-min-confidence-threshold-for-calling 30 \
  -O variants_gatk.vcf.gz

# Joint calling (múltiplas amostras)
gatk GenotypeGVCFs \
  -R reference.fa \
  -V sample1.g.vcf.gz \
  -V sample2.g.vcf.gz \
  -O cohort_variants.vcf.gz
```

## Filtragem de Variantes

### Filtragem Básica
```bash
# Filtros de qualidade com bcftools
bcftools filter -i 'QUAL>=30 && DP>=10 && DP<=250' \
  -Oz -o variants_filtered.vcf.gz variants_raw.vcf.gz

# Filtros específicos para SNPs
bcftools filter -i 'TYPE="snp" && QUAL>=30 && DP>=10 && MQ>=40' \
  -Oz -o snps_filtered.vcf.gz variants_raw.vcf.gz

# Filtros específicos para INDELs
bcftools filter -i 'TYPE="indel" && QUAL>=30 && DP>=10' \
  -Oz -o indels_filtered.vcf.gz variants_raw.vcf.gz
```

### Filtragem com GATK
```bash
# Variant Quality Score Recalibration (VQSR)
gatk VariantRecalibrator \
  -R reference.fa \
  -V input.vcf.gz \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 omni.vcf \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O snps.recal \
  --tranches-file snps.tranches
```

## Anotação de Variantes

### Anotação com VEP
```bash
# Anotação básica
vep -i variants_filtered.vcf.gz \
  --cache --dir_cache /path/to/vep/cache \
  --assembly GRCh38 \
  --fasta reference.fa \
  --vcf \
  --symbol \
  --terms SO \
  --canonical \
  --protein \
  --biotype \
  --af_1kg \
  --af_gnomad \
  --clin_sig \
  -o variants_annotated.vcf
```

### Anotação com SnpEff
```bash
# Download da base de dados
snpEff download GRCh38.99

# Anotação
snpEff ann GRCh38.99 variants_filtered.vcf.gz > variants_snpeff.vcf
```

## Controle de Qualidade

### Estatísticas das Variantes
```bash
# Contagem por tipo
bcftools stats variants_filtered.vcf.gz > variant_stats.txt

# Distribuição de qualidade
bcftools query -f '%QUAL\n' variants_filtered.vcf.gz | \
  sort -n > quality_distribution.txt

# Ti/Tv ratio
bcftools stats variants_filtered.vcf.gz | grep "TSTV"
```

### Visualização
```R
# Script R para plots de QC
library(vcfR)
library(ggplot2)

# Carregar VCF
vcf <- read.vcfR("variants_filtered.vcf.gz")

# Plot de distribuição de qualidade
qual_df <- data.frame(QUAL = getQUAL(vcf))
ggplot(qual_df, aes(x = QUAL)) + 
  geom_histogram(bins = 50) + 
  theme_minimal()
```

## Parâmetros de Configuração

### config/variant_calling.yaml
```yaml
variant_calling:
  tools:
    - bcftools
    - freebayes
    - gatk
  
  bcftools:
    min_base_quality: 20
    min_mapping_quality: 20
    max_depth: 250
    annotations: ["AD", "DP", "SP"]
  
  freebayes:
    min_base_quality: 20
    min_mapping_quality: 20
    min_coverage: 5
    min_alternate_fraction: 0.2
  
  gatk:
    min_mapping_quality: 20
    base_quality_threshold: 18
    min_confidence_threshold: 30
  
  filtering:
    min_quality: 30
    min_depth: 10
    max_depth: 250
    min_mq: 40
```

## Boas Práticas

### Planejamento
1. **Definir objetivos**: SNPs, INDELs, CNVs?
2. **Escolher ferramenta**: Considerar tipo de dados e recursos
3. **Validação**: Sempre validar subset com Sanger
4. **Backup**: Manter cópias dos dados originais

### Execução
1. **Paralelização**: Usar scatter-gather para grandes datasets
2. **Monitoramento**: Acompanhar uso de recursos
3. **Checkpoints**: Salvar resultados intermediários
4. **Logs**: Documentar todos os parâmetros utilizados

### Validação
1. **Ti/Tv ratio**: ~2.0-2.1 para exomas, ~2.8-3.0 para genomas
2. **Concordância**: Comparar múltiplos callers
3. **Mendelian violations**: Para trios familiares
4. **Hardy-Weinberg**: Para estudos populacionais

## Troubleshooting

### Problemas Comuns
- **Baixa concordância**: Verificar qualidade dos alinhamentos
- **Excesso de falsos positivos**: Ajustar filtros de qualidade
- **Missing calls**: Reduzir stringência dos filtros
- **Problemas de memória**: Processar por cromossomo

### Otimização de Performance
```bash
# Processamento paralelo por região
parallel -j 8 "bcftools mpileup -r chr{} -f ref.fa input.bam | \
  bcftools call -mv > chr{}.vcf" ::: {1..22} X Y

# Concatenar resultados
bcftools concat chr*.vcf -Oz -o final.vcf.gz
```

## Links e Referências

### Documentação Interna
- [Protocolo de Alinhamento](../alignment/alignment_protocol.md)
- [Protocolo de QC](../quality_control/qc_protocol.md)
- [Configuração do Pipeline](../../config/README.md)

### Ferramentas
- [bcftools](http://samtools.github.io/bcftools/)
- [FreeBayes](https://github.com/freebayes/freebayes)
- [GATK](https://gatk.broadinstitute.org/)
- [VEP](https://ensembl.org/info/docs/tools/vep/)

### Bancos de Dados
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)

---

**Versão**: 1.0
**Última atualização**: September 2025
**Autor**: Genomic Data Analysis Pipeline Team
