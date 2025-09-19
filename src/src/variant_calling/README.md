# Variant Calling Module

## 📋 Visão Geral

Este é o módulo de chamada de variantes do pipeline genômico, responsável por identificar e anotar variantes genômicas (SNPs, indels, CNVs, SVs) a partir de dados de sequenciamento DNA-seq e RNA-seq. O módulo implementa múltiplas ferramentas de variant calling com benchmarking comparativo e integração downstream completa.

## 🗂️ Estrutura do Módulo

### Arquivos Atuais
- `variant_caller_benchmark.md` - Documentação de benchmarking de variant callers

### ⚠️ Componentes Ausentes Identificados

#### 1. Scripts Executáveis (`scripts/`)
- `run_gatk_haplotypecaller.sh` - Pipeline GATK
- `run_freebayes.sh` - Pipeline FreeBayes  
- `run_bcftools.sh` - Pipeline bcftools
- `run_somatic_calling.sh` - Chamada de variantes somáticas
- `batch_variant_calling.sh` - Processamento em lote
- `variant_calling_qc.sh` - Controle de qualidade

#### 2. Configurações (`config/`)
- `gatk_config.yaml` - Configurações GATK
- `freebayes_config.yaml` - Configurações FreeBayes
- `bcftools_config.yaml` - Configurações bcftools
- `filtering_rules.yaml` - Regras de filtragem
- `annotation_config.yaml` - Configurações de anotação

#### 3. Workflows (`workflows/`)
- `nextflow/variant_calling.nf` - Pipeline principal Nextflow
- `snakemake/Snakefile` - Pipeline Snakemake
- `cwl/variant_calling.cwl` - Workflow CWL

#### 4. Containers (`containers/`)
- `Dockerfile.gatk` - Container GATK
- `Dockerfile.freebayes` - Container FreeBayes
- `Dockerfile.annotation` - Container anotação

#### 5. Utilitários Python (`utils/`)
- `variant_caller.py` - Classes de variant callers
- `vcf_processor.py` - Processamento de VCFs
- `annotation_parser.py` - Parser de anotações
- `metrics_calculator.py` - Cálculo de métricas

#### 6. Análises e Relatórios (`analysis/`)
- `generate_reports.R` - Geração de relatórios
- `variant_stats.py` - Estatísticas de variantes
- `comparison_plots.R` - Gráficos comparativos

#### 7. Dados de Teste (`test_data/`)
- `sample.bam` - Dados de teste
- `reference.fa` - Genoma de referência
- `truth_set.vcf.gz` - Conjunto verdade

#### 8. Troubleshooting (`troubleshooting/`)
- `common_issues.md` - Problemas comuns
- `debugging_guide.md` - Guia de depuração
- `performance_tuning.md` - Otimização

#### 9. Integração Downstream (`integration/`)
- `to_annotation.sh` - Integração com anotação
- `to_visualization.sh` - Integração com visualização
- `export_formats.py` - Conversão de formatos

## 🚀 Instalação e Configuração

### Pré-requisitos
```bash
# Ferramentas obrigatórias
- GATK >= 4.2.0
- FreeBayes >= 1.3.0
- bcftools >= 1.15
- VEP >= 105
- SnpEff >= 5.0

# Dependências Python
pip install -r requirements.txt

# Dependências R
Rscript install_packages.R
```

### Configuração do Ambiente
```bash
# Ativar ambiente Conda
conda activate genomic-pipeline

# Configurar paths
export GATK_PATH=/path/to/gatk
export VEP_CACHE=/path/to/vep_cache
export REFERENCE_GENOME=/path/to/reference.fa
```

## 📊 Uso

### 1. Chamada de Variantes Individual

#### GATK HaplotypeCaller
```bash
./scripts/run_gatk_haplotypecaller.sh \
  --input sample.bam \
  --reference reference.fa \
  --output variants_gatk.vcf.gz \
  --config config/gatk_config.yaml
```

#### FreeBayes
```bash
./scripts/run_freebayes.sh \
  --input sample.bam \
  --reference reference.fa \
  --output variants_freebayes.vcf \
  --config config/freebayes_config.yaml
```

#### bcftools
```bash
./scripts/run_bcftools.sh \
  --input sample.bam \
  --reference reference.fa \
  --output variants_bcftools.vcf.gz \
  --config config/bcftools_config.yaml
```

### 2. Processamento em Lote
```bash
# Processar múltiplas amostras
./scripts/batch_variant_calling.sh \
  --input-dir /path/to/bam_files/ \
  --output-dir /path/to/results/ \
  --caller gatk,freebayes,bcftools \
  --threads 16
```

### 3. Workflows Gerenciados

#### Nextflow
```bash
nextflow run workflows/nextflow/variant_calling.nf \
  --input_dir "data/bam/*.bam" \
  --reference_genome "data/reference/genome.fa" \
  --output_dir "results/variant_calling" \
  --callers "gatk,freebayes,bcftools"
```

#### Snakemake
```bash
snakemake -s workflows/snakemake/Snakefile \
  --configfile config/variant_calling_config.yaml \
  --cores 16
```

### 4. Análise Comparativa
```bash
# Executar benchmarking
python utils/run_benchmark.py \
  --callers gatk,freebayes,bcftools \
  --truth-set data/truth_sets/NA12878.vcf.gz \
  --confident-regions data/truth_sets/confident_regions.bed \
  --output results/benchmark/
```

## 📈 Análise de Resultados

### Métricas de Qualidade
```bash
# Gerar estatísticas
python utils/variant_stats.py \
  --input variants.vcf.gz \
  --output stats/variant_statistics.txt

# Gerar relatório HTML
Rscript analysis/generate_reports.R \
  --input results/benchmark/ \
  --output reports/variant_calling_report.html
```

### Visualizações
- Distribuições de qualidade por caller
- Curvas ROC comparativas
- Heatmaps de concordância
- Gráficos de performance computacional

## 🔧 Troubleshooting

### Problemas Comuns

#### 1. Erro de Memória (OutOfMemoryError)
```bash
# Aumentar heap Java
export JAVA_OPTS="-Xmx32g"

# Usar processamento por regiões
gatk HaplotypeCaller \
  --intervals chr1:1-1000000 \
  --native-pair-hmm-threads 4
```

#### 2. Baixa Concordância entre Callers
- Verificar qualidade dos alignments
- Ajustar parâmetros de filtragem
- Normalizar variantes com bcftools norm

#### 3. Performance Lenta
- Usar paralelização por cromossomo
- Implementar scatter-gather
- Utilizar containers otimizados

### Logs e Debugging
```bash
# Habilitar logs detalhados
export LOG_LEVEL=DEBUG

# Verificar métricas de sistema
./troubleshooting/monitor_resources.sh
```

## 🧪 Testes

### Suite de Testes
```bash
# Executar testes unitários
python -m pytest tests/

# Teste de integração
bash tests/integration_test.sh

# Teste de benchmark
python tests/test_benchmark.py
```

### Validação de Resultados
```bash
# Comparar com resultados esperados
python tests/validate_results.py \
  --test-results results/test/ \
  --expected-results tests/expected/
```

## 🔄 Integração Downstream

### 1. Anotação de Variantes
```bash
./integration/to_annotation.sh \
  --input variants.vcf.gz \
  --output annotated_variants.vcf.gz
```

### 2. Visualização
```bash
./integration/to_visualization.sh \
  --input variants.vcf.gz \
  --output visualization_data/
```

### 3. Exportação de Formatos
```bash
# Converter para formatos específicos
python integration/export_formats.py \
  --input variants.vcf.gz \
  --output-format tsv,json,xml \
  --output-dir exports/
```

## 📋 Benchmarking

### Conjuntos de Dados de Referência
- **Genome in a Bottle (GIAB)**: NA12878, NA24385, NA24631
- **Platinum Genomes**: Illumina high-confidence variants
- **1000 Genomes**: Variabilidade populacional

### Métricas Avaliadas
- **Sensibilidade (Recall)**: TP / (TP + FN)
- **Precisão (Precision)**: TP / (TP + FP)  
- **F1-Score**: 2 × (Precisão × Sensibilidade) / (Precisão + Sensibilidade)
- **Performance Computacional**: Tempo, memória, CPU

### Relatórios de Benchmark
Ver `variant_caller_benchmark.md` para análises detalhadas comparativas.

## 🔒 Boas Práticas

### 1. Controle de Qualidade
- Filtrar variantes de baixa qualidade (QUAL < 30)
- Remover regiões de baixa complexidade
- Validar anotações funcionais

### 2. Reprodutibilidade
- Versionamento de ferramentas e bancos de dados
- Containerização com Docker/Singularity
- Rastreabilidade completa de parâmetros

### 3. Escalabilidade
- Paralelização por cromossomo/região
- Uso de sistemas de workflow (Nextflow/Snakemake)
- Otimização para HPC e nuvem

## 🔗 Referências

1. McKenna, A. et al. "The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data." Genome Research (2010)
2. Garrison, E. & Marth, G. "Haplotype-based variant detection from short-read sequencing." arXiv:1207.3907 (2012)
3. Li, H. et al. "The Sequence Alignment/Map format and SAMtools." Bioinformatics (2009)
4. Zook, J.M. et al. "Extensive sequencing of seven human genomes to characterize benchmark reference materials." Scientific Data (2016)

## 📞 Suporte

Para questões técnicas ou sugestões:
- **Issues**: [GitHub Issues](https://github.com/galafis/genomic-data-analysis-pipeline/issues)
- **Documentação**: Ver `docs/` no diretório raiz
- **Wiki**: [Project Wiki](https://github.com/galafis/genomic-data-analysis-pipeline/wiki)

---

**Última atualização**: $(date '+%Y-%m-%d')  
**Versão do módulo**: 1.0.0  
**Compatibilidade**: Pipeline genomic-data-analysis-pipeline v2.0+
