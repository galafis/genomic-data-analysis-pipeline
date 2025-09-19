# Análises e Relatórios do Variant Calling

Este diretório contém scripts e utilitários (Python e R) para análise de resultados, geração de estatísticas descritivas, visualização de dados e relatórios de controle de qualidade do módulo de variant calling. O objetivo é padronizar análises downstream, benchmarking comparativo e geração de relatórios profissionais, mantendo o mesmo padrão visual/organizacional do pipeline.

## 📋 Visão Geral

Componentes típicos (exemplos sugeridos):
- generate_reports.R — geração de relatórios HTML/PDF com estatísticas
- variant_stats.py — estatísticas descritivas de variantes (contagens, distribuições, métricas)
- comparison_plots.R — gráficos comparativos entre callers (ROC, precisão-recall, concordância)
- qc_dashboard.R — dashboard interativo Shiny para exploração de resultados
- benchmark_analysis.py — análise detalhada de benchmarking contra truth sets
- variant_visualization.py — plots específicos de variantes (Manhattan, circos, heatmaps)
- integration_plots.R — visualizações de integração com módulos downstream
- export_summaries.py — exportação de sumários para formatos específicos

Estrutura recomendada:
```
analysis/
├── generate_reports.R
├── variant_stats.py
├── comparison_plots.R
├── qc_dashboard.R
├── benchmark_analysis.py
├── variant_visualization.py
├── integration_plots.R
└── export_summaries.py
```

## 🔧 Principais Scripts e Funcionalidades

### 1) generate_reports.R
- render_html_report(input_dir, output_file, template="variant_report.Rmd", include_plots=TRUE)
- render_pdf_report(input_dir, output_file, executive_summary=TRUE)
- batch_report_generation(sample_dirs, output_dir, parallel=TRUE)
  - parallel: processamento paralelo com furrr/parallel
  - executive_summary: resumo executivo para stakeholders

### 2) variant_stats.py
- calculate_basic_stats(vcf_file, output_json=None, output_txt=None)
- variant_type_distribution(vcf_file, plot_output="variant_types.png")
- quality_metrics_summary(vcf_file, caller_comparison=False)
- depth_af_distributions(vcf_file, bins=50, output_plots="depth_af_plots.png")

### 3) comparison_plots.R
- plot_caller_venn(vcf_list, sample_name, output="venn_diagram.png")
- roc_curves_comparison(truth_vcf, caller_vcfs, confident_bed=NULL)
- precision_recall_plots(benchmark_results, facet_by="caller")
- concordance_heatmap(caller_matrix, cluster_rows=TRUE, cluster_cols=TRUE)

### 4) qc_dashboard.R (Shiny)
- launch_qc_dashboard(data_dir, port=3838, host="localhost")
- Componentes: métricas por amostra, distribuições interativas, filtros dinâmicos
- Suporte para múltiplos callers e comparações lado a lado

### 5) benchmark_analysis.py
- comprehensive_benchmark(truth_vcf, caller_vcfs, confident_regions, output_dir)
- stratified_analysis(benchmark_results, strata=["snp", "indel", "complex"])
- performance_metrics(tp, fp, fn, runtime_stats=None)
- generate_benchmark_report(results_dict, output_html="benchmark_report.html")

### 6) variant_visualization.py
- manhattan_plot(vcf_file, output="manhattan.png", highlight_genes=None)
- circos_plot(vcf_file, reference_fasta, output="circos.svg")
- variant_heatmap(vcf_matrix, sample_metadata=None, clustering="ward")
- af_vs_depth_scatter(vcf_file, color_by="caller", alpha=0.7)

## 🧪 Exemplos Práticos de Uso

### 1) Geração de relatório completo
```bash
Rscript analysis/generate_reports.R \
  --input results/variant_calling/ \
  --output reports/variant_calling_report.html \
  --template templates/comprehensive_template.Rmd
```

### 2) Estatísticas básicas de variantes
```bash
python analysis/variant_stats.py \
  --vcf results/variants.vcf.gz \
  --output-json stats/variant_stats.json \
  --output-txt stats/variant_summary.txt
```

### 3) Análise comparativa entre callers
```bash
Rscript analysis/comparison_plots.R \
  --gatk results/variants_gatk.vcf.gz \
  --freebayes results/variants_freebayes.vcf.gz \
  --bcftools results/variants_bcftools.vcf.gz \
  --sample NA12878 \
  --output plots/caller_comparison/
```

### 4) Benchmarking detalhado
```bash
python analysis/benchmark_analysis.py \
  --truth data/truth_sets/NA12878.vcf.gz \
  --callers results/gatk.vcf.gz,results/freebayes.vcf.gz \
  --confident-regions data/truth_sets/confident_regions.bed \
  --output benchmark/detailed_analysis/
```

### 5) Dashboard interativo de QC
```bash
Rscript analysis/qc_dashboard.R \
  --data-dir results/variant_calling/ \
  --port 3838 \
  --host 0.0.0.0
```

### 6) Visualizações genômicas avançadas
```bash
python analysis/variant_visualization.py \
  --vcf results/variants.vcf.gz \
  --reference data/reference/genome.fa \
  --output-dir visualization/ \
  --plots manhattan,circos,heatmap
```

## 🔌 Integração com Outputs do Pipeline

- **Entrada padrão**: outputs de utils/ (VCFs normalizados, métricas calculadas, anotações processadas)
- **Integração com scripts/**: cada script de análise aceita --config para YAMLs de config/
- **Workflows**: análises podem ser incorporadas como etapas finais em workflows/nextflow/ e workflows/snakemake/
- **Outputs para visualization/**: plots e dados processados são exportados para o módulo de visualização
- **Formato de saída padronizado**: JSON/TSV/HTML para compatibilidade com downstream

## 📊 Dicas de Plots e Reporting

### Visualizações Essenciais
1. **Distribuições de Qualidade**: QUAL, DP, GQ por caller e tipo de variante
2. **Métricas de Concordância**: Venn diagrams, upset plots para intersecções
3. **Performance**: ROC curves, precision-recall, F1-scores
4. **Genômicas**: Manhattan plots para GWAS-style, circos para visão global
5. **Controle de Qualidade**: Ti/Tv ratios, distribuições de AF, depth coverage

### Boas Práticas de Reporting
- Templates R Markdown padronizados (HTML, PDF, Word)
- Seções organizadas: Sumário Executivo, Métodos, Resultados, QC, Apêndices
- Gráficos interativos com plotly/DT para exploração
- Integração de tabelas descritivas com kable/kableExtra
- Versionamento de resultados e rastreabilidade de parâmetros

### Customização Visual
- Paletas de cores consistentes (RColorBrewer, viridis)
- Temas ggplot2 padronizados para uniformidade
- Logos e branding institucional em relatórios
- Exportação em múltiplos formatos (PNG, SVG, PDF)

## 🔗 Links para Workflows, Scripts e Utils

### Relacionados no Módulo
- **Scripts**: ../scripts/README.md (scripts executáveis de variant calling)
- **Utils**: ../utils/README.md (utilitários Python/R para processamento)
- **Workflows**: ../workflows/README.md (definições Nextflow/Snakemake/CWL)
- **Config**: ../config/README.md (arquivos de configuração YAML)

### Integração com Outros Módulos
- **Preprocessing**: ../../preprocessing/README.md (QC inicial, trimming)
- **Alignment**: ../../alignment/README.md (alinhamento e processamento BAM)
- **Annotation**: ../../annotation/README.md (anotação funcional de variantes)
- **Visualization**: ../../visualization/README.md (visualizações genômicas avançadas)

### Documentação Principal
- **README do Pipeline**: ../../../README.md (documentação completa)
- **Workflows Globais**: ../../../workflows/README.md (workflows do pipeline completo)
- **Configurações Globais**: ../../../config/README.md (configurações centralizadas)

## 🧩 Requisitos e Dependências

### Linguagens e Ambientes
- **R** >= 4.1 com Bioconductor >= 3.14
- **Python** >= 3.9 com pandas, matplotlib, seaborn, plotly
- **Pandoc** >= 2.11 para renderização de relatórios
- **LaTeX** (opcional) para PDFs de alta qualidade

### Pacotes R Essenciais
```r
# Bioconductor
VariantAnnotation, GenomicRanges, Biostrings, rtracklayer
# CRAN
ggplot2, plotly, DT, kableExtra, rmarkdown, shiny, shinydashboard
```

### Pacotes Python Essenciais
```python
# Processamento
pandas, numpy, scipy, pysam
# Visualização
matplotlib, seaborn, plotly, bokeh
# Genômica
pyvcf, cyvcf2, pybedtools
```

## 📦 Instalação e Setup Rápido

### Ambiente R
```bash
# Instalar dependências R
Rscript -e "install.packages(c('ggplot2', 'plotly', 'DT', 'rmarkdown', 'shiny'))"
Rscript -e "BiocManager::install(c('VariantAnnotation', 'GenomicRanges'))"
```

### Ambiente Python
```bash
# Instalar dependências Python
pip install pandas matplotlib seaborn plotly cyvcf2 pybedtools

# Ou via requirements.txt
pip install -r requirements.txt
```

### Verificação de Setup
```bash
# Testar funcionalidades principais
Rscript -e "library(VariantAnnotation); print('R setup OK')"
python -c "import pandas, matplotlib, cyvcf2; print('Python setup OK')"
```

## 🔒 Boas Práticas de Análise

### Controle de Qualidade
- Validar inputs antes do processamento (formato VCF, índices)
- Implementar checkpoints em análises longas
- Logs estruturados para rastreamento de erros
- Backup automático de resultados críticos

### Reprodutibilidade
- Versionamento de scripts e templates
- Seeds fixos para análises estocásticas
- Containerização com renv/conda para R/Python
- Documentação de versões de ferramentas e bancos

### Escalabilidade
- Paralelização com furrr (R) e multiprocessing (Python)
- Processamento por chunks para datasets grandes
- Otimização de memória para análises genome-wide
- Integração com sistemas de grid/cluster

## 📞 Suporte e Contribuição

- **Issues**: reportar bugs ou solicitar novas funcionalidades
- **Pull Requests**: contribuições com exemplos de entrada/saída
- **Templates**: manter templates R Markdown atualizados
- **Documentação**: atualizar este README ao adicionar novos scripts

---

**Última atualização**: $(date '+%Y-%m-%d')  
**Versão do módulo**: 1.0.0  
**Compatibilidade**: Pipeline genomic-data-analysis-pipeline v2.0+
