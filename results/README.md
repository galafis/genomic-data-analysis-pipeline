# 📊 Resultados / Results

Este diretório é destinado ao armazenamento de todos os resultados gerados pelo pipeline de análise de dados genômicos.

## 📁 Estrutura dos Resultados

```
results/
├── dna_seq/              # Resultados de análise DNA-seq
│   ├── qc/               # Controle de qualidade
│   ├── alignment/        # Arquivos de alinhamento
│   ├── variants/         # Variantes identificadas
│   ├── annotation/       # Anotações funcionais
│   └── reports/          # Relatórios finais
├── rna_seq/              # Resultados de análise RNA-seq
│   ├── qc/               # Controle de qualidade
│   ├── alignment/        # Alinhamentos
│   ├── quantification/   # Quantificação de expressão
│   ├── differential/     # Análise diferencial
│   └── reports/          # Relatórios finais
├── scrna_seq/            # Resultados de single-cell RNA-seq
│   ├── qc/               # Controle de qualidade
│   ├── filtering/        # Filtragem de células/genes
│   ├── normalization/    # Normalização
│   ├── clustering/       # Clustering de células
│   ├── markers/          # Genes marcadores
│   └── reports/          # Relatórios finais
├── chip_seq/             # Resultados de ChIP-seq
│   ├── qc/               # Controle de qualidade
│   ├── alignment/        # Alinhamentos
│   ├── peaks/            # Identificação de picos
│   ├── motifs/           # Análise de motivos
│   └── reports/          # Relatórios finais
├── multi_omics/          # Análises integradas
│   ├── integration/      # Integração de dados
│   ├── machine_learning/ # Modelos ML
│   └── visualization/    # Visualizações integradas
└── README.md             # Este arquivo
```

## 📋 Tipos de Arquivos Gerados

### Arquivos de Controle de Qualidade
- **FastQC Reports**: `.html`, `.zip`
- **MultiQC Reports**: `multiqc_report.html`
- **Estatísticas de alinhamento**: `.txt`, `.log`

### Arquivos de Alinhamento
- **BAM/SAM files**: `.bam`, `.sam`
- **Índices**: `.bai`, `.csi`
- **Estatísticas**: `.flagstat`, `.idxstats`

### Arquivos de Variantes
- **VCF files**: `.vcf`, `.vcf.gz`
- **Anotações**: `.annotated.vcf`
- **Relatórios de impacto**: `.html`, `.txt`

### Arquivos de Expressão
- **Matrizes de contagem**: `.tsv`, `.csv`
- **Resultados de expressão diferencial**: `.xlsx`, `.csv`
- **Gráficos**: `.png`, `.pdf`, `.svg`

### Visualizações
- **Heatmaps**: `.png`, `.pdf`
- **PCA plots**: `.png`, `.pdf`
- **Volcano plots**: `.png`, `.pdf`
- **Dashboards interativos**: `.html`

## 🔍 Navegação por Tipo de Análise

### DNA-seq
- **Controle de Qualidade**: [results/dna_seq/qc/](dna_seq/qc/)
- **Variantes Identificadas**: [results/dna_seq/variants/](dna_seq/variants/)
- **Anotações Funcionais**: [results/dna_seq/annotation/](dna_seq/annotation/)
- **Relatórios Finais**: [results/dna_seq/reports/](dna_seq/reports/)

### RNA-seq
- **Quantificação**: [results/rna_seq/quantification/](rna_seq/quantification/)
- **Análise Diferencial**: [results/rna_seq/differential/](rna_seq/differential/)
- **Relatórios**: [results/rna_seq/reports/](rna_seq/reports/)

### Single-cell RNA-seq
- **Clustering**: [results/scrna_seq/clustering/](scrna_seq/clustering/)
- **Genes Marcadores**: [results/scrna_seq/markers/](scrna_seq/markers/)
- **Relatórios**: [results/scrna_seq/reports/](scrna_seq/reports/)

### ChIP-seq
- **Picos Identificados**: [results/chip_seq/peaks/](chip_seq/peaks/)
- **Motivos**: [results/chip_seq/motifs/](chip_seq/motifs/)
- **Relatórios**: [results/chip_seq/reports/](chip_seq/reports/)

## 📈 Relatórios Principais

Cada tipo de análise gera relatórios padronizados:

### Relatório de Controle de Qualidade
- **multiqc_report.html**: Relatório consolidado de QC
- **sample_statistics.txt**: Estatísticas por amostra
- **quality_metrics.csv**: Métricas de qualidade

### Relatório de Análise Principal
- **analysis_summary.html**: Sumário da análise
- **results_table.xlsx**: Tabela de resultados principais
- **visualization_dashboard.html**: Dashboard interativo

### Relatório Técnico
- **technical_report.pdf**: Relatório técnico detalhado
- **parameters.yaml**: Parâmetros utilizados
- **software_versions.txt**: Versões do software

## 🔒 Gestão de Arquivos

### Nomenclatura Padrão
```
{sample_id}_{analysis_type}_{step}_{timestamp}.{extension}

Exemplos:
- SAMPLE001_dna_seq_variants_20231218_1430.vcf
- SAMPLE002_rna_seq_counts_20231218_1430.tsv
- SAMPLE003_scrna_clustering_20231218_1430.pdf
```

### Organização Temporal
- **Por data**: Resultados organizados por data de execução
- **Por versão**: Versionamento de análises
- **Por projeto**: Separação por projeto/estudo

### Compactação
- Arquivos grandes (>100MB) são automaticamente compactados
- Formatos: `.gz`, `.bz2`, `.zip`
- Manter arquivos originais por 30 dias

## 📊 Monitoramento de Espaço

### Tamanhos Esperados
- **DNA-seq**: 5-20 GB por amostra
- **RNA-seq**: 1-5 GB por amostra
- **Single-cell**: 500MB-2GB por amostra
- **ChIP-seq**: 1-3 GB por amostra

### Limpeza Automática
- Arquivos temporários removidos após conclusão
- Logs antigos (>90 dias) arquivados
- Backups regulares dos resultados principais

## 🔧 Utilitários

### Scripts de Gestão
```bash
# Listar resultados por tipo
./scripts/list_results.sh --type dna_seq

# Gerar sumário de espaço
./scripts/storage_summary.sh

# Arquivar resultados antigos
./scripts/archive_results.sh --older-than 90

# Verificar integridade
./scripts/check_integrity.sh
```

### Visualização
```bash
# Abrir dashboard principal
./scripts/open_dashboard.sh

# Gerar relatório consolidado
./scripts/generate_report.sh --all

# Comparar resultados
./scripts/compare_runs.sh run1 run2
```

## 🔄 Backup e Restauração

### Estratégia de Backup
- **Diario**: Resultados do dia
- **Semanal**: Arquivos principais
- **Mensal**: Backup completo
- **Anual**: Arquivo permanente

### Localização dos Backups
- **Local**: `/backup/genomic_results/`
- **Rede**: `//server/backup/genomic/`
- **Nuvem**: `s3://genomic-backup/results/`

## 📞 Suporte

Para questões sobre os resultados:
- **Documentação**: [docs/user-guide/](../docs/user-guide/)
- **FAQ**: [docs/troubleshooting/faq.md](../docs/troubleshooting/faq.md)
- **Issues**: Abra uma issue no GitHub
- **Email**: suporte@genomic-pipeline.org

---

**Nota**: Este diretório é gerado automaticamente pelo pipeline. Não edite manualmente os arquivos de resultado.
