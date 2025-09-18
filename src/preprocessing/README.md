# 🧬 Módulos de Pré-processamento e Controle de Qualidade

## 📋 Visão Geral

Este diretório contém módulos especializados para pré-processamento de dados de sequenciamento de próxima geração (NGS) e controle de qualidade abrangente. Os módulos implementam as etapas iniciais críticas do pipeline de análise genômica, garantindo que os dados brutos sejam adequadamente processados antes das análises subsequentes.

## 🎯 Objetivo

O objetivo principal desta pasta é fornecer ferramentas robustas e modulares para:

- **Controle de Qualidade**: Avaliação da qualidade dos dados brutos de sequenciamento
- **Filtragem de Qualidade**: Remoção de reads de baixa qualidade
- **Trimagem de Adaptadores**: Remoção de sequências de adaptadores e primers
- **Normalização**: Padronização de formatos de dados
- **Relatórios de QC**: Geração de relatórios detalhados de qualidade

## 📁 Estrutura de Conteúdo

### Módulos Principais

```
src/preprocessing/
├── quality_control/
│   ├── fastqc_wrapper.py      # Wrapper para FastQC
│   ├── multiqc_reporter.py    # Agregação de relatórios QC
│   └── quality_metrics.py     # Métricas customizadas de qualidade
├── trimming/
│   ├── adapter_trimming.py    # Trimagem de adaptadores
│   ├── quality_trimming.py    # Trimagem baseada em qualidade
│   └── primer_removal.py      # Remoção de primers
├── filtering/
│   ├── read_filtering.py      # Filtragem de reads
│   ├── duplicate_removal.py   # Remoção de duplicatas
│   └── contamination_check.py # Verificação de contaminação
└── normalization/
    ├── format_converter.py    # Conversão entre formatos
    ├── read_normalizer.py     # Normalização de reads
    └── batch_processor.py     # Processamento em lote
```

### Ferramentas Integradas

- **FastQC**: Controle de qualidade de sequências
- **MultiQC**: Agregação de relatórios
- **Trimmomatic**: Trimagem de adaptadores e qualidade
- **fastp**: Pré-processamento rápido de FASTQ
- **cutadapt**: Remoção de adaptadores
- **bbduk**: Filtragem e trimagem avançada

## 🔧 Funcionalidades

### Controle de Qualidade
- Análise de qualidade base por posição
- Distribuição de qualidade por sequência
- Detecção de sequências sobrerrepresentadas
- Análise de conteúdo GC
- Detecção de duplicatas
- Verificação de adaptadores

### Pré-processamento
- Trimagem adaptativa de qualidade
- Remoção de adaptadores automática
- Filtragem de reads curtas
- Correção de erros de sequenciamento
- Normalização de profundidade de cobertura

### Relatórios
- Dashboards interativos de QC
- Métricas estatísticas detalhadas
- Gráficos de distribuição de qualidade
- Comparações antes/depois do processamento
- Alertas automáticos para problemas de qualidade

## 📊 Tipos de Dados Suportados

- **DNA-seq**: Sequenciamento de DNA genômico
- **RNA-seq**: Sequenciamento de RNA (bulk e single-cell)
- **ChIP-seq**: Sequenciamento de imunoprecipitação
- **ATAC-seq**: Sequenciamento de acessibilidade da cromatina
- **Bisulfite-seq**: Sequenciamento de metilação
- **Amplicon-seq**: Sequenciamento direcionado

## 🛠️ Tecnologias Utilizadas

### Linguagens de Programação
- **Python**: Processamento de dados e automação
- **R**: Análise estatística e visualização
- **Bash**: Scripts de automação

### Bibliotecas e Frameworks
- **BioPython**: Manipulação de sequências biológicas
- **pysam**: Processamento de arquivos SAM/BAM
- **pandas**: Manipulação de dados
- **matplotlib/seaborn**: Visualização
- **Nextflow/Snakemake**: Integração de workflows

## 🚀 Como Usar

### Execução Básica
```bash
# Controle de qualidade básico
python quality_control/fastqc_wrapper.py --input raw_data/ --output qc_results/

# Trimagem de adaptadores
python trimming/adapter_trimming.py --input raw_reads.fastq --output trimmed_reads.fastq

# Filtragem de qualidade
python filtering/read_filtering.py --input trimmed_reads.fastq --min-quality 20
```

### Integração com Workflows
```bash
# Via Nextflow
nextflow run preprocess.nf --reads "data/*.fastq.gz" --outdir results/

# Via Snakemake
snakemake --configfile config/preprocess_config.yaml
```

## 📈 Métricas de Qualidade

### Métricas Principais
- **Q30**: Porcentagem de bases com qualidade ≥ 30
- **GC Content**: Distribuição de conteúdo GC
- **Duplication Rate**: Taxa de duplicação
- **Adapter Content**: Presença de adaptadores
- **N Content**: Porcentagem de bases indeterminadas

### Critérios de Aceitação
- Q30 > 80%
- Conteúdo GC entre 40-60% (dependendo da espécie)
- Taxa de duplicação < 20%
- Conteúdo de adaptadores < 5%

## 📋 Dependências

### Ferramentas Externas
```bash
# Instalação via conda
conda install -c bioconda fastqc multiqc trimmomatic fastp cutadapt
```

### Bibliotecas Python
```bash
pip install biopython pysam pandas matplotlib seaborn plotly
```

## 🧪 Testes

Todos os módulos incluem testes unitários abrangentes:

```bash
# Executar testes
pytest tests/test_preprocessing.py -v

# Testes de integração
pytest tests/integration/test_workflow.py
```

## 📚 Documentação

- Cada módulo possui documentação detalhada
- Exemplos de uso em `examples/`
- Tutoriais em `docs/tutorials/`
- API reference em `docs/api/`

## 🔗 Integração

Este módulo se integra com:
- `src/alignment/`: Dados pré-processados são enviados para alinhamento
- `src/visualization/`: Relatórios de QC são visualizados
- `config/`: Parâmetros de configuração
- `workflows/`: Integração com sistemas de workflow

## 🤝 Contribuição

Para contribuir com novos módulos de pré-processamento:
1. Siga as convenções de nomenclatura
2. Inclua testes unitários
3. Documente todas as funções
4. Adicione exemplos de uso
5. Atualize este README

## 📞 Suporte

Para questões específicas sobre pré-processamento:
- Consulte a documentação em `docs/preprocessing/`
- Verifique issues conhecidas no repositório
- Entre em contato com a equipe de desenvolvimento
