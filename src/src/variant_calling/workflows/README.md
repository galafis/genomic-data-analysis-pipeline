# Workflows

Este diretório contém os workflows de análise de variantes genômicas implementados com diferentes sistemas de gerenciamento de workflow.

## Visão Geral

Os workflows implementam o pipeline de variant calling usando três sistemas principais:
- **Nextflow**: Workflow moderno e escalável com suporte nativo a containers
- **Snakemake**: Sistema baseado em Python com sintaxe intuitiva
- **CWL (Common Workflow Language)**: Padrão aberto para portabilidade entre plataformas

## Estrutura de Arquivos

```
workflows/
├── nextflow/
│   ├── main.nf              # Workflow principal
│   ├── nextflow.config      # Configurações
│   └── modules/             # Módulos reutilizáveis
├── snakemake/
│   ├── Snakefile           # Workflow principal
│   ├── config.yaml         # Configurações
│   └── rules/              # Regras individuais
├── cwl/
│   ├── workflow.cwl        # Workflow principal
│   ├── tools/              # Ferramentas CWL
│   └── workflows/          # Sub-workflows
└── README.md               # Este arquivo
```

## Nextflow

### Descrição
Workflow moderno que oferece paralelização automática, suporte a containers e execução em diferentes ambientes (local, cluster, cloud).

### Exemplo Básico
```bash
# Execução local
nextflow run main.nf --input samples.csv --reference genome.fa --outdir results/

# Com perfil Docker
nextflow run main.nf --input samples.csv --reference genome.fa --outdir results/ -profile docker

# Em cluster SLURM
nextflow run main.nf --input samples.csv --reference genome.fa --outdir results/ -profile slurm
```

### Argumentos Principais
- `--input`: Arquivo CSV com amostras (obrigatório)
- `--reference`: Genoma de referência FASTA (obrigatório)
- `--outdir`: Diretório de saída (padrão: results/)
- `--skip_qc`: Pular controle de qualidade
- `--caller`: Caller de variantes (gatk, freebayes, varscan)
- `--intervals`: Arquivo BED com regiões de interesse

### Configuração
Edite `nextflow.config` para:
- Configurar recursos computacionais
- Definir perfis de execução
- Configurar containers/modules

### Integração
- Scripts: Referenciados via `params.scripts_dir`
- Configs: Carregados automaticamente do diretório config/
- Módulos: Importados do diretório modules/

## Snakemake

### Descrição
Sistema baseado em Python que usa sintaxe similar ao Make, oferecendo flexibilidade e facilidade de desenvolvimento.

### Exemplo Básico
```bash
# Execução local
snakemake --cores 8 --configfile config.yaml

# Com Conda
snakemake --cores 8 --configfile config.yaml --use-conda

# Em cluster
snakemake --cores 8 --configfile config.yaml --cluster "sbatch -p normal -t 24:00:00"
```

### Argumentos Principais
- `--cores`: Número de núcleos (obrigatório)
- `--configfile`: Arquivo de configuração
- `--use-conda`: Usar ambientes Conda
- `--use-singularity`: Usar containers Singularity
- `--dry-run`: Simular execução
- `--forceall`: Forçar re-execução de todas as regras

### Configuração
Edite `config.yaml` para:
- Definir caminhos de entrada e saída
- Configurar parâmetros dos tools
- Especificar recursos por regra

### Integração
- Scripts: Chamados via `script:` ou `shell:`
- Configs: Carregados via `configfile:`
- Módulos: Importados via `include:`

## CWL (Common Workflow Language)

### Descrição
Padrão aberto que garante portabilidade entre diferentes plataformas e executores CWL.

### Exemplo Básico
```bash
# Com cwltool
cwltool workflow.cwl inputs.yml

# Com Toil
toil-cwl-runner workflow.cwl inputs.yml

# Com Cromwell
java -jar cromwell.jar run workflow.cwl --inputs inputs.yml
```

### Argumentos Principais
Definidos no arquivo `inputs.yml`:
- `input_samples`: Array de arquivos FASTQ
- `reference_genome`: Arquivo FASTA de referência
- `output_directory`: Diretório de saída
- `variant_caller`: Ferramenta de variant calling

### Configuração
- `workflow.cwl`: Define o workflow principal
- `inputs.yml`: Especifica dados de entrada
- `tools/`: Contém definições CWL das ferramentas

### Integração
- Scripts: Referenciados via `baseCommand` e `arguments`
- Configs: Incorporados nos arquivos CWL ou inputs.yml
- Módulos: Importados via `$import` ou `run`

## Recomendações de Execução

### Execução Local
```bash
# Para datasets pequenos (<10 amostras)
# Nextflow
nextflow run main.nf --input small_dataset.csv --reference genome.fa -profile local

# Snakemake
snakemake --cores 4 --configfile config_local.yaml

# CWL
cwltool workflow.cwl inputs_small.yml
```

### Execução em Cluster
```bash
# SLURM com Nextflow
nextflow run main.nf --input dataset.csv --reference genome.fa -profile slurm

# SLURM com Snakemake
snakemake --cores 100 --cluster "sbatch -p compute -t 24:00:00 -c {threads} --mem={resources.mem_mb}"

# SLURM com CWL (usando Toil)
toil-cwl-runner --batchSystem slurm workflow.cwl inputs.yml
```

### Execução em Cloud
```bash
# AWS com Nextflow
nextflow run main.nf --input s3://bucket/samples.csv --reference s3://bucket/genome.fa -profile aws

# Google Cloud com Snakemake
snakemake --cores 50 --kubernetes --default-remote-provider GS --default-remote-prefix gs://bucket/

# Cloud com CWL (usando Toil)
toil-cwl-runner --provisioner aws --nodeType m5.large workflow.cwl inputs_cloud.yml
```

## Considerações de Performance

### Recursos Recomendados
- **CPU**: 8-32 cores por amostra para variant calling
- **Memória**: 16-64 GB dependendo do tamanho do genoma
- **Storage**: SSD para I/O intensivo, especialmente para alignment
- **Rede**: Baixa latência para sistemas distribuídos

### Otimizações
- Use containers/conda para reprodutibilidade
- Implemente checkpoints para workflows longos
- Configure cache para evitar re-execuções desnecessárias
- Use paralelização por cromossomo quando possível

## Monitoramento e Logs

### Nextflow
- Logs: `.nextflow.log`
- Reports: `timeline.html`, `report.html`
- Trace: `trace.txt`

### Snakemake
- Logs: `logs/` directory
- Reports: `snakemake --report report.html`
- Benchmarks: `benchmarks/` directory

### CWL
- Logs: Dependem do executor (cwltool, Toil, etc.)
- Outputs: Definidos nos arquivos CWL

## Links Relacionados

- [Scripts README](../scripts/README.md) - Documentação dos scripts utilizados
- [Config README](../config/README.md) - Configurações e parâmetros
- [README Principal](../README.md) - Visão geral do módulo variant_calling
- [Documentação do Pipeline](../../../../README.md) - Documentação geral do projeto

## Suporte e Contribuição

Para questões específicas sobre workflows:
1. Consulte a documentação oficial de cada sistema
2. Verifique os logs de execução
3. Teste em ambiente menor antes de execuções em produção
4. Contribua com melhorias via pull requests

---

**Nota**: Este README deve ser atualizado conforme novos workflows são adicionados ou modificações são feitas nos existentes.
