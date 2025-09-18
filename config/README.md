# Config

Este diretório contém arquivos de configuração para o pipeline de análise genômica.

## Estrutura

O diretório será organizado com:

• **Nextflow/**: Perfis e configurações específicas do Nextflow
• **Snakemake/**: Arquivos de configuração YAML para workflows Snakemake
• **CWL/**: Configurações para Common Workflow Language
• **HPC/**: Perfis para diferentes sistemas HPC (Slurm, PBS, SGE)
• **Cloud/**: Configurações para execução em nuvem (AWS, GCP, Azure)
• **Resources/**: Especificações de recursos computacionais

## Uso

Os arquivos de configuração permitem personalizar:

• Parâmetros de entrada e saída dos workflows
• Recursos computacionais (CPU, memória, tempo)
• Caminhos para ferramentas e bancos de dados
• Configurações específicas de ambiente
• Perfis de execução (local, HPC, nuvem)

## Exemplos de Configuração

### Nextflow

```groovy
profiles {
    standard {
        process.executor = 'local'
        process.cpus = 2
        process.memory = '4 GB'
    }
    
    slurm {
        process.executor = 'slurm'
        process.queue = 'genomics'
        process.time = '24.h'
    }
}
```

### Snakemake

```yaml
# config/dna_seq_config.yaml
samples: "data/samples"
reference: "data/reference/genome.fa"
output: "results/dna_seq"

resources:
  bwa_mem:
    cpus: 8
    mem_gb: 16
  gatk_hc:
    cpus: 4
    mem_gb: 32
```

## Arquivos de Configuração Principais

• **nextflow.config**: Configuração principal do Nextflow (já presente na raiz)
• **dna_seq_config.yaml**: Parâmetros para análise DNA-seq
• **rna_seq_config.yaml**: Parâmetros para análise RNA-seq
• **chip_seq_config.yaml**: Parâmetros para análise ChIP-seq
• **scrna_seq_config.yaml**: Parâmetros para análise single-cell
• **resources.config**: Especificações de recursos por processo

Este diretório faz parte da estrutura organizacional do pipeline de análise genômica.
