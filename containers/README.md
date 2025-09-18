# Containers

Este diretório contém definições de containers Docker e Singularity para o pipeline de análise genômica.

## Estrutura

O diretório será organizado com:

• **Docker/**: Dockerfiles e configurações para containers Docker
• **Singularity/**: Definições de containers Singularity
• **Registry/**: Scripts para push/pull de imagens de registries públicos
• **Build/**: Scripts automatizados de construção de containers

## Uso

Os containers garantem reprodutibilidade e portabilidade do pipeline, encapsulando todas as dependências necessárias para cada etapa de análise.

### Docker

Para construir containers Docker localmente:

```bash
cd containers/docker
docker build -t genomic-pipeline:latest .
```

### Singularity

Para containers Singularity em ambientes HPC:

```bash
cd containers/singularity
singularity build genomic-pipeline.sif genomic-pipeline.def
```

## Containers Disponíveis

• **Base**: Container base com ferramentas fundamentais de bioinformática
• **Alignment**: Container especializado em alinhamento (BWA, Bowtie2, STAR)
• **Variant-calling**: Container para chamada de variantes (GATK, FreeBayes)
• **RNA-seq**: Container para análises de RNA-seq (Salmon, DESeq2)
• **Single-cell**: Container para análises single-cell (Seurat, Scanpy)
• **Visualization**: Container para visualizações (R Shiny, IGV)

Este diretório faz parte da estrutura organizacional do pipeline de análise genômica.
