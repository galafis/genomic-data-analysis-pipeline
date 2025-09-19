# README_batch.md — Execução Batch Multi-amostra

## Finalidade

Orientar processamento em lote (batch) de chamada, filtragem, anotação, QC e visualização de variantes em ambientes produtivos ou pesquisa.

## Organização Recomendada de Diretórios

project_root/
├── data/ # BAMs de entrada
│ ├── sample1.bam
│ └── sample2.bam
├── reference/ # FASTA, .fai, .dict
│ └── genome.fa
├── results/
│ ├── variants/
│ │ ├── sample1.vcf.gz
│ │ └── sample2.vcf.gz
│ ├── qc/
│ │ └── sample1_qc.txt
│ └── batch_summary/
│ └── batch_summary.txt


## Nomeação de Arquivos

- `{amostra}.bam`, `{amostra}.vcf.gz`, `{amostra}_qc.txt`, `{amostra}_report.txt`
- Sufixos especiais para: `_filtered`, `_annotated`, `_multi`, `_summary`, etc.

## Execução Batch Automatizada

1. Chame variantes com `batch_variant_calling.py` para todos os BAMs do diretório.
2. Filtre todos os VCFs com `variant_filter.py`.
3. Anote com `variant_annotation.py`.
4. Gere métricas com `qc_variant_stats.py`, relatórios com `vcf_report_generator.py`.
5. Sumarize resultados em `variant_batch_summary.py` e visualize em `vcf_visualization.py`.

### Exemplo Python

from batch_variant_calling import rodar_batch_variant_calling

bams = ["data/sample1.bam", "data/sample2.bam"]
resultados = rodar_batch_variant_calling(
lista_bams=bams,
referencia_fasta="reference/genome.fa",
caller="freebayes",
diretorio_saida="results/variants/",
threads=8
)
for log in resultados:
print(log)


## Requisitos Computacionais

- Para batch > 20 amostras: 16+ threads (CPU), 64GB RAM, preferir orquestradores de cluster (GNU Parallel, SLURM, PBS, Nextflow)
- Armazenar logs e histórico, preferencialmente com versionamento automático dos parâmetros

## Práticas Recomendadas e Troubleshooting

- Sempre faça QC dos BAMs antes do batch
- Log completo por amostra/processo facilita troubleshooting
- Verifique existência dos arquivos `.bai` e integridade do VCF antes da anotação
- Analise batch_summary antes de integração downstream/omics
- Documente bugs e lições aprendidas ao longo do reprocessamento




