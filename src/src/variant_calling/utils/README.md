# Utilitários do Variant Calling

Este diretório contém utilitários (Python e R) para análise auxiliar do módulo de variant calling: processamento e normalização de VCF, cálculo de métricas, parsing de anotações e rotinas de análise. O objetivo é padronizar operações comuns, reutilizáveis entre scripts e workflows (Nextflow, Snakemake, CWL), mantendo o mesmo padrão visual/organizacional do pipeline.

## 📋 Visão Geral

Componentes típicos (exemplos sugeridos):
- variant_caller.py — classes/helper para integração com GATK/FreeBayes/bcftools
- vcf_processor.py — normalização, filtragem e operações em VCF/BCF
- annotation_parser.py — parsing de anotações (VEP/SnpEff, dbNSFP, ClinVar, gnomAD)
- metrics_calculator.py — métricas de qualidade e concordância (PREC, REC, F1)
- run_benchmark.py — benchmarking sintético entre callers
- variant_stats.py — estatísticas descritivas de variantes

Estrutura recomendada:
utils/
├── variant_caller.py
├── vcf_processor.py
├── annotation_parser.py
├── metrics_calculator.py
├── run_benchmark.py
└── variant_stats.py

## 🔧 Principais Funções e Parâmetros

1) vcf_processor.py
- normalize(input_vcf, reference_fasta, output_vcf, left_align=True, split_multiallelics=True)
  - left_align: alinhamento à esquerda (bcftools norm)
  - split_multiallelics: divide variantes multialélicas
- filter_variants(input_vcf, output_vcf, min_qual=30, min_dp=10, max_dp=1000, min_gq=20, regions_bed=None, blacklist_bed=None)
- subset_regions(input_vcf, bed, output_vcf)

2) annotation_parser.py
- parse_vep(input_vcf, fields=["Consequence","SYMBOL","AF","CADD"], output_tsv=None)
- parse_snpeff(input_vcf, effects=["HIGH","MODERATE"], output_tsv=None)
- merge_databases(input_vcf, gnomad_vcf, clinvar_vcf, dbnsfp_tsv, output_vcf)

3) metrics_calculator.py
- confusion_matrix(truth_vcf, call_vcf, confident_bed=None, match_strict=True)
- compute_metrics(tp, fp, fn)
- compare_callers(callers_vcfs: dict, truth_vcf, out_report)

4) variant_caller.py
- run_gatk_hc(bam, ref, intervals=None, threads=8, out_vcf="variants_gatk.vcf.gz", extra_args=None)
- run_freebayes(bam, ref, regions_bed=None, ploidy=2, min_af=0.2, out_vcf="variants_freebayes.vcf.gz")
- run_bcftools(bam, ref, mpileup_opts=None, call_opts=None, out_vcf="variants_bcftools.vcf.gz")

5) run_benchmark.py
- main(args): orquestra rodada sintética/real e gera métricas (precisão, sensibilidade, F1, tempo)

6) variant_stats.py
- summarize(input_vcf, out_txt=None, out_json=None): contagens por tipo (SNV/INDEL), por impacto, AF, DP

## 🧪 Exemplos Sintéticos de Uso

1) Normalização e filtragem de VCF
```
python -c "from utils.vcf_processor import normalize, filter_variants; \
normalize('results/variants.vcf.gz','data/reference/genome.fa','results/variants.norm.vcf.gz'); \
filter_variants('results/variants.norm.vcf.gz','results/variants.filt.vcf.gz', min_qual=30, min_dp=10)"
```

2) Parsing de anotações (VEP)
```
python -c "from utils.annotation_parser import parse_vep; \
parse_vep('results/annotated.vcf.gz', fields=['Consequence','SYMBOL','CADD'], output_tsv='results/annot.tsv')"
```

3) Métricas contra truth set
```
python -c "from utils.metrics_calculator import confusion_matrix, compute_metrics; \
cm = confusion_matrix('data/truth_sets/NA12878.vcf.gz','results/variants.filt.vcf.gz', confident_bed='data/truth_sets/confident_regions.bed'); \
print(compute_metrics(**cm))"
```

4) Execução programática de callers
```
python -c "from utils.variant_caller import run_gatk_hc; \
run_gatk_hc('data/sample.bam','data/reference/genome.fa', threads=8, out_vcf='results/variants_gatk.vcf.gz')"
```

5) Benchmark sintético
```
python utils/run_benchmark.py \
  --callers gatk,freebayes,bcftools \
  --truth-set data/truth_sets/NA12878.vcf.gz \
  --confident-regions data/truth_sets/confident_regions.bed \
  --output results/benchmark/
```

6) Estatísticas rápidas
```
python utils/variant_stats.py \
  --input results/variants.filt.vcf.gz \
  --output results/stats/variant_statistics.txt
```

## 🔌 Integração com Scripts e Workflows

- Scripts (scripts/): cada script aceita --config para arquivos YAML em config/ e pode importar funções deste diretório (ex.: from utils.vcf_processor import normalize).
- Nextflow (workflows/nextflow): mapear params a partir de YAML/JSON e chamar módulos Python via process, com perfis definidos em nextflow.config.
  Ex.: nextflow run workflows/nextflow/variant_calling.nf --config_path config/gatk_config.yaml -profile slurm
- Snakemake (workflows/snakemake): usar configfile e importar funções via script/python no Snakefile (ou usar scripts: e conda:).
  Ex.: snakemake -s workflows/snakemake/Snakefile --configfile config/variant_calling_config.yaml --cores 16
- CWL (workflows/cwl): referenciar utilitários como CommandLineTool ou via dockerPull, com inputs de annotation_config.yaml e filtering_rules.yaml.

## 🔒 Boas Práticas

- Versões de ferramentas e bancos fixadas nos YAMLs de config/
- Nome de funções e parâmetros em snake_case; docstrings claros; logs no nível INFO/DEBUG
- Testes unitários para funções críticas (pytest)
- Entrada/saída comprimida quando possível (.vcf.gz/.bcf) e índices (.tbi/.csi)
- Checagem de integridade (sha256) de referências e bancos

## 🔗 Links Relacionados

- Scripts README: ../scripts/README.md
- Config README: ../config/README.md
- Workflows README: ../workflows/README.md
- README do módulo: ../README.md
- Documentação do Pipeline (raiz): ../../../README.md

## 🧩 Requisitos

- Python >= 3.9, bcftools >= 1.15, GATK >= 4.2, VEP/SnpEff conforme config/annotation_config.yaml
- R >= 4.1 para análises/relatórios quando aplicável

## 📦 Instalação rápida (opcional)

```
pip install -r requirements.txt
```

## 📞 Suporte e Contribuição

- Abra issues para novos utilitários ou melhorias
- Submeta PRs com exemplos mínimos de entrada/saída e testes
- Mantenha esta documentação atualizada ao adicionar/alterar funções
