# Configurações do Variant Calling

Este diretório contém arquivos de configuração em YAML e JSON que centralizam parâmetros do pipeline de chamada e anotação de variantes, bem como recursos computacionais. O objetivo é padronizar execuções entre scripts, workflows (Nextflow, Snakemake, CWL) e ambientes (local, cluster, cloud), mantendo reprodutibilidade e versionamento claro.

## Visão Geral

Os arquivos de configuração cobrem:
- Parâmetros de chamada de variantes (por ferramenta)
- Regras de filtragem e normalização
- Configuração de anotação funcional e bases de dados
- Recursos computacionais (threads, memória, filas/partições)
- Caminhos de entrada/saída e referências

Estrutura sugerida do diretório:

config/
├── gatk_config.yaml           # Parâmetros GATK HaplotypeCaller/GenotypeGVCFs
├── freebayes_config.yaml      # Parâmetros FreeBayes
├── bcftools_config.yaml       # Parâmetros bcftools mpileup/call
├── filtering_rules.yaml       # Regras de filtragem pós-chamada
├── annotation_config.yaml     # Parâmetros de anotação (VEP/SnpEff)
├── resources.yaml             # Recursos computacionais por etapa
├── references.yaml            # Caminhos para genoma e índices
└── profiles.json              # Perfis de ambiente (local, slurm, cloud)

## Exemplos de Blocos YAML

### 1) Parâmetros GATK (gatk_config.yaml)

caller:
  name: gatk
  haplotypecaller:
    emit_ref_confidence: GVCF
    standard_min_confidence_threshold_for_calling: 30.0
    dbsnp: data/db/dbsnp.vcf.gz
    intervals: null  # ex.: data/targets/panel.bed
    threads: 8
  genotypegvcfs:
    stand_call_conf: 30.0
    stand_emit_conf: 10.0
    max_alternate_alleles: 3

io:
  input_bam: null
  reference_fasta: data/reference/genome.fa
  output_vcf: results/variants/variants_gatk.vcf.gz

### 2) Parâmetros FreeBayes (freebayes_config.yaml)

caller:
  name: freebayes
  ploidy: 2
  min_alternate_fraction: 0.2
  min_mapping_quality: 30
  min_base_quality: 20
  regions_bed: null

io:
  input_bam: null
  reference_fasta: data/reference/genome.fa
  output_vcf: results/variants/variants_freebayes.vcf.gz

### 3) Parâmetros bcftools (bcftools_config.yaml)

caller:
  name: bcftools
  mpileup:
    q: 20     # base quality
    Q: 30     # mapping quality
    d: 10000  # max depth
  call:
    m: 1      # multiallelic calling
    v: true   # output variants only

io:
  input_bam: null
  reference_fasta: data/reference/genome.fa
  output_vcf: results/variants/variants_bcftools.vcf.gz

### 4) Filtragem (filtering_rules.yaml)

filters:
  snv:
    min_qual: 30
    min_dp: 10
    max_dp: 1000
    min_gq: 20
  indel:
    min_qual: 30
    min_dp: 10
    max_dp: 1000
  region_masks:
    low_complexity_bed: data/masks/low_complexity.bed
    blacklist_bed: data/masks/blacklist.bed
  normalization:
    left_align: true
    split_multiallelics: true

### 5) Anotação (annotation_config.yaml)

annotation:
  tool: vep  # vep | snpeff
  vep:
    cache_dir: /path/to/vep_cache
    assembly: GRCh38
    plugins:
      - LoF
      - CADD
    extra_args: "--everything --fork 4"
  snpeff:
    genome: GRCh38.99
    data_dir: /path/to/snpeff
  databases:
    gnomad: data/db/gnomad.genomes.vcf.gz
    clinvar: data/db/clinvar.vcf.gz
    dbnsfp: data/db/dbNSFP.tsv.gz

### 6) Recursos (resources.yaml)

resources:
  alignment:
    threads: 8
    mem_gb: 32
    time: "12:00:00"
    partition: normal
  variant_calling:
    threads: 16
    mem_gb: 64
    time: "24:00:00"
    partition: compute
  annotation:
    threads: 8
    mem_gb: 32
    time: "08:00:00"

### 7) Referências (references.yaml)

references:
  fasta: data/reference/genome.fa
  fai: data/reference/genome.fa.fai
  dict: data/reference/genome.dict
  known_sites:
    - data/db/known_sites/dbsnp.vcf.gz
    - data/db/known_sites/Mills_and_1000G.vcf.gz

### 8) Perfis (profiles.json)

{
  "local": {
    "executor": "local",
    "container": "docker",
    "work_dir": "work/local"
  },
  "slurm": {
    "executor": "slurm",
    "queue": "compute",
    "work_dir": "work/slurm"
  },
  "aws": {
    "executor": "aws",
    "instance_type": "m5.large",
    "work_dir": "s3://bucket/work"
  }
}

## Integração com Scripts e Workflows

- Scripts bash (scripts/): cada script aceita --config para apontar para o YAML correspondente.
  - Ex.: ./scripts/run_gatk_haplotypecaller.sh --config config/gatk_config.yaml
- Nextflow (workflows/nextflow): mapear params a partir de arquivos YAML/JSON via parsing no main.nf e perfis no nextflow.config.
  - Ex.: nextflow run workflows/nextflow/variant_calling.nf --config_path config/gatk_config.yaml -profile slurm
- Snakemake (workflows/snakemake): definir configfile: "config/variant_calling_config.yaml" e propagar chaves via config["..."] em regras.
  - Ex.: snakemake -s workflows/snakemake/Snakefile --configfile config/variant_calling_config.yaml --cores 16
- CWL (workflows/cwl): referenciar parâmetros em inputs.yml e dividir ferramentas em tools/ com argumentos vindos de annotation_config.yaml e filtering_rules.yaml.
  - Ex.: cwltool workflows/cwl/workflow.cwl inputs.yml

## Boas Práticas de Versionamento

- Fixe versões de ferramentas e bancos (dbSNP, gnomAD, ClinVar) nos YAMLs.
- Use controle semântico nos arquivos (campo version: x.y.z) e mantenha CHANGELOG.
- Evite alterar valores diretamente em produção; crie variantes por ambiente: config/envs/{local,slurm,cloud}.
- Valide mudanças com --dry-run (Snakemake) ou -resume/timeline (Nextflow) e registre hashes de config em logs.
- Inclua checksums (sha256) de referências em references.yaml para auditoria.

## Convenções de Estilo (padrão visual/organizacional do pipeline)

- Títulos e seções em português, com ícones quando aplicável (📋, 🗂️, 🔧, 🔒).
- Estrutura de tópicos concisa, blocos de código curtos e comentados.
- Caminhos relativos sob data/, results/, config/, workflows/ conforme demais READMEs.
- Nomes de chaves em snake_case; valores booleanos explícitos; comentários indicando unidades.

## Links Relacionados

- Scripts README: ../scripts/README.md
- Workflows README: ../workflows/README.md
- README Principal do módulo: ../README.md
- Documentação do Pipeline (raiz): ../../../README.md

## Suporte e Contribuição

- Abra issues para propor novos parâmetros ou bancos.
- Submeta PRs com exemplos mínimos de inputs e resultados.
- Atualize este README ao adicionar novos arquivos de configuração.
