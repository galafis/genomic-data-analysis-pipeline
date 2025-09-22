# 🧬 Pré-processamento de Dados (QC, Trimming, Normalização, Controle de Lote)

## 📋 Visão Geral
Este módulo prepara dados de NGS para as etapas posteriores do pipeline. Abrange controle de qualidade (QC), remoção de adaptadores/trimagem, normalização inicial e verificações/mitigações de efeito de lote (batch), gerando saídas consistentes e reprodutíveis para alinhamento, chamada de variantes e análises downstream.

## 🎯 Funções do Módulo
- Controle de Qualidade (QC): FastQC/MultiQC, métricas customizadas, relatórios integrados.
- Trimagem/Adaptadores: Remoção de adaptadores/primers, trimagem por qualidade e tamanho mínimo.
- Filtragem Inicial: Remoção de reads de baixa qualidade, contaminantes, duplicatas técnicas.
- Normalização Inicial: Padronização de formatos, renomeação de amostras, layout de diretórios, checagens de consistência.
- Controle de Lote (Batch): Detecção de vieses por run/plataforma/lote e estratégias de mitigação inicial.

## 📁 Estrutura de Diretórios
```
src/preprocessing/
├── quality_control/
│   ├── fastqc_wrapper.py
│   ├── multiqc_reporter.py
│   └── quality_metrics.py
├── trimming/
│   ├── adapter_trimming.py
│   ├── quality_trimming.py
│   └── primer_removal.py
├── filtering/
│   ├── read_filtering.py
│   ├── duplicate_removal.py
│   └── contamination_check.py
└── normalization/
    ├── format_converter.py
    ├── read_normalizer.py
    └── batch_processor.py
```

## 🔧 Ferramentas Recomendadas e Pipelines
- FastQC: QC por amostra (relatórios HTML/ZIP).
- MultiQC: Agregação de relatórios QC em nível de projeto.
- fastp: Trimagem rápida (qualidade/adaptadores), relatórios JSON/HTML.
- Trimmomatic: Trimagem flexível com parâmetros finos.
- cutadapt: Remoção de adaptadores/primers personalizados.
- bbduk (BBTools): Filtragem/trimagem avançadas e remoção de contaminantes.

### Exemplos de Execução (scripts)
- QC básico (FastQC + MultiQC):
```
# QC por amostra
fastqc -t 8 -o qc/ fastq/*.fastq.gz

# Agregar
multiqc -o qc/multiqc/ qc/
```
- Trimagem com fastp (pareado):
```
fastp \
  -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
  -o trimmed/sample_R1.trim.fastq.gz -O trimmed/sample_R2.trim.fastq.gz \
  --detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 \
  --length_required 50 --thread 8 \
  --json qc/fastp/sample.json --html qc/fastp/sample.html
```
- Trimmomatic (pareado) com remoção de adaptadores:
```
trimmomatic PE -threads 8 \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  trimmed/R1_paired.fq.gz trimmed/R1_unpaired.fq.gz \
  trimmed/R2_paired.fq.gz trimmed/R2_unpaired.fq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
```
- Cutadapt para primers específicos (amplicon):
```
cutadapt -g ^PRIMER_FWD -G ^PRIMER_REV \
  -o trimmed/R1.fastq.gz -p trimmed/R2.fastq.gz \
  sample_R1.fastq.gz sample_R2.fastq.gz
```
- Agregação final de QC:
```
multiqc -o qc/multiqc_final/ qc/ trimmed/
```

## ⚙️ Orientações de Uso e Configuração
- Instalação via conda (recomendado):
```
conda install -c bioconda fastqc multiqc fastp trimmomatic cutadapt bbmap
```
- Parâmetros sugeridos:
  - Qualidade mínima (Phred): 20–30
  - Tamanho mínimo de read: 50–75 bp (ajuste conforme aplicação)
  - Remoção de adaptadores: habilitar detecção automática ou fornecer FASTA de adaptadores
  - Threads: utilize de acordo com o ambiente (ex.: 8–32)
- Layout de entradas/saídas:
  - Input: data/raw/<sample>_R{1,2}.fastq.gz
  - Output: data/processed/trimmed/, qc/, qc/multiqc/
- Validação pós-processamento:
  - Conferir taxas de sobrevivência de reads (>70% recomendável)
  - Verificar redução/eliminação de adaptadores
  - Checar distribuição de qualidade/GC após trimagem

## 🔗 Integração com Alinhamento e Workflows
- Alinhamento (src/alignment/): utilizar arquivos trimados pareados como entrada do alinhador (ex.: BWA-MEM2, STAR).
  - Ex.: `bwa-mem2 index ref.fa` e `bwa-mem2 mem -t 16 ref.fa R1_paired.fq.gz R2_paired.fq.gz | samtools sort -o sample.bam`
- Workflows: integração por Nextflow/Snakemake.
```
# Nextflow
nextflow run workflows/preprocess.nf \
  --reads "data/raw/*_{R1,R2}.fastq.gz" --outdir data/processed/

# Snakemake
snakemake -j 16 --configfile config/preprocess_config.yaml
```
- Passagem de metadados: garantir preservação de sample IDs/lanes para merge posterior.

## 🧪 Controle de Lote (Batch)
- Detecção inicial: comparar métricas de QC por run/máquina/kit (MultiQC sections por grupo).
- Mitigação no pré-processamento:
  - Homogeneizar parâmetros de trimagem entre lotes.
  - Remover contaminantes específicos do lote (bbduk/k-mer filtering).
- Registros: manter tabela de metadados (sample, run_date, instrument, kit, operador) em `metadata/samples.tsv`.
- Observação: correções estatísticas avançadas (ex.: ComBat) ocorrem em módulos downstream (expressão/contagem), mas a consistência inicial é crítica aqui.

## ✅ Boas Práticas
- Fixar versões de ferramentas (environment.yml/lockfiles).
- Registrar command-lines e logs por amostra (reprodutibilidade).
- Usar controles (spike-ins, phiX) e checar taxa de contaminação.
- Evitar over-trimming: balancear qualidade x cobertura.
- Reexecutar MultiQC após cada etapa crítica.
- Automatizar com testes em dados de exemplo (tests/test_preprocessing.py).

## 📈 Métricas-Chave de QC
- Q30 (≥80%), GC esperado por espécie, taxa de duplicação (<20%), conteúdo de adaptador (<5%), N content (<1%).
- fastp report JSON/HTML e FastQC: avaliar antes/depois do trimming.

## 🔌 Configuração (exemplo de YAML)
```yaml
preprocessing:
  threads: 16
  min_len: 50
  min_qual: 20
  adapter_fasta: adapters.fa
  tools:
    trim: fastp
    qc: [fastqc, multiqc]
  outputs:
    trimmed_dir: data/processed/trimmed
    qc_dir: qc/
```

## 🔗 Documentação Relacionada
- Main Pipeline README: src/src/README.md
- Variant Calling README: src/src/variant_calling/README.md
- Annotation README: src/src/variant_calling/annotation/README.md
- Filtering README: src/src/variant_calling/filtering/README.md
- Visualization README: src/src/variant_calling/visualization/README.md

## 🖼️ Padrão Visual do Projeto
- Cabeçalhos com emojis temáticos, listas claras, blocos de código com prompts mínimos, paleta consistente com demais READMEs.
- Estrutura: Visão Geral → Funções → Ferramentas/Exemplos → Configuração → Integração → Boas Práticas → Métricas → Links.

## 📞 Suporte
- Consulte docs/preprocessing/ e issues do repositório.
- Contato: equipe de desenvolvimento do projeto.
