# 🧬 Módulo de Chamada de Variantes (Variant Calling)

## 📋 Descrição

Esta pasta abriga scripts e módulos para a chamada de variantes (SNPs, indels, CNVs, SVs) a partir de dados alinhados. Integra ferramentas como GATK, FreeBayes, Strelka2 e pipelines dedicados à detecção, filtragem e anotação de variantes.

## 🛠️ Ferramentas Principais

### Ferramentas de Chamada de Variantes
- **GATK HaplotypeCaller**: Chamada de variantes germinativas com alta sensibilidade
- **GATK Mutect2**: Detecção de variantes somáticas em amostras tumorais
- **FreeBayes**: Chamador de variantes bayesiano para múltiplos alelos
- **Strelka2**: Chamada de variantes somáticas e germinativas otimizada
- **VarScan**: Detecção de SNPs, indels e CNVs
- **Platypus**: Chamada de variantes com foco em regiões complexas

### Ferramentas de Processamento Estrutural
- **Manta**: Detecção de variantes estruturais (SVs)
- **LUMPY**: Identificação de rearranjos estruturais
- **BreakDancer**: Detecção de breakpoints e translocações
- **CNVnator**: Análise de variações no número de cópias (CNVs)

## 📁 Estrutura de Diretórios

```
variant_calling/
├── germline/              # Scripts para variantes germinativas
│   ├── gatk_haplotypecaller.py
│   ├── freebayes_caller.py
│   └── joint_genotyping.py
├── somatic/               # Scripts para variantes somáticas
│   ├── mutect2_caller.py
│   ├── strelka2_caller.py
│   └── varscan_somatic.py
├── structural/            # Detecção de variantes estruturais
│   ├── manta_sv.py
│   ├── lumpy_sv.py
│   └── cnv_analysis.py
├── filtering/             # Módulos de filtragem
│   ├── variant_filter.py
│   ├── quality_control.py
│   └── population_filter.py
├── annotation/            # Anotação de variantes
│   ├── vep_annotator.py
│   ├── snpeff_annotator.py
│   └── functional_impact.py
├── utils/                 # Utilitários auxiliares
│   ├── vcf_tools.py
│   ├── bed_operations.py
│   └── variant_stats.py
└── workflows/             # Workflows integrados
    ├── germline_pipeline.nf
    ├── somatic_pipeline.nf
    └── sv_pipeline.nf
```

## 🔬 Tipos de Variantes Suportadas

### Variantes de Pequena Escala
- **SNPs (Single Nucleotide Polymorphisms)**: Substituições de nucleotídeos únicos
- **Indels**: Inserções e deleções de pequeno porte (< 50bp)
- **MNPs (Multi-Nucleotide Polymorphisms)**: Substituições de múltiplos nucleotídeos

### Variantes Estruturais
- **CNVs (Copy Number Variants)**: Duplicações e deleções de segmentos
- **SVs (Structural Variants)**: Inversões, translocações, rearranjos
- **Repetições**: Expansões de sequências repetitivas

## ⚙️ Fluxo de Trabalho Principal

### 1. Preparação dos Dados
```bash
# Preparação de arquivos BAM
python prepare_bam.py --input aligned.bam --output ready.bam

# Criação de índices necessários
samtools index ready.bam
gatk CreateSequenceDictionary -R reference.fasta
```

### 2. Chamada de Variantes Germinativas
```bash
# GATK HaplotypeCaller
python gatk_haplotypecaller.py \
    --input ready.bam \
    --reference reference.fasta \
    --output variants.vcf

# FreeBayes (alternativa)
python freebayes_caller.py \
    --input ready.bam \
    --reference reference.fasta \
    --output variants_fb.vcf
```

### 3. Chamada de Variantes Somáticas
```bash
# Mutect2 para variantes somáticas
python mutect2_caller.py \
    --tumor tumor.bam \
    --normal normal.bam \
    --reference reference.fasta \
    --output somatic.vcf

# Strelka2 (alternativa otimizada)
python strelka2_caller.py \
    --tumor tumor.bam \
    --normal normal.bam \
    --reference reference.fasta \
    --output strelka_results/
```

### 4. Detecção de Variantes Estruturais
```bash
# Manta para SVs
python manta_sv.py \
    --input aligned.bam \
    --reference reference.fasta \
    --output manta_results/

# CNVnator para CNVs
python cnv_analysis.py \
    --input aligned.bam \
    --reference reference.fasta \
    --bin_size 100 \
    --output cnv_results.txt
```

### 5. Filtragem e Controle de Qualidade
```bash
# Filtragem por qualidade
python variant_filter.py \
    --input variants.vcf \
    --quality 30 \
    --depth 10 \
    --output filtered_variants.vcf

# Controle de qualidade
python quality_control.py \
    --input filtered_variants.vcf \
    --output qc_report.html
```

### 6. Anotação Funcional
```bash
# Anotação com VEP
python vep_annotator.py \
    --input filtered_variants.vcf \
    --database ensembl \
    --output annotated_variants.vcf

# Análise de impacto funcional
python functional_impact.py \
    --input annotated_variants.vcf \
    --output impact_analysis.txt
```

## 📊 Métricas de Qualidade

### Métricas Gerais
- **Ti/Tv Ratio**: Proporção transições/transversões (esperado: ~2.0-2.1)
- **Het/Hom Ratio**: Proporção heterozigoto/homozigoto
- **Singleton Rate**: Taxa de variantes únicas na população
- **Mendelian Errors**: Erros de herança mendeliana (trios)

### Métricas por Variante
- **QUAL**: Qualidade da chamada de variante
- **DP**: Profundidade de cobertura total
- **GQ**: Qualidade do genótipo
- **AF**: Frequência alélica
- **FS**: Fisher Strand Bias
- **MQ**: Qualidade de mapeamento

## 🔧 Configurações Recomendadas

### Parâmetros GATK HaplotypeCaller
```bash
--standard-min-confidence-threshold-for-calling 30.0
--standard-min-confidence-threshold-for-emitting 10.0
--max-alternate-alleles 3
--max-genotype-count 1024
```

### Parâmetros FreeBayes
```bash
--min-base-quality 20
--min-mapping-quality 30
--min-coverage 10
--min-alternate-fraction 0.2
```

### Critérios de Filtragem
```bash
# Hard filters para SNPs
QUAL >= 30.0
QD >= 2.0
FS <= 60.0
MQ >= 40.0
MQRankSum >= -12.5
ReadPosRankSum >= -8.0

# Hard filters para Indels
QUAL >= 30.0
QD >= 2.0
FS <= 200.0
ReadPosRankSum >= -20.0
```

## 🚀 Execução de Workflows

### Workflow Completo - Variantes Germinativas
```bash
nextflow run workflows/germline_pipeline.nf \
    --input 'data/bam/*.bam' \
    --reference 'reference/genome.fasta' \
    --known_sites 'reference/dbsnp.vcf' \
    --output 'results/germline/'
```

### Workflow Completo - Variantes Somáticas
```bash
nextflow run workflows/somatic_pipeline.nf \
    --tumor 'data/tumor/*.bam' \
    --normal 'data/normal/*.bam' \
    --reference 'reference/genome.fasta' \
    --panel_of_normals 'reference/pon.vcf' \
    --output 'results/somatic/'
```

### Workflow - Variantes Estruturais
```bash
nextflow run workflows/sv_pipeline.nf \
    --input 'data/bam/*.bam' \
    --reference 'reference/genome.fasta' \
    --exclude_regions 'reference/blacklist.bed' \
    --output 'results/structural/'
```

## 📈 Análises Downstream

### Análise de Associação (GWAS)
- Preparação de dados para análise de associação genômica
- Controle de qualidade populacional
- Correção de estratificação populacional
- Análises de associação caso-controle

### Análise de Segregação Familiar
- Análise de trios (pai-mãe-filho)
- Identificação de variantes de novo
- Análise de segregação co-dominante
- Identificação de variantes patogênicas

### Análise Farmacogenômica
- Identificação de variantes farmacogenômicas
- Predição de resposta a medicamentos
- Análise de metabolismo de drogas
- Relatórios clínicos personalizados

## 🔍 Validação e Benchmarking

### Conjuntos de Dados de Referência
- **Genome in a Bottle (GIAB)**: Padrões de referência para validação
- **1000 Genomes**: Variantes populacionais de referência
- **ExAC/gnomAD**: Frequências alélicas populacionais
- **ClinVar**: Variantes de significado clínico

### Métricas de Performance
- **Sensibilidade**: Taxa de variantes verdadeiras detectadas
- **Especificidade**: Taxa de não-variantes corretamente identificadas
- **Precisão**: Proporção de chamadas corretas
- **F1-Score**: Média harmônica entre precisão e sensibilidade

## 📚 Recursos Adicionais

### Documentação
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- [FreeBayes Manual](https://github.com/freebayes/freebayes)
- [Strelka2 Documentation](https://github.com/Illumina/strelka)
- [VEP User Guide](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Bases de Dados
- **dbSNP**: Variantes conhecidas e frequências
- **COSMIC**: Mutações somáticas em câncer
- **ClinVar**: Interpretação clínica de variantes
- **gnomAD**: Agregado de variantes populacionais

## ⚠️ Considerações Importantes

### Limitações Técnicas
- Regiões de baixa complexidade podem gerar falsos positivos
- Variantes em regiões repetitivas requerem validação adicional
- CNVs grandes podem não ser detectados adequadamente
- Mosaicismo pode afetar a detecção de variantes somáticas

### Boas Práticas
- Sempre use controles negativos e positivos
- Valide variantes críticas por métodos independentes
- Considere a qualidade da amostra e preparação da biblioteca
- Implemente controles de qualidade rigorosos
- Documente todos os parâmetros e versões de software utilizados

## 🤝 Contribuição

Para contribuir com melhorias neste módulo:
1. Adicione novos algoritmos de chamada de variantes
2. Otimize parâmetros para diferentes tipos de amostras
3. Implemente novos métodos de filtragem
4. Adicione suporte para novos formatos de dados
5. Melhore a documentação e exemplos de uso

---

**Nota**: Este módulo está em constante desenvolvimento. Sempre verifique a documentação mais recente e as melhores práticas atualizadas para cada ferramenta utilizada.
