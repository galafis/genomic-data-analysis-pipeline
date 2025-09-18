# 🧬 Módulo de Anotação Funcional e Clínica de Variantes

## 📋 Descrição

Esta pasta reúne scripts e módulos para anotação funcional e clínica de variantes, integrando ferramentas como VEP, SnpEff, ANNOVAR e bancos de dados relevantes (dbSNP, ClinVar, gnomAD) para enriquecer os resultados gerados após a chamada de variantes.

## 🛠️ Ferramentas Principais

### Ferramentas de Anotação Funcional

• **VEP (Variant Effect Predictor)**: Predição de efeitos funcionais e consequências biológicas de variantes
• **SnpEff**: Anotação genômica e predição de efeitos de variantes genéticas
• **ANNOVAR**: Anotação funcional de variantes genéticas identificadas
• **ClinVar**: Anotação de significado clínico e patogenicidade
• **InterVar**: Interpretação automática de variantes baseada em diretrizes ACMG

### Bancos de Dados Integrados

• **dbSNP**: Base de dados de polimorfismos de nucleotídeo único
• **gnomAD**: Frequências alélicas populacionais globais
• **ClinVar**: Interpretações clínicas de variantes
• **COSMIC**: Catálogo de mutações somáticas em câncer
• **1000 Genomes**: Variações genéticas populacionais
• **ExAC**: Agregador de dados de exoma

## 📁 Estrutura de Diretórios

```
annotation/
├── functional/           # Anotação funcional
│   ├── vep_annotator.py
│   ├── snpeff_annotator.py
│   └── annovar_wrapper.py
├── clinical/            # Anotação clínica
│   ├── clinvar_annotator.py
│   ├── intervar_classifier.py
│   └── pathogenicity_predictor.py
├── population/          # Frequências populacionais
│   ├── gnomad_frequencies.py
│   ├── thousand_genomes.py
│   └── exac_frequencies.py
├── databases/           # Gestão de bancos de dados
│   ├── db_downloader.py
│   ├── db_updater.py
│   └── db_validator.py
├── utils/               # Utilitários auxiliares
│   ├── vcf_annotator.py
│   ├── tsv_converter.py
│   └── report_generator.py
└── workflows/           # Workflows integrados
    ├── comprehensive_annotation.nf
    ├── clinical_pipeline.nf
    └── population_annotation.nf
```

## 🎯 Tipos de Anotação Suportados

### Anotação Funcional
• **Consequência no Gene**: Impacto da variante na estrutura e função do gene
• **Predição de Patogenicidade**: Scores de predição (SIFT, PolyPhen, CADD)
• **Conservação Evolutiva**: Scores de conservação (PhyloP, PhastCons)
• **Domínios Proteicos**: Anotação de domínios funcionais afetados

### Anotação Clínica
• **Significado Clínico**: Classificação ACMG/AMP (Patogênica, Benigna, VUS)
• **Associações com Doenças**: Correlações fenotípicas conhecidas
• **Farmacogenômica**: Impacto na resposta a medicamentos
• **Penetrância e Expressividade**: Características de herança

### Anotação Populacional
• **Frequências Alélicas**: Distribuição em diferentes populações
• **Análise de Raridade**: Classificação baseada em frequência
• **Desequilíbrio de Ligação**: Padrões de herança em blocos
• **Estrutura Populacional**: Estratificação étnica e geográfica

## ⚙️ Fluxo de Trabalho Principal

### 1. Preparação dos Dados

```bash
# Validação do arquivo VCF de entrada
python utils/vcf_annotator.py --validate --input variants.vcf

# Normalização das variantes
bcftools norm -m-both -f reference.fasta variants.vcf > normalized.vcf
```

### 2. Anotação Funcional com VEP

```bash
# Anotação com VEP (Variant Effect Predictor)
python functional/vep_annotator.py \
    --input normalized.vcf \
    --cache /path/to/vep/cache \
    --fasta reference.fasta \
    --output vep_annotated.vcf
```

### 3. Anotação com SnpEff (Alternativa)

```bash
# Anotação com SnpEff
python functional/snpeff_annotator.py \
    --input normalized.vcf \
    --database GRCh38.99 \
    --output snpeff_annotated.vcf
```

### 4. Anotação de Frequências Populacionais

```bash
# Anotação com gnomAD
python population/gnomad_frequencies.py \
    --input vep_annotated.vcf \
    --database gnomad.v3.1.2 \
    --output gnomad_annotated.vcf
```

### 5. Anotação Clínica

```bash
# Anotação com ClinVar
python clinical/clinvar_annotator.py \
    --input gnomad_annotated.vcf \
    --database clinvar.vcf.gz \
    --output clinvar_annotated.vcf

# Classificação ACMG com InterVar
python clinical/intervar_classifier.py \
    --input clinvar_annotated.vcf \
    --output acmg_classified.vcf
```

### 6. Geração de Relatórios

```bash
# Geração de relatório resumido
python utils/report_generator.py \
    --input acmg_classified.vcf \
    --format html \
    --output annotation_report.html
```

## 📊 Campos de Anotação Principais

### Campos VEP
• **Consequence**: Consequência predita da variante
• **IMPACT**: Severidade do impacto (HIGH, MODERATE, LOW, MODIFIER)
• **SYMBOL**: Símbolo do gene
• **Gene**: ID Ensembl do gene
• **Feature_type**: Tipo de feature afetada (Transcript, RegulatoryFeature)
• **SIFT**: Score de predição SIFT
• **PolyPhen**: Score de predição PolyPhen

### Campos gnomAD
• **AF**: Frequência alélica global
• **AF_popmax**: Frequência máxima em qualquer população
• **AC**: Contagem de alelos
• **AN**: Número total de alelos
• **nhomalt**: Número de homozigotos alternativos

### Campos ClinVar
• **CLNSIG**: Significado clínico
• **CLNREVSTAT**: Status de revisão
• **CLNDN**: Nome da doença
• **CLNDISDB**: Base de dados da doença

## 🔧 Configurações Recomendadas

### Parâmetros VEP
```bash
--everything
--cache
--offline
--format vcf
--vcf
--symbol
--terms SO
--tsl
--hgvs
--fasta reference.fasta
--canonical
--protein
--biotype
--af
--af_1kg
--af_gnomad
--max_af
--pubmed
--variant_class
```

### Parâmetros SnpEff
```bash
-v          # Verbose
-stats      # Gerar estatísticas
-csvStats   # Estatísticas em CSV
-canon      # Usar apenas transcritos canônicos
```

### Filtros de Qualidade
```bash
# Critérios mínimos para anotação
QUAL >= 20.0
DP >= 8
GQ >= 15
AF >= 0.01  # Para variantes somáticas
```

## 🚀 Execução de Workflows

### Workflow Completo - Anotação Abrangente
```bash
nextflow run workflows/comprehensive_annotation.nf \
    --input 'variants/*.vcf' \
    --reference 'reference/genome.fasta' \
    --vep_cache '/opt/vep/cache' \
    --databases 'databases/' \
    --output 'results/annotation/'
```

### Workflow - Anotação Clínica Focada
```bash
nextflow run workflows/clinical_pipeline.nf \
    --input 'variants/*.vcf' \
    --clinvar 'databases/clinvar.vcf.gz' \
    --acmg_classification true \
    --output 'results/clinical/'
```

### Workflow - Anotação Populacional
```bash
nextflow run workflows/population_annotation.nf \
    --input 'variants/*.vcf' \
    --gnomad 'databases/gnomad.vcf.gz' \
    --thousand_genomes 'databases/1000g.vcf.gz' \
    --output 'results/population/'
```

## 📈 Métricas de Qualidade da Anotação

### Cobertura de Anotação
• **Gene Coverage**: Percentual de variantes com anotação gênica
• **Functional Impact**: Percentual com predição de impacto funcional
• **Clinical Annotation**: Percentual com anotação clínica disponível
• **Population Frequency**: Percentual com dados de frequência populacional

### Qualidade das Predições
• **SIFT Score Distribution**: Distribuição dos scores SIFT
• **PolyPhen Predictions**: Classificações PolyPhen (probably/possibly damaging)
• **CADD Score Range**: Faixa de scores CADD observados
• **Conservation Metrics**: Métricas de conservação evolutionary

## 🔍 Interpretação dos Resultados

### Classificação ACMG/AMP
• **Pathogenic (P)**: Evidência forte de patogenicidade
• **Likely Pathogenic (LP)**: Evidência moderada de patogenicidade
• **Uncertain Significance (VUS)**: Significado clínico incerto
• **Likely Benign (LB)**: Evidência moderada de benignidade
• **Benign (B)**: Evidência forte de benignidade

### Priorização de Variantes
• **High Impact**: Variantes com consequências severas (stop_gained, frameshift)
• **Rare Variants**: Frequência alélica < 0.1% em populações de referência
• **Known Pathogenic**: Variantes já descritas como patogênicas
• **Novel Variants**: Variantes não descritas em bases de dados públicas

## 📚 Recursos Adicionais

### Documentação
• [VEP User Guide](https://www.ensembl.org/info/docs/tools/vep/index.html)
• [SnpEff Manual](http://snpeff.sourceforge.net/SnpEff_manual.html)
• [ANNOVAR Documentation](https://annovar.openbioinformatics.org/en/latest/)
• [ClinVar Submission Guidelines](https://www.ncbi.nlm.nih.gov/clinvar/docs/submit/)

### Bases de Dados
• **Ensembl**: Anotações genômicas e transcritos
• **RefSeq**: Sequências de referência curadas
• **UniProt**: Informações sobre proteínas
• **OMIM**: Herança mendeliana online em humanos

## ⚠️ Considerações Importantes

### Limitações da Anotação
• **Completude**: Nem todas as variantes têm anotação funcional completa
• **Atualização**: Bases de dados requerem atualização regular
• **Contexto Clínico**: Anotação deve ser interpretada no contexto clínico
• **Validação**: Predições computacionais requerem validação experimental

### Boas Práticas
• **Múltiplas Ferramentas**: Use diferentes algoritmos de predição
• **Controle de Versão**: Documente versões de bases de dados utilizadas
• **Filtros Apropriados**: Aplique filtros adequados ao contexto de análise
• **Revisão Manual**: Revise manualmente variantes de interesse clínico
• **Atualização Regular**: Mantenha bases de dados atualizadas

## 🔄 Manutenção e Atualizações

### Atualização de Bases de Dados
```bash
# Download e atualização automática
python databases/db_updater.py \
    --databases gnomad,clinvar,dbsnp \
    --target_dir databases/ \
    --schedule monthly
```

### Validação de Integridade
```bash
# Verificação da integridade das bases de dados
python databases/db_validator.py \
    --database_dir databases/ \
    --check_integrity \
    --report validation_report.txt
```

## 🤝 Contribuição

Para contribuir com melhorias neste módulo:

1. **Adicione novas ferramentas de anotação**
2. **Integre bases de dados adicionais**
3. **Otimize workflows para diferentes casos de uso**
4. **Implemente novos algoritmos de priorização**
5. **Melhore a geração de relatórios e visualizações**

---

**Nota**: Este módulo está em constante evolução. Sempre verifique a documentação mais recente e as melhores práticas atualizadas para cada ferramenta e base de dados utilizadas.
