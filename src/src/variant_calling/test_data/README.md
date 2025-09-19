# 📊 Test Data - Variant Calling

## 📋 Propósito

Este diretório contém dados de teste essenciais para validação, integração contínua e controle de qualidade do módulo de variant calling. Os dados de teste garantem reprodutibilidade, detectam regressões e facilitam benchmarking comparativo entre diferentes algoritmos de chamada de variantes.

### 🎯 Finalidades dos Dados de Teste

- **Validação**: Verificar funcionamento correto de pipelines e algoritmos
- **Integração**: Testes automatizados em CI/CD (GitHub Actions, Jenkins)
- **QA/QC**: Controle de qualidade contínuo e detecção de regressões
- **Benchmarking**: Comparação de performance entre callers (GATK, FreeBayes, bcftools)
- **Desenvolvimento**: Testes rápidos durante implementação de novas features
- **Treinamento**: Dados para documentação e tutoriais

## 🗂️ Estrutura dos Dados de Teste

```
test_data/
├── small/                      # Datasets pequenos para testes rápidos
│   ├── chr20_subset.bam       # BAM pequeno (chr20:1-5000000)
│   ├── chr20_subset.bam.bai   # Índice do BAM
│   ├── reference_chr20.fa     # Referência chr20 subset
│   ├── reference_chr20.fa.fai # Índice FASTA
│   └── reference_chr20.dict   # Dicionário de sequência
├── medium/                     # Datasets médios para validação
│   ├── NA12878_chr20.bam      # BAM completo chr20
│   ├── NA12878_chr20.bam.bai
│   ├── reference_hg38_chr20.fa
│   └── truth_sets/
│       ├── NA12878_GIAB_chr20.vcf.gz    # Truth set GIAB
│       ├── NA12878_GIAB_chr20.vcf.gz.tbi
│       └── confident_regions_chr20.bed   # Regiões confiáveis
├── synthetic/                  # Dados sintéticos controlados
│   ├── simulated_reads.bam     # Reads simulados com variantes conhecidas
│   ├── known_variants.vcf.gz   # Variantes ground truth
│   └── simulation_params.json  # Parâmetros de simulação
└── expected_outputs/           # Resultados esperados para validação
    ├── gatk/
    │   ├── small_test_output.vcf.gz
    │   └── medium_test_output.vcf.gz
    ├── freebayes/
    └── bcftools/
```

## 📁 Tipos de Arquivos Esperados

### 1. **Arquivos de Entrada Primários**

#### BAM Files
- **small/chr20_subset.bam**: ~100MB, chr20:1-5M, ~30x coverage
- **medium/NA12878_chr20.bam**: ~2GB, chr20 completo, ~30x coverage
- **synthetic/simulated_reads.bam**: Dados sintéticos com SNPs/indels conhecidos

#### Referências FASTA
- **reference_chr20.fa**: Genoma referência hg38 chr20 subset
- **reference_hg38_chr20.fa**: hg38 chr20 completo
- Todos com índices (.fai) e dicionários (.dict) correspondentes

### 2. **Truth Sets e Validação**

#### Truth Sets VCF
- **NA12878_GIAB_chr20.vcf.gz**: Genome in a Bottle high-confidence calls
- **known_variants.vcf.gz**: Variantes simuladas ground truth
- **1000G_chr20_common.vcf.gz**: Variantes comuns 1000 Genomes

#### Regiões Confiáveis BED
- **confident_regions_chr20.bed**: Regiões GIAB alta confiança
- **callable_regions.bed**: Regiões com cobertura adequada
- **exclude_regions.bed**: Regiões problemáticas (repetições, gaps)

### 3. **Arquivos de Configuração**

#### Parâmetros de Teste
- **test_config.yaml**: Configurações para diferentes cenários
- **benchmark_params.json**: Parâmetros de benchmarking
- **simulation_params.json**: Setup de simulação de dados

### 4. **Outputs Esperados**

#### VCF Results
- Resultados esperados para cada caller em cenários padronizados
- Métricas de performance (runtime, memória, precisão/recall)
- Relatórios HTML de comparação

## 🏆 Datasets Recomendados/Referência

### Datasets Públicos Padrão

#### 1. **Genome in a Bottle (GIAB)**
```bash
# NA12878 (CEU Female)
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/...

# Truth sets
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/...
```

#### 2. **1000 Genomes Project**
```bash
# Amostras diversas para teste de população
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/...
```

#### 3. **Platinum Genomes**
```bash
# High-quality Illumina calls
wget https://github.com/Illumina/PlatinumGenomes/...
```

### Datasets Sintéticos

#### 1. **Simulação com DWGSIM**
```bash
# Gerar reads com variantes conhecidas
dwgsim -N 1000000 -1 150 -2 150 reference.fa simulated
```

#### 2. **Simulação com ART**
```bash
# Simulação mais realista com modelos de erro
art_illumina -ss HS25 -i reference.fa -l 150 -f 30 -o simulated
```

## 🔄 Versionamento e Composição

### Estratégia de Versionamento

#### 1. **Semantic Versioning**
```
test_data_v2.1.0/
├── CHANGELOG.md     # Histórico de mudanças
├── VERSION.txt      # Versão atual
└── checksums.md5    # Verificação de integridade
```

#### 2. **Controle de Integridade**
```bash
# Gerar checksums para validação
find . -type f -name "*.bam" -o -name "*.vcf.gz" -o -name "*.fa" | \
  xargs md5sum > checksums.md5

# Verificar integridade
md5sum -c checksums.md5
```

### Composição Estratégica

#### 1. **Tamanhos Escalonados**
- **XS (< 50MB)**: Testes unitários rápidos (< 30s)
- **S (50MB-500MB)**: Testes de integração (< 5min)
- **M (500MB-2GB)**: Testes de validação (< 30min)
- **L (> 2GB)**: Benchmarking completo (horas)

#### 2. **Cenários de Cobertura**
- **Alta qualidade**: Q30+, DP 20-100x, mapeamento único
- **Baixa qualidade**: Q10-20, DP 5-15x, multi-mapping
- **Casos extremos**: Regiões repetitivas, homopolímeros, CNVs

## 🔧 Integração com Scripts/Workflows

### Scripts de Validação

#### 1. **validate_test_data.sh**
```bash
#!/bin/bash
# Validar integridade e completude dos dados de teste

./utils/validate_test_data.sh --input test_data/ --checksum --completeness
```

#### 2. **run_quick_test.sh**
```bash
#!/bin/bash
# Teste rápido com dados small/
./scripts/run_gatk_haplotypecaller.sh \
  --input test_data/small/chr20_subset.bam \
  --reference test_data/small/reference_chr20.fa \
  --output results/quick_test.vcf.gz \
  --test-mode
```

### Integração CI/CD

#### 1. **GitHub Actions**
```yaml
# .github/workflows/test_variant_calling.yml
name: Test Variant Calling
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Download test data
        run: ./test_data/download_test_data.sh
      - name: Run tests
        run: ./test_data/run_integration_tests.sh
```

#### 2. **Nextflow Testing**
```bash
# Teste com perfil específico
nextflow run workflows/variant_calling.nf \
  -profile test \
  --input 'test_data/small/*.bam' \
  --reference test_data/small/reference_chr20.fa
```

### Integração com Workflows

#### 1. **Test Profile Configuration**
```yaml
# config/test.config
params {
  input_dir = 'test_data/small/'
  reference = 'test_data/small/reference_chr20.fa'
  truth_set = 'test_data/small/truth_chr20.vcf.gz'
  output_dir = 'results/test/'
  quick_mode = true
}
```

## ✅ Boas Práticas

### 1. **Organização e Estrutura**
- Separar dados por tamanho e propósito (small/medium/large)
- Incluir sempre índices e arquivos auxiliares necessários
- Documentar origem, processamento e características de cada dataset
- Manter estrutura de diretórios consistente com pipeline principal

### 2. **Qualidade e Validação**
- Validar integridade com checksums MD5/SHA256
- Incluir metadados detalhados (coverage, qualidade, demographics)
- Testar todos os paths críticos do pipeline
- Documentar resultados esperados e tolerâncias aceitáveis

### 3. **Performance e Escalabilidade**
- Datasets pequenos para desenvolvimento iterativo (< 1min)
- Datasets médios para validação completa (< 30min)
- Automatizar download/setup de dados grandes
- Usar caching inteligente para evitar re-downloads

### 4. **Reprodutibilidade**
- Fixar seeds aleatórios em simulações
- Versionamento claro de ferramentas e referências
- Containerização de ambientes de teste
- Logs detalhados de geração e processamento

### 5. **Manutenção**
- Revisar e atualizar datasets regularmente
- Remover dados obsoletos ou redundantes
- Documentar mudanças em CHANGELOG.md
- Testar compatibilidade com versões de ferramentas

## 🧪 Casos de Uso Específicos

### 1. **Desenvolvimento Rápido**
```bash
# Teste rápido durante desenvolvimento
./test_data/quick_test.sh gatk
./test_data/quick_test.sh freebayes
./test_data/quick_test.sh bcftools
```

### 2. **Validação de Release**
```bash
# Suite completa de testes
./test_data/full_validation_suite.sh \
  --callers all \
  --datasets small,medium \
  --generate-report
```

### 3. **Benchmarking Comparativo**
```bash
# Comparação de performance
python test_data/benchmark_callers.py \
  --input test_data/medium/NA12878_chr20.bam \
  --truth test_data/medium/truth_sets/NA12878_GIAB_chr20.vcf.gz \
  --output benchmark_results/
```

## 🔗 Links Úteis

### Documentação Relacionada
- **Pipeline principal**: ../README.md
- **Scripts de execução**: ../scripts/README.md
- **Configurações**: ../config/README.md
- **Análises**: ../analysis/README.md
- **Workflows**: ../workflows/README.md
- **Troubleshooting**: ../troubleshooting/README.md

### Recursos Externos
- **GIAB**: https://www.nist.gov/programs-projects/genome-bottle
- **1000 Genomes**: http://www.internationalgenome.org/
- **Platinum Genomes**: https://github.com/Illumina/PlatinumGenomes
- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651
- **FreeBayes**: https://github.com/freebayes/freebayes
- **bcftools**: http://samtools.github.io/bcftools/

### Simuladores de Dados
- **DWGSIM**: https://github.com/nh13/DWGSIM
- **ART**: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/
- **NEAT**: https://github.com/zstephens/neat-genreads

## 📦 Requisitos e Dependências

### Espaço de Armazenamento
- **Mínimo**: 5GB (datasets small + expected outputs)
- **Recomendado**: 20GB (datasets small + medium + sintéticos)
- **Completo**: 50GB+ (todos os datasets + histórico de testes)

### Tempo de Download/Setup
- **Datasets pequenos**: < 10min (download direto)
- **Datasets médios**: 30min-2h (dependendo da conexão)
- **Setup inicial completo**: 2-4h (primeira vez)

### Ferramentas Necessárias
```bash
# Validação e processamento
samtools >= 1.15
bcftools >= 1.15  
tabix >= 1.15
bedtools >= 2.30
python >= 3.8
R >= 4.0

# Opcional para simulação
dwgsim
art_illumina
neat-genreads
```

---

**Última atualização**: 2025-09-19  
**Versão**: 1.0.0  
**Compatibilidade**: Pipeline genomic-data-analysis-pipeline v2.0+
