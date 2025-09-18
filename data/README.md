# Data

Este diretório contém os dados utilizados pelo pipeline de análise genômica, incluindo dados de entrada, referências e dados de exemplo para testes e validação.

## Estrutura do Diretório

```
data/
├── samples/           # Dados de amostras de sequenciamento
│   ├── dna_seq/       # Amostras de DNA sequencing
│   ├── rna_seq/       # Amostras de RNA sequencing
│   ├── scrna_seq/     # Amostras de single-cell RNA sequencing
│   └── chip_seq/      # Amostras de ChIP sequencing
├── reference/         # Genomas e anotações de referência
│   ├── genome.fa      # Genoma de referência (FASTA)
│   ├── genes.gtf      # Anotação de genes (GTF/GFF)
│   └── indices/       # Índices pré-computados para alinhadores
├── databases/         # Bancos de dados para anotação
│   ├── dbsnp/         # Variantes conhecidas (dbSNP)
│   ├── cosmic/        # Mutações somáticas (COSMIC)
│   └── clinvar/       # Variantes clínicas (ClinVar)
└── examples/          # Dados de exemplo para testes
    ├── small_dataset/ # Dataset pequeno para testes rápidos
    └── tutorial/      # Dados para tutoriais
```

## Tipos de Dados Suportados

### Dados de Sequenciamento
- **FASTQ**: Arquivos de sequenciamento bruto (single-end e paired-end)
- **BAM/SAM**: Arquivos de alinhamento
- **VCF**: Arquivos de variantes
- **BED**: Arquivos de regiões genômicas
- **GTF/GFF**: Arquivos de anotação genômica

### Formatos de Entrada
- Illumina (HiSeq, NovaSeq, MiSeq)
- Oxford Nanopore (MinION, PromethION)
- PacBio (Sequel, RSII)
- 10x Genomics (Chromium)

## Organização de Amostras

### Convenção de Nomenclatura
```
{sample_id}_{condition}_{replicate}_{read}.fastq.gz
Exemplo: SAMPLE001_treatment_rep1_R1.fastq.gz
```

### Metadados de Amostras
Cada projeto deve incluir um arquivo `samples.tsv` com:
- sample_id: Identificador único da amostra
- condition: Condição experimental
- batch: Lote de sequenciamento
- library_type: Tipo de biblioteca (DNA-seq, RNA-seq, etc.)
- paired_end: Boolean indicando se é paired-end
- read_length: Comprimento das reads
- insert_size: Tamanho do fragmento (para paired-end)

## Genomas de Referência

### Humano
- **GRCh38**: Genoma de referência principal
- **GRCh37**: Versão anterior (para compatibilidade)
- **T2T-CHM13**: Genoma telômero-a-telômero completo

### Outros Organismos
- **Mouse (GRCm39)**: Para estudos modelo
- **Drosophila (dm6)**: Para genética funcional
- **C. elegans (ce11)**: Para estudos de desenvolvimento

## Bancos de Dados de Anotação

### Variantes
- **dbSNP**: Polimorfismos de nucleotídeo único
- **COSMIC**: Mutações somáticas em câncer
- **ClinVar**: Significado clínico de variantes
- **gnomAD**: Frequências populacionais

### Anotação Funcional
- **RefSeq**: Sequências de referência curadas
- **Ensembl**: Anotação genômica abrangente
- **GENCODE**: Anotação de genes codificadores e não-codificadores

## Requisitos de Armazenamento

### Espaço em Disco
- **Pequeno projeto** (< 10 amostras): ~100 GB
- **Projeto médio** (10-100 amostras): ~1 TB
- **Projeto grande** (> 100 amostras): > 10 TB

### Performance
- **SSD recomendado** para dados de trabalho ativo
- **Storage de rede** adequado para arquivamento
- **Backup regular** essencial para dados críticos

## Controle de Qualidade

### Verificações Automáticas
- Integridade de arquivos (checksums MD5/SHA256)
- Formato de arquivos válido
- Presença de metadados obrigatórios
- Convenções de nomenclatura

### Métricas de Qualidade
- Qualidade média das bases (Phred score)
- Duplication rate
- GC content
- Distribuição de comprimentos

## Segurança e Privacidade

### Dados Sensíveis
- **Dados humanos**: Requerem aprovação ética
- **Anonimização**: Obrigatória para dados clínicos
- **Criptografia**: Para dados em trânsito e em repouso

### Controle de Acesso
- Permissões baseadas em projeto
- Log de acesso auditável
- Políticas de retenção de dados

## Scripts Utilitários

### Preparação de Dados
```bash
# Validar integridade
./scripts/validate_data.sh data/samples/

# Gerar metadados
./scripts/generate_metadata.py data/samples/

# Verificar qualidade
./scripts/qc_check.sh data/samples/
```

### Download de Referências
```bash
# Genoma humano GRCh38
./scripts/download_reference.sh GRCh38

# Anotações Ensembl
./scripts/download_annotation.sh ensembl_homo_sapiens
```

## Exemplos de Uso

### Dataset de Teste
```bash
# Baixar dados de exemplo
wget -P data/examples/ \
  ftp://example.com/test_data.tar.gz

# Extrair e organizar
tar -xzf data/examples/test_data.tar.gz -C data/examples/
```

### Validação de Pipeline
```bash
# Executar com dados de exemplo
nextflow run main.nf \
  --input data/examples/small_dataset/ \
  --genome data/reference/genome.fa \
  --outdir results/test/
```

## Troubleshooting

### Problemas Comuns
- **Espaço insuficiente**: Verificar disponibilidade antes de downloads
- **Arquivos corrompidos**: Validar checksums após transferências
- **Permissões**: Verificar acesso de leitura/escrita nos diretórios
- **Formato incompatível**: Converter usando ferramentas apropriadas

### Contato para Suporte
Para questões relacionadas aos dados:
- Issues no GitHub: [github.com/galafis/genomic-data-analysis-pipeline/issues](https://github.com/galafis/genomic-data-analysis-pipeline/issues)
- Email: suporte-dados@projeto.org

---

**Nota**: Este diretório é parte integrante do pipeline de análise genômica e deve ser mantido organizado conforme as diretrizes estabelecidas. Todos os dados devem ser validados antes da utilização em análises de produção.
