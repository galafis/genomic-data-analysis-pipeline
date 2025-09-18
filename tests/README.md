# Testes Automatizados - Pipeline de Análise Genômica

## 🧪 Visão Geral

Este diretório contém a estrutura de testes automatizados para o pipeline de análise de dados genômicos. Os testes são fundamentais para garantir a confiabilidade, reprodutibilidade e qualidade dos resultados das análises genômicas.

## 📋 Estrutura de Testes

```
tests/
├── unit/                    # Testes unitários para módulos individuais
│   ├── test_preprocessing/  # Testes para módulos de pré-processamento
│   ├── test_alignment/      # Testes para módulos de alinhamento
│   ├── test_variant_calling/# Testes para chamada de variantes
│   └── test_annotation/     # Testes para anotação de variantes
├── integration/             # Testes de integração de workflows
│   ├── test_dna_seq/       # Teste completo do pipeline DNA-seq
│   ├── test_rna_seq/       # Teste completo do pipeline RNA-seq
│   └── test_scrna_seq/     # Teste completo do pipeline single-cell
├── data/                    # Dados de teste sintéticos e referência
│   ├── reference/          # Genomas e anotações de referência reduzidos
│   ├── fastq/              # Arquivos FASTQ sintéticos
│   └── expected_outputs/    # Resultados esperados para validação
├── fixtures/                # Fixtures e configurações de teste
├── benchmarks/              # Testes de performance e benchmarking
└── README.md               # Este arquivo
```

## 🎯 Tipos de Testes

### Testes Unitários
- **Módulos de QC**: Verificação de métricas de qualidade
- **Funções de parsing**: Validação de parsers de formatos (SAM/BAM, VCF, GTF)
- **Algoritmos de processamento**: Testes de funções específicas
- **Utilities**: Testes de funções auxiliares e utilitários

### Testes de Integração
- **Workflows completos**: Execução end-to-end dos pipelines
- **Compatibilidade de formatos**: Verificação de inputs/outputs
- **Configurações múltiplas**: Teste com diferentes parâmetros
- **Ambientes diversos**: Validação em diferentes sistemas

### Testes de Regressão
- **Outputs determinísticos**: Verificação de resultados consistentes
- **Métricas de qualidade**: Validação de estatísticas esperadas
- **Formatos de saída**: Verificação de integridade dos arquivos

## 🚀 Execução dos Testes

### Pré-requisitos
```bash
# Instalar dependências de teste
conda install pytest pytest-cov pytest-xdist
pip install snakemake-test nextflow-test
```

### Testes Unitários
```bash
# Executar todos os testes unitários
pytest tests/unit/ -v

# Executar testes específicos
pytest tests/unit/test_preprocessing/ -v

# Executar com cobertura
pytest tests/unit/ --cov=src/ --cov-report=html
```

### Testes de Integração
```bash
# Teste do pipeline DNA-seq
nextflow run tests/integration/test_dna_seq/main.nf -profile test

# Teste do pipeline RNA-seq
snakemake -s tests/integration/test_rna_seq/Snakefile --configfile tests/fixtures/rna_config.yaml

# Teste do pipeline single-cell
pytest tests/integration/test_scrna_seq/ -v --slow
```

### Testes de Performance
```bash
# Benchmarking de performance
pytest tests/benchmarks/ --benchmark-only

# Profiling de memória
pytest tests/benchmarks/ --profile
```

## 📊 Métricas de Qualidade

### Cobertura de Código
- **Objetivo**: Mínimo 80% de cobertura
- **Crítico**: 95% para módulos de análise central
- **Relatórios**: Gerados automaticamente em HTML/XML

### Testes de Performance
- **Tempo de execução**: Monitoramento de regressões de performance
- **Uso de memória**: Validação de limites de RAM
- **Throughput**: Medição de amostras processadas por hora

### Validação de Resultados
- **Acurácia**: Comparação com truth sets conhecidos
- **Precisão**: Validação de métricas estatísticas
- **Reprodutibilidade**: Verificação de outputs determinísticos

## 🔧 Configuração de Ambiente de Teste

### Dados de Teste
```bash
# Download de dados de referência reduzidos
wget -O tests/data/reference/chr22.fa.gz https://example.com/chr22.fa.gz
gunzip tests/data/reference/chr22.fa.gz

# Geração de dados sintéticos
python scripts/generate_test_data.py --output tests/data/fastq/
```

### Configuração de CI/CD
```yaml
# .github/workflows/tests.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup conda
        uses: conda-incubator/setup-miniconda@v2
      - name: Install dependencies
        run: conda env create -f environment.yml
      - name: Run tests
        run: |
          conda activate genomic-pipeline
          pytest tests/ -v --cov=src/
```

## 🎮 Boas Práticas para Testes

### Nomenclatura
- **Arquivos de teste**: `test_*.py` ou `*_test.py`
- **Funções de teste**: `test_funcionalidade_condicao()`
- **Classes de teste**: `TestModuleName`
- **Fixtures**: Nome descritivo da funcionalidade

### Estrutura de Testes
```python
# Exemplo de teste unitário
import pytest
from src.preprocessing import quality_control

class TestQualityControl:
    @pytest.fixture
    def sample_fastq(self):
        return "tests/data/fastq/sample_R1.fastq.gz"
    
    def test_fastqc_basic_stats(self, sample_fastq):
        """Teste básico de estatísticas FastQC"""
        stats = quality_control.run_fastqc(sample_fastq)
        assert stats['total_sequences'] > 0
        assert 0 <= stats['gc_content'] <= 100
    
    def test_quality_filtering(self, sample_fastq):
        """Teste de filtro de qualidade"""
        filtered = quality_control.filter_by_quality(
            sample_fastq, min_quality=20
        )
        assert filtered.count() <= stats['total_sequences']
```

### Testes Parametrizados
```python
@pytest.mark.parametrize("aligner,expected_format", [
    ("bwa", "SAM"),
    ("bowtie2", "SAM"),
    ("star", "BAM")
])
def test_alignment_output_format(aligner, expected_format):
    """Teste de formatos de saída dos alinhadores"""
    result = run_alignment(aligner, "sample.fastq")
    assert result.format == expected_format
```

### Mocking e Fixtures
```python
@pytest.fixture
def mock_genome_reference():
    """Mock de referência genômica para testes rápidos"""
    return MockGenome(chromosomes=["chr1", "chr2"], length=1000000)

@patch('src.alignment.bwa.run_bwa_mem')
def test_alignment_pipeline(mock_bwa, mock_genome_reference):
    """Teste de pipeline de alinhamento com mock"""
    mock_bwa.return_value = "aligned.bam"
    result = align_sequences("sample.fastq", mock_genome_reference)
    assert result == "aligned.bam"
    mock_bwa.assert_called_once()
```

## 📈 Integração Contínua

### Estratégia de Testes
1. **Pull Request**: Testes unitários + testes rápidos de integração
2. **Merge to main**: Suite completa de testes
3. **Release**: Testes de regressão + benchmarks
4. **Nightly**: Testes extensivos + validação de datasets grandes

### Notificação de Falhas
- **Slack/Teams**: Alertas imediatos para falhas críticas
- **Email**: Resumo diário de status dos testes
- **Dashboard**: Visualização de métricas de qualidade

## 🐛 Debugging e Troubleshooting

### Logs de Teste
```bash
# Executar com logs detalhados
pytest tests/ -v -s --log-cli-level=DEBUG

# Capturar saída de stdout/stderr
pytest tests/ --capture=no

# Parar no primeiro erro
pytest tests/ -x
```

### Teste Específico em Debug
```bash
# Executar teste específico com pdb
pytest tests/unit/test_variant_calling.py::test_gatk_haplotypecaller -vs --pdb

# Executar com profiling
pytest tests/integration/ --profile-svg
```

## 📚 Recursos Adicionais

### Documentação
- [Pytest Documentation](https://docs.pytest.org/)
- [Snakemake Testing](https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html)
- [Nextflow Testing Best Practices](https://www.nextflow.io/docs/latest/testing.html)

### Ferramentas Recomendadas
- **pytest**: Framework principal de testes
- **pytest-cov**: Cobertura de código
- **pytest-xdist**: Execução paralela
- **pytest-benchmark**: Testes de performance
- **hypothesis**: Testes baseados em propriedades
- **factory_boy**: Geração de dados de teste

## 🤝 Contribuições

Ao adicionar novos testes:
1. Seguir as convenções de nomenclatura
2. Incluir docstrings descritivas
3. Manter cobertura de código acima de 80%
4. Adicionar testes para casos edge
5. Documentar dependências específicas

---

**Nota**: Este framework de testes é essencial para manter a qualidade e confiabilidade do pipeline genômico. Todos os desenvolvedores devem executar os testes localmente antes de submeter pull requests.
