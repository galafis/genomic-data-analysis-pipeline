# 🧪 Testes e Cobertura de Código

[![Test Coverage](https://github.com/galafis/genomic-data-analysis-pipeline/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/galafis/genomic-data-analysis-pipeline/actions/workflows/test-coverage.yml)
[![codecov](https://codecov.io/gh/galafis/genomic-data-analysis-pipeline/branch/master/graph/badge.svg)](https://codecov.io/gh/galafis/genomic-data-analysis-pipeline)

## 📋 Índice

- [Visão Geral](#visão-geral)
- [CI/CD Pipeline](#cicd-pipeline)
- [Executar Testes Localmente](#executar-testes-localmente)
- [Gerar Relatório de Cobertura](#gerar-relatório-de-cobertura)
- [Módulos Testados](#módulos-testados)
- [Estrutura de Testes](#estrutura-de-testes)

## 🔍 Visão Geral

Este projeto implementa testes automatizados com cobertura de código integrada ao pipeline CI/CD do GitHub Actions. Os testes garantem a qualidade e confiabilidade do código, especialmente dos módulos críticos de visualização e exportação de dados.

## 🚀 CI/CD Pipeline

O workflow de testes automatizados `.github/workflows/test-coverage.yml` é executado automaticamente em:

- **Push** para os branches `master`, `main` ou `develop`
- **Pull Requests** para os branches `master`, `main` ou `develop`
- **Execução manual** via GitHub Actions UI

### Funcionalidades do Pipeline

✅ Executa testes com `pytest` e gera relatórios de cobertura (HTML/XML/Terminal)  
✅ Upload automático de relatórios para Codecov  
✅ Comentários automáticos em PRs com métricas de cobertura  
✅ Geração de badge de cobertura no push para master  
✅ Artifacts HTML disponíveis para visualização detalhada (30 dias)

## 💻 Executar Testes Localmente

### Pré-requisitos

```bash
# Instalar dependências de teste
pip install pytest pytest-cov pytest-mock

# Instalar dependências do módulo
pip install pandas openpyxl
```

### Executar Todos os Testes

```bash
# Executar todos os testes no diretório tests/
pytest tests/ -v
```

### Executar Testes Específicos

```bash
# Executar apenas testes do módulo vcf_export_tools
pytest tests/test_vcf_export_tools.py -v

# Executar um teste específico
pytest tests/test_vcf_export_tools.py::test_export_to_csv -v
```

## 📈 Gerar Relatório de Cobertura

### Cobertura no Terminal

```bash
# Gerar relatório de cobertura no terminal
pytest tests/test_vcf_export_tools.py \
  --cov=src/visualization/interactive/vcf_export_tools \
  --cov-report=term \
  -v
```

### Relatório HTML Interativo

```bash
# Gerar relatório HTML detalhado
pytest tests/test_vcf_export_tools.py \
  --cov=src/visualization/interactive/vcf_export_tools \
  --cov-report=html \
  --cov-report=term \
  -v

# Abrir o relatório no navegador
# O arquivo está em: htmlcov/index.html
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
start htmlcov/index.html  # Windows
```

### Relatório XML (para CI/CD)

```bash
# Gerar relatório XML (usado pelo Codecov e outros serviços)
pytest tests/test_vcf_export_tools.py \
  --cov=src/visualization/interactive/vcf_export_tools \
  --cov-report=xml \
  --cov-report=term \
  -v
```

### Todos os Formatos Simultaneamente

```bash
# Gerar HTML, XML e exibir no terminal
pytest tests/test_vcf_export_tools.py \
  --cov=src/visualization/interactive/vcf_export_tools \
  --cov-report=html \
  --cov-report=xml \
  --cov-report=term \
  -v
```

## 📚 Módulos Testados

### `src/visualization/interactive/vcf_export_tools.py`

Módulo para exportação de dados VCF em múltiplos formatos.

**Funcionalidades Testadas:**
- ✅ Exportação para CSV
- ✅ Exportação para Excel (.xlsx)
- ✅ Exportação para JSON
- ✅ Tratamento de erros e validação de entrada
- ✅ Formatação de dados e tipos de coluna

**Cobertura Atual:**
Consulte o [badge de cobertura](https://codecov.io/gh/galafis/genomic-data-analysis-pipeline) para ver a cobertura atualizada.

## 📁 Estrutura de Testes

```
tests/
├── test_vcf_export_tools.py    # Testes do módulo vcf_export_tools
├── unit/                      # Testes unitários (futuro)
├── integration/               # Testes de integração (futuro)
└── README.md                  # Documentação dos testes
```

### Exemplo de Teste

```python
import pytest
from src.visualization.interactive.vcf_export_tools import export_to_csv

def test_export_to_csv(tmp_path):
    """Testa exportação de dados para CSV."""
    # Preparar dados de teste
    data = [
        {'chr': '1', 'pos': 1000, 'ref': 'A', 'alt': 'T'},
        {'chr': '2', 'pos': 2000, 'ref': 'G', 'alt': 'C'}
    ]
    
    # Exportar para CSV
    output_file = tmp_path / "output.csv"
    export_to_csv(data, str(output_file))
    
    # Verificar que o arquivo foi criado
    assert output_file.exists()
    
    # Verificar conteúdo
    content = output_file.read_text()
    assert 'chr,pos,ref,alt' in content
    assert '1,1000,A,T' in content
```

## 📊 Visualizar Cobertura no GitHub

### Via GitHub Actions

1. Acesse a aba **[Actions](https://github.com/galafis/genomic-data-analysis-pipeline/actions)**
2. Selecione o workflow **"Test Coverage"**
3. Clique em uma execução recente
4. Baixe o artifact **"coverage-report-html"** para visualizar o relatório completo

### Via Codecov

Visite o [dashboard do Codecov](https://codecov.io/gh/galafis/genomic-data-analysis-pipeline) para:
- Ver gráficos de tendência de cobertura
- Comparar cobertura entre commits
- Identificar linhas não cobertas por testes
- Gerar relatórios detalhados por arquivo/função

## 🛠️ Próximos Passos

- [ ] Adicionar testes para módulos de alinhamento
- [ ] Implementar testes de integração
- [ ] Aumentar cobertura para >90%
- [ ] Adicionar testes de performance/benchmarking
- [ ] Configurar testes para diferentes plataformas (Windows/macOS/Linux)

## 👤 Autor

**Gabriel** - Implementação do sistema de testes e CI/CD

---

📌 **Nota**: Para contribuir com novos testes, consulte o guia em `tests/README.md`
