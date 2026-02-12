# 📊 Genomic Data Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

[English](#english) | [Português](#português)

---

## English

### 🎯 Overview

**Genomic Data Analysis Pipeline** — Advanced data science project: genomic-data-analysis-pipeline

Total source lines: **9,209** across **46** files in **3** languages.

### ✨ Key Features

- **Production-Ready Architecture**: Modular, well-documented, and following best practices
- **Comprehensive Implementation**: Complete solution with all core functionality
- **Clean Code**: Type-safe, well-tested, and maintainable codebase
- **Easy Deployment**: Docker support for quick setup and deployment

### 🚀 Quick Start

#### Prerequisites
- Python 3.12+


#### Installation

1. **Clone the repository**
```bash
git clone https://github.com/galafis/genomic-data-analysis-pipeline.git
cd genomic-data-analysis-pipeline
```

2. **Create virtual environment**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```





### 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov --cov-report=html

# Run with verbose output
pytest -v
```

### 📁 Project Structure

```
genomic-data-analysis-pipeline/
├── assets/
├── config/
│   └── README.md
├── containers/
│   └── README.md
├── data/
│   └── README.md
├── docs/
│   └── README.md
├── results/
│   └── README.md
├── scripts/
│   └── README.md
├── src/
│   ├── alignment/
│   │   ├── benchmarks/
│   │   ├── protocols/
│   │   ├── tools/
│   │   ├── README.md
│   │   └── bwa_mem2_align.py
│   ├── annotation/
│   │   └── README.md
│   ├── preprocessing/
│   │   ├── filtering/
│   │   ├── normalization/
│   │   ├── quality_control/
│   │   ├── trimming/
│   │   ├── README.md
│   │   └── quality_control.py
│   ├── scripts/
│   │   └── README.md
│   ├── src/
│   │   └── variant_calling/
│   ├── variant_calling/
│   │   ├── README.md
│   │   ├── README_batch.md
│   │   ├── batch_troubleshooting.md
│   │   ├── batch_variant_calling.py
│   │   ├── qc_variant_stats.py
│   │   ├── variant_annotation.py
│   │   ├── variant_batch_summary.py
│   │   ├── variant_caller.py
│   │   ├── variant_calling_protocol.md
│   │   ├── variant_filter.py
│   │   ├── vcf_report_generator.py
│   │   └── vcf_visualization.py
│   ├── visualization/
│   │   ├── interactive/
│   │   ├── README.md
│   │   └── plot_variants.py
│   ├── workflows/
│   │   ├── nextflow/
│   │   └── README.md
│   └── README.md
├── tests/
│   ├── README.md
│   └── test_vcf_export_tools.py
├── workflows/
│   ├── nextflow/
│   └── README.md
├── README.md
├── TESTING.md
└── environment.yml
```

### 🛠️ Tech Stack

| Technology | Usage |
|------------|-------|
| Python | 35 files |
| Shell | 10 files |
| R | 1 files |

### 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### 👤 Author

**Gabriel Demetrios Lafis**

- GitHub: [@galafis](https://github.com/galafis)
- LinkedIn: [Gabriel Demetrios Lafis](https://linkedin.com/in/gabriel-demetrios-lafis)

---

## Português

### 🎯 Visão Geral

**Genomic Data Analysis Pipeline** — Advanced data science project: genomic-data-analysis-pipeline

Total de linhas de código: **9,209** em **46** arquivos em **3** linguagens.

### ✨ Funcionalidades Principais

- **Arquitetura Pronta para Produção**: Modular, bem documentada e seguindo boas práticas
- **Implementação Completa**: Solução completa com todas as funcionalidades principais
- **Código Limpo**: Type-safe, bem testado e manutenível
- **Fácil Implantação**: Suporte Docker para configuração e implantação rápidas

### 🚀 Início Rápido

#### Pré-requisitos
- Python 3.12+


#### Instalação

1. **Clone the repository**
```bash
git clone https://github.com/galafis/genomic-data-analysis-pipeline.git
cd genomic-data-analysis-pipeline
```

2. **Create virtual environment**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```




### 🧪 Testes

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov --cov-report=html

# Run with verbose output
pytest -v
```

### 📁 Estrutura do Projeto

```
genomic-data-analysis-pipeline/
├── assets/
├── config/
│   └── README.md
├── containers/
│   └── README.md
├── data/
│   └── README.md
├── docs/
│   └── README.md
├── results/
│   └── README.md
├── scripts/
│   └── README.md
├── src/
│   ├── alignment/
│   │   ├── benchmarks/
│   │   ├── protocols/
│   │   ├── tools/
│   │   ├── README.md
│   │   └── bwa_mem2_align.py
│   ├── annotation/
│   │   └── README.md
│   ├── preprocessing/
│   │   ├── filtering/
│   │   ├── normalization/
│   │   ├── quality_control/
│   │   ├── trimming/
│   │   ├── README.md
│   │   └── quality_control.py
│   ├── scripts/
│   │   └── README.md
│   ├── src/
│   │   └── variant_calling/
│   ├── variant_calling/
│   │   ├── README.md
│   │   ├── README_batch.md
│   │   ├── batch_troubleshooting.md
│   │   ├── batch_variant_calling.py
│   │   ├── qc_variant_stats.py
│   │   ├── variant_annotation.py
│   │   ├── variant_batch_summary.py
│   │   ├── variant_caller.py
│   │   ├── variant_calling_protocol.md
│   │   ├── variant_filter.py
│   │   ├── vcf_report_generator.py
│   │   └── vcf_visualization.py
│   ├── visualization/
│   │   ├── interactive/
│   │   ├── README.md
│   │   └── plot_variants.py
│   ├── workflows/
│   │   ├── nextflow/
│   │   └── README.md
│   └── README.md
├── tests/
│   ├── README.md
│   └── test_vcf_export_tools.py
├── workflows/
│   ├── nextflow/
│   └── README.md
├── README.md
├── TESTING.md
└── environment.yml
```

### 🛠️ Stack Tecnológica

| Tecnologia | Uso |
|------------|-----|
| Python | 35 files |
| Shell | 10 files |
| R | 1 files |

### 📄 Licença

Este projeto está licenciado sob a Licença MIT - veja o arquivo [LICENSE](LICENSE) para detalhes.

### 👤 Autor

**Gabriel Demetrios Lafis**

- GitHub: [@galafis](https://github.com/galafis)
- LinkedIn: [Gabriel Demetrios Lafis](https://linkedin.com/in/gabriel-demetrios-lafis)
