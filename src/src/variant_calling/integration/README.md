# 🔄 Integração Downstream — Variant Calling

Este diretório reúne utilitários e guias para integrar os resultados de chamada de variantes com etapas downstream do pipeline: anotação funcional, visualização e exportação para formatos interoperáveis. Mantém o mesmo padrão visual e organizacional dos demais diretórios deste módulo.

## 📌 Propósito
- Anotação: enriquecer VCF/BCF com efeitos funcionais (impacto, genes, consequências) usando VEP/SnpEff e regras do projeto.
- Visualização: preparar artefatos e tracks para IGV/IGV.js, Shiny dashboards e relatórios HTML/PDF.
- Exportação: converter VCFs em formatos tabulares e intercambiáveis (TSV/CSV, JSON, BED, gVCF, MAF, XML), preservando metadados essenciais.

## 🗂️ Estrutura
- to_annotation.sh — orquestra anotações com VEP/SnpEff e pós-processamento (bcftools, vt).
- to_visualization.sh — gera tracks (BED, bedGraph), subsets por cromossomo/genes e índices para browsers.
- export_formats.py — conversões e normalizações de VCF para múltiplos formatos com validações.
- README.md — este arquivo.

## 🚀 Fluxos de Execução Típicos
1) Anotação funcional
```
./integration/to_annotation.sh \
  --input results/variant_calling/sample.vcf.gz \
  --reference data/reference/genome.fa \
  --annotation-config src/src/variant_calling/config/annotation_config.yaml \
  --out annotated/sample.annotated.vcf.gz
```
Saída: VCF anotado + índices (.tbi) e logs. Parâmetros do VEP/SnpEff são herdados do arquivo de configuração do módulo.

2) Preparação para visualização
```
./integration/to_visualization.sh \
  --input annotated/sample.annotated.vcf.gz \
  --tracks-out viz/tracks/ \
  --regions data/reference/regions/genes_of_interest.bed
```
Saídas: BED/BEDPE, listas de variantes prioritárias, subsets por região, índices para IGV.

3) Exportação para formatos
```
python integration/export_formats.py \
  --input annotated/sample.annotated.vcf.gz \
  --output-dir exports/ \
  --output-format tsv,json,maf,bed
```
Saídas: arquivos normalizados, com colunas canônicas e validações básicas (POS, REF/ALT, INFO-chave).

## 🧩 Exemplos de Scripts
- to_annotation.sh
  - Normaliza (bcftools norm), anota (VEP/SnpEff), reindexa (tabix) e registra versões.
  - Flags úteis: `--vep-cache`, `--snpeff-db`, `--threads`, `--min-qual`.
- to_visualization.sh
  - Gera BED por impacto/genes, prepara session files para IGV, exporta TSV leve para dashboards.
  - Flags úteis: `--regions`, `--split-by-chrom`, `--max-variants`.
- export_formats.py
  - Converte para TSV/CSV/JSON/MAF/BED/XML; aplica schema mínimo e verifica duplicatas.
  - Flags úteis: `--fields`, `--compress`, `--schema`, `--drop-missing`.

## ✅ Recomendações de Uso
- Sempre normalize VCF antes de qualquer etapa (left-align, split multialélicos).
- Congele versões de referência e bancos (VEP cache, SnpEff DB) no annotation_config.yaml.
- Gere índices (.tbi) após cada transformação para acesso aleatório eficiente.
- Para grandes coortes, pagine por cromossomo/região e use `--threads`.

## 🔁 Interoperabilidade e Boas Práticas
- Preserve metadados VCF nas exportações (INFO-chave: AC, AF, DP, ANN/CSQ).
- Use nomenclatura consistente de amostras e cohort IDs em todos os formatos.
- Valide MAF/TSV com um schema declarativo (p. ex., JSON Schema) antes de ingestão.
- Para integração com R/Bioconductor, exporte SummarizedExperiment-parity TSV/CSV (wide/long) conforme docs/.
- Mantenha hash/manifest das entradas-saídas para reprodutibilidade.

## 🔗 Links Relacionados
- README do módulo variant_calling: src/src/variant_calling/README.md
- Configurações de anotação: src/src/variant_calling/config/annotation_config.yaml
- Workflows: src/src/variant_calling/workflows/ (Nextflow/Snakemake/CWL)
- Módulo de anotação: src/annotation/
- Módulo de visualização: src/visualization/
- README do pipeline principal: README.md na raiz do repositório

## 🧪 Testes e Validação
- Compare contagens e checksums entre VCF de entrada e exportações derivadas.
- Execute validadores (vcf-validator, maf-tools) nas saídas.
- Amostras de teste em: src/src/variant_calling/test_data/

## 🗒️ Notas de Compatibilidade
- Requer: bcftools ≥1.15, htslib/tabix, VEP ≥105 ou SnpEff ≥5.0.
- Compatível com o perfil de containers definidos em containers/ do projeto.

## 📄 Licença
Este diretório segue a licença MIT do projeto. Veja LICENSE na raiz.
