# Validation record / Registro de validação

Review date / Data da revisão: 2026-09-10.

## Change and reason / Mudança e motivo

EN: Handled missing depth and empty datasets; rejected malformed records and invalid export names; corrected two test expectations to match actual fixture depth and pandas integer types.

PT: Tratados profundidade ausente e conjuntos vazios; rejeitados registros malformados e nomes de exportação inválidos; corrigidas duas expectativas de testes para os dados da amostra e inteiros pandas.

## Reproduce / Reproduzir

```sh
python -m venv .venv
# Activate .venv for your shell / Ative .venv no seu terminal
python -m pip install -r requirements-validation.txt
python -m pytest -q
python -m examples.review_demo
```

## Example evidence / Evidência do exemplo

At QUAL ≥ 30 and DP ≥ 10, demo-1 and demo-4 remain. / Com QUAL ≥ 30 e DP ≥ 10, permanecem demo-1 e demo-4.

EN: The suite includes normal operations and regression cases for the corrected behavior. The example checks values produced by the implementation. Use the linked workflow to inspect the result for a specific commit; no performance benchmark is inferred from a passing build.

PT: A suíte inclui operações válidas e regressões dos comportamentos corrigidos. O exemplo verifica valores produzidos pela implementação. Consulte a automação para conferir o resultado de um commit específico; aprovação de compilação não implica benchmark de desempenho.

## Limits / Limites

EN: The reviewed path is VCF parsing, filtering and export. Other directories contain external-tool wrappers and workflow experiments; full sequencing, alignment and variant-calling pipelines were not executed in this review. This parser does not interpret genotype FORMAT fields, normalize variants or assign clinical significance.

PT: O caminho revisado é leitura, filtragem e exportação de VCF. Outros diretórios contêm integrações com ferramentas externas e experimentos de workflow; sequenciamento, alinhamento e chamada de variantes completos não foram executados nesta revisão. O leitor não interpreta FORMAT de genótipos, normaliza variantes ou atribui significado clínico.

[Return to README / Voltar ao README](../README.md)

## Verified suite / Suíte verificada

**17 software tests passed / testes de software aprovados.**

README Mermaid syntax and local documentation links were checked. / A sintaxe Mermaid do README e os links locais da documentação foram conferidos.
