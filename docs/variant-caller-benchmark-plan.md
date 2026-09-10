# Variant caller benchmark plan / Plano de comparação de chamadas de variantes

## English

This is an unexecuted comparison protocol. The repository does not contain measured caller-ranking results. The earlier illustrative result table has been removed to avoid presenting invented measurements as evidence.

1. Select a reference genome, truth VCF and confident-region BED with recorded versions and checksums.
2. Use the same aligned input, reference build, sample identity and analysis regions for each caller.
3. Record caller version, command, threads, machine configuration and exit status.
4. Normalize representations against the same reference before comparing variants. Treat SNPs and indels separately.
5. Measure precision and recall within the declared confident regions; record wall time and peak memory separately.
6. Preserve commands, logs, input hashes and output hashes with the result table. Label unsuccessful runs explicitly.

No result should be published until this protocol is executed and independently checked. The VCF exporter exercised by this repository is not a variant caller or a truth-set comparison tool.

## Português

Este é um protocolo de comparação ainda não executado. O repositório não contém resultados medidos que classifiquem os programas. A tabela ilustrativa anterior foi removida para evitar apresentar medições inventadas como evidência.

1. Escolha genoma de referência, VCF de verdade e BED de regiões confiáveis, registrando versões e hashes.
2. Use a mesma entrada alinhada, referência, amostra e regiões de análise em cada programa.
3. Registre versão, comando, threads, configuração da máquina e código de saída.
4. Normalize as representações com a mesma referência antes da comparação. Separe SNPs e indels.
5. Meça precisão e sensibilidade nas regiões declaradas; registre tempo e pico de memória separadamente.
6. Preserve comandos, logs e hashes das entradas e saídas com a tabela de resultados. Identifique execuções malsucedidas.

Não publique resultados antes de executar e conferir o protocolo. O exportador VCF exercitado neste repositório não realiza chamada de variantes nem comparação com conjunto de verdade.

[Reviewed executable path / Caminho executável revisado](../README.md)
