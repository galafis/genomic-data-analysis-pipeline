# 🔧 Troubleshooting do Variant Calling

Este README documenta diagnóstico, soluções para problemas comuns e tuning/performance para o módulo de Variant Calling, mantendo o mesmo padrão visual e organizacional do pipeline. Foco: respostas rápidas, checklist de causas prováveis, comandos prontos e integração com logs, scripts e configs.

## 📋 Propósito
- Guia de diagnóstico estruturado por sintomas e contexto
- Soluções para problemas comuns (memória, falhas de execução, baixa sensibilidade, workflow)
- Boas práticas de tuning e performance (CPU, RAM, IO, paralelismo)
- Integração com logs, scripts utilitários e arquivos de configuração

## 🗂️Âmbito e Tipos de Issues
- Recursos: memória insuficiente, saturação de CPU/IO, espaço em disco
- Falhas de execução: erros GATK/FreeBayes/bcftools, índices ausentes, referência incompatível
- Qualidade/resultado: baixa sensibilidade, alta taxa de falsos positivos, baixa concordância entre callers
- Workflow: dependências/containers, variáveis de ambiente, paths, orquestradores (Nextflow/Snakemake)

---

## 🧭 Fluxo Rápido de Diagnóstico
1) Identifique o SINTOMA principal (erro de memória, falha de etapa, qualidade baixa)
2) Reúna EVIDÊNCIAS: últimos 200–500 linhas de log, uso de recursos, versões
3) Valide PRÉ-REQUISITOS: referências, índices, checagens de integridade
4) Aplique AÇÕES corretivas mínimas e re-rodar em subset pequeno (intervalos/chr)
5) Escale a correção para todo o dataset (scatter-gather/parallel)

Comandos úteis:
- Ver recursos: ./troubleshooting/monitor_resources.sh
- Checar versões: ./utils/print_versions.sh
- Validar VCF: bcftools +fixploidy; bcftools norm; vt validate

---

## 💥 Problemas Comuns e Soluções

### 1) Erros de Memória (Java OOM, killed by OOM)
Sintomas: "OutOfMemoryError", processo encerrado com código 137, kernel OOM killer.
Soluções:
- Aumentar heap Java (GATK):
  export JAVA_OPTS="-Xmx32g -XX:+UseG1GC"
- Reduzir janela de processamento com intervalos:
  gatk HaplotypeCaller \
    -I sample.bam -R reference.fa -O out.vcf.gz \
    -L chr1:1-50,000,000 --native-pair-hmm-threads 4
- Ativar scatter-gather por cromossomo/interval list
- Certificar índices BAM/VCF/FAI presentes: samtools index, bcftools index, samtools faidx
- Validar referência e dicionário: gatk CreateSequenceDictionary

### 2) Falhas de Execução (parâmetros/arquivos)
Sintomas: "file not found", "Input files reference and BAM have different contigs", "Index not found".
Soluções:
- Conferir compatibilidade de contigs (chr vs no-chr):
  picard LiftoverVcf ou reindexar referência consistente
- Garantir índices:
  samtools index sample.bam; bcftools index variants.vcf.gz; samtools faidx reference.fa
- Atualizar paths de cache (VEP/SnpEff) em config/*.yaml
- Rodar utilitário de validação:
  python utils/check_inputs.py --bam sample.bam --ref reference.fa --vcf variants.vcf.gz

### 3) Baixa Sensibilidade / Alta Taxa de FP
Sintomas: menos variantes que o esperado ou muitas variantes espúrias.
Soluções:
- Normalizar e decompor variantes:
  bcftools norm -m -both -f reference.fa in.vcf.gz -Oz -o norm.vcf.gz
- Ajustar filtros de qualidade (caller-specific e pós-processamento)
- Usar regiões confiáveis (BED) e remover low-complexity
- Ajustar profundidade mínima (DP), qualidade de genótipo (GQ) e QUAL
- Benchmark rápido:
  python analysis/benchmark_analysis.py --truth data/truth_sets/NA12878.vcf.gz \
    --callers results/gatk.vcf.gz,results/freebayes.vcf.gz \
    --confident-regions data/truth_sets/confident_regions.bed \
    --output benchmark/quick_check/

### 4) Baixa Concordância entre Callers
Soluções:
- Harmonizar parâmetros de filtragem pós-chamada e normalização (bcftools norm)
- Comparar interseções/uniões com Venn/upset:
  Rscript analysis/comparison_plots.R --gatk A.vcf.gz --freebayes B.vcf.gz \
    --bcftools C.vcf.gz --sample SAMPLE --output plots/caller_comparison/
- Investigar regiões e tipos específicos (INDEL vs SNP, homopolímeros)

### 5) Performance Lenta
Soluções:
- Paralelizar por cromossomo/intervalos; usar --threads de forma equilibrada
- Habilitar IO eficiente (bgzip/tabix, discos SSD, local scratch)
- Containers otimizados com dependências nativas
- Analisar gargalos com ./troubleshooting/monitor_resources.sh e logs de tempo

### 6) Problemas em Workflows (Nextflow/Snakemake)
Soluções:
- Limpar cache e reprocessar steps falhos: nextflow -resume seletivo ou snakemake --rerun-incomplete
- Checar perfis/config (CPU/RAM por processo) e canais de entrada
- Validar ambientes/containers: singularity cache, docker pull atualizado

---

## 🧪 Exemplos de Solução (copy-paste)
- Reindexar e padronizar referência:
  samtools faidx data/reference/genome.fa
  gatk CreateSequenceDictionary -R data/reference/genome.fa -O data/reference/genome.dict
- Normalizar VCF e remover variantes inválidas:
  bcftools norm -f data/reference/genome.fa -m -both in.vcf.gz -Oz -o out.norm.vcf.gz
  bcftools view -e 'QUAL<30 || DP<10' out.norm.vcf.gz -Oz -o out.filtered.vcf.gz
- Rodar subset por intervalos para depuração:
  gatk HaplotypeCaller -I sample.bam -R reference.fa -L chr20 \
    -O debug_chr20.vcf.gz --native-pair-hmm-threads 4

---

## 🔗 Links Úteis
- Documentação do módulo: ../README.md
- Análises e relatórios: ../analysis/README.md
- Utils (validadores e métricas): ../utils/README.md
- Workflows: ../workflows/README.md
- GATK Best Practices: https://gatk.broadinstitute.org/
- bcftools/samtools: http://www.htslib.org/
- FreeBayes: https://github.com/freebayes/freebayes

---

## 🧾 Integração com Logs, Scripts e Config
- Logs: export LOG_LEVEL=DEBUG para verbosidade; checar logs em results/logs/ e .nextflow.log/.snakemake/log
- Scripts de suporte:
  - troubleshooting/monitor_resources.sh (CPU/RAM/IO)
  - utils/print_versions.sh (versões de ferramentas)
  - utils/run_benchmark.py (checagem comparativa rápida)
- Config: parâmetros centrais em config/*.yaml; use --config em scripts e workflows para reproducibilidade

---

## 📦 Requisitos e Ambiente
- Consistência de versões (R/Python/htslib, GATK/FreeBayes/bcftools)
- Índices presentes: .bai/.crai, .tbi/.csi, .fai e .dict
- Armazenamento: espaço para intermediários e relatórios (HTML/PDF/plots)

---

## ✅ Boas Práticas
- Rodadas de teste em cromossomos pequenos antes do lote completo
- Seeds fixos onde aplicável; versionamento de configs e containers
- Documentar parâmetros no relatório final e salvar JSON/TSV de métricas

Última atualização: $(date '+%Y-%m-%d')
Versão da seção: 1.0.0
Compatibilidade: Pipeline genomic-data-analysis-pipeline v2.0+
