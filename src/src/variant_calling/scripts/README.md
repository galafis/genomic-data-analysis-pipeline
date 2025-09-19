# Variant Calling Scripts

Este diretório contém scripts para chamada e processamento de variantes genômicas utilizando diferentes ferramentas e algoritmos.

## Scripts Disponíveis

### `run_gatk_haplotypecaller.sh`
**Descrição:** Script para chamada de variantes utilizando GATK HaplotypeCaller, baseado em algoritmo de assembly local para detecção precisa de SNPs e indels.

**Uso básico:**
```bash
./run_gatk_haplotypecaller.sh -i input.bam -r reference.fasta -o output.vcf
```

**Parâmetros principais:**
- `-i, --input`: Arquivo BAM de entrada (obrigatório)
- `-r, --reference`: Genoma de referência FASTA (obrigatório)
- `-o, --output`: Arquivo VCF de saída (obrigatório)
- `-t, --threads`: Número de threads (padrão: 4)
- `--intervals`: Regiões específicas para análise
- `--emit-ref-confidence`: Modo GVCF para análise conjunta

### `run_freebayes.sh`
**Descrição:** Script para chamada de variantes usando FreeBayes, caller bayesiano para detecção de variantes complexas e poliploidias.

**Uso básico:**
```bash
./run_freebayes.sh -i input.bam -r reference.fasta -o output.vcf
```

**Parâmetros principais:**
- `-i, --input`: Arquivo BAM de entrada (obrigatório)
- `-r, --reference`: Genoma de referência FASTA (obrigatório)
- `-o, --output`: Arquivo VCF de saída (obrigatório)
- `-p, --ploidy`: Ploidia do organismo (padrão: 2)
- `-q, --min-base-quality`: Qualidade mínima da base (padrão: 20)
- `-m, --min-mapping-quality`: Qualidade mínima do mapeamento (padrão: 30)
- `--targets`: Arquivo BED com regiões alvo

### `run_bcftools.sh`
**Descrição:** Script para chamada de variantes utilizando BCFtools mpileup/call, adequado para análises rápidas e datasets grandes.

**Uso básico:**
```bash
./run_bcftools.sh -i input.bam -r reference.fasta -o output.vcf
```

**Parâmetros principais:**
- `-i, --input`: Arquivo BAM de entrada (obrigatório)
- `-r, --reference`: Genoma de referência FASTA (obrigatório)
- `-o, --output`: Arquivo VCF de saída (obrigatório)
- `-q, --min-BQ`: Qualidade mínima da base (padrão: 20)
- `-Q, --min-MQ`: Qualidade mínima do mapeamento (padrão: 30)
- `--multiallelic-caller`: Usar caller multialélico
- `--variants-only`: Reportar apenas posições variantes

### `run_variant_filter.sh`
**Descrição:** Script para filtragem e controle de qualidade de variantes chamadas, aplicando filtros customizáveis por qualidade, profundidade e frequência alélica.

**Uso básico:**
```bash
./run_variant_filter.sh -i input.vcf -o filtered.vcf
```

**Parâmetros principais:**
- `-i, --input`: Arquivo VCF de entrada (obrigatório)
- `-o, --output`: Arquivo VCF filtrado de saída (obrigatório)
- `--min-qual`: Qualidade mínima da variante (padrão: 30)
- `--min-depth`: Profundidade mínima (padrão: 10)
- `--max-depth`: Profundidade máxima (padrão: 1000)
- `--min-af`: Frequência alélica mínima (padrão: 0.1)
- `--hard-filters`: Aplicar filtros GATK hard-filtering
- `--vqsr`: Usar GATK VQSR (Variant Quality Score Recalibration)

### `run_variant_annotation.sh`
**Descrição:** Script para anotação funcional de variantes utilizando ferramentas como SnpEff/VEP, adicionando informações sobre efeitos biológicos e consequências funcionais.

**Uso básico:**
```bash
./run_variant_annotation.sh -i input.vcf -o annotated.vcf -d database
```

**Parâmetros principais:**
- `-i, --input`: Arquivo VCF de entrada (obrigatório)
- `-o, --output`: Arquivo VCF anotado de saída (obrigatório)
- `-d, --database`: Base de dados de anotação (obrigatório)
- `--tool`: Ferramenta de anotação (snpeff|vep, padrão: snpeff)
- `--genome`: Versão do genoma (hg38|hg19|mm10, etc.)
- `--fields`: Campos específicos para anotação
- `--clinical`: Incluir anotações clínicas (ClinVar, COSMIC)

## Integração em Workflows

### Pipeline Sequencial
```bash
# 1. Chamada de variantes
./run_gatk_haplotypecaller.sh -i sample.bam -r ref.fa -o raw_variants.vcf

# 2. Filtragem
./run_variant_filter.sh -i raw_variants.vcf -o filtered_variants.vcf

# 3. Anotação
./run_variant_annotation.sh -i filtered_variants.vcf -o final_variants.vcf -d hg38
```

### Integração com Workflow Managers

#### Snakemake
```python
rule variant_calling:
    input:
        bam="{sample}.bam",
        ref="reference.fasta"
    output:
        vcf="{sample}_variants.vcf"
    shell:
        "./run_gatk_haplotypecaller.sh -i {input.bam} -r {input.ref} -o {output.vcf}"
```

#### Nextflow
```groovy
process VARIANT_CALLING {
    input:
    tuple val(sample_id), path(bam)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}_variants.vcf")
    
    script:
    """
    ./run_gatk_haplotypecaller.sh -i ${bam} -r ${reference} -o ${sample_id}_variants.vcf
    """
}
```

### Processamento em Lote
```bash
# Script para processar múltiplas amostras
for sample in samples/*.bam; do
    base=$(basename ${sample} .bam)
    ./run_gatk_haplotypecaller.sh -i ${sample} -r reference.fasta -o ${base}_variants.vcf
    ./run_variant_filter.sh -i ${base}_variants.vcf -o ${base}_filtered.vcf
    ./run_variant_annotation.sh -i ${base}_filtered.vcf -o ${base}_annotated.vcf -d hg38
done
```

## Configuração e Dependências

### Ferramentas Necessárias
- GATK (>= 4.0)
- FreeBayes (>= 1.3)
- BCFtools (>= 1.9)
- SnpEff ou VEP
- SAMtools (>= 1.9)

### Variáveis de Ambiente
```bash
export GATK_PATH=/path/to/gatk
export SNPEFF_PATH=/path/to/snpeff
export REFERENCE_DIR=/path/to/references
```

## Logs e Monitoramento

Todos os scripts geram logs detalhados em `logs/` com informações sobre:
- Estatísticas de processamento
- Tempo de execução
- Uso de recursos
- Métricas de qualidade

## Integração com o Pipeline Principal

Estes scripts são integrados automaticamente no pipeline principal através do:
- **Config**: `config/variant_calling.yaml`
- **Workflow**: `workflows/variant_calling.snakefile`
- **Documentação**: Consulte o [README principal](../../README.md) para configuração completa

## Suporte e Troubleshooting

Para problemas específicos:
1. Verifique os logs em `logs/`
2. Consulte a documentação das ferramentas
3. Execute com `--dry-run` para validar parâmetros
4. Use `--help` em qualquer script para detalhes

## Referências

- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- [FreeBayes Manual](https://github.com/freebayes/freebayes)
- [BCFtools Documentation](http://samtools.github.io/bcftools/)
- [SnpEff Manual](http://pcingola.github.io/SnpEff/)
