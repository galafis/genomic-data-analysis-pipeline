# Guia de Troubleshooting - Pipeline Batch Multicase

## Índice

1. [Checklist Pré-Execução](#checklist-pré-execução)
2. [Erros Comuns e Soluções](#erros-comuns-e-soluções)
3. [Problemas com Índices](#problemas-com-índices)
4. [Problemas com VCF](#problemas-com-vcf)
5. [Problemas com Recursos](#problemas-com-recursos)
6. [Problemas com Anotações](#problemas-com-anotações)
7. [VCFs Vazios](#vcfs-vazios)
8. [Exemplos de Comandos](#exemplos-de-comandos)
9. [Análise de Logs e Debugging](#análise-de-logs-e-debugging)
10. [Contribuição](#contribuição)

## Checklist Pré-Execução

### Verificações Essenciais
- [ ] **Arquivos de entrada**: BAM/CRAM indexados e acessíveis
- [ ] **Genoma de referência**: FASTA indexado (.fai) disponível
- [ ] **Recursos**: VCF de variantes conhecidas indexados
- [ ] **Espaço em disco**: Suficiente para outputs intermediários e finais
- [ ] **Memória**: Adequada para número de amostras simultâneas
- [ ] **Permissões**: Leitura nos inputs e escrita no diretório de output
- [ ] **Configuração**: Parâmetros de qualidade e filtros definidos

### Estrutura de Diretórios
```
workdir/
├── inputs/
│   ├── bams/
│   └── references/
├── outputs/
│   ├── variant_calling/
│   ├── filtering/
│   ├── annotation/
│   └── qc/
└── logs/
```

## Erros Comuns e Soluções

### 1. Erro de Memória Insuficiente
**Sintoma**: `java.lang.OutOfMemoryError`

**Soluções**:
```bash
# Aumentar heap do Java
export JAVA_OPTS="-Xmx32g -Xms8g"

# Reduzir número de threads simultâneas
--threads 4  # ao invés de 8 ou 16

# Processar em lotes menores
--batch-size 10  # ao invés de 50
```

### 2. Arquivo Não Encontrado
**Sintoma**: `FileNotFoundException` ou `No such file or directory`

**Verificações**:
```bash
# Verificar existência e permissões
ls -la /path/to/file
file /path/to/file

# Verificar índices
ls -la reference.fasta*
ls -la input.bam*
```

### 3. Erro de Formato de Arquivo
**Sintoma**: `Malformed VCF/BAM header`

**Validação**:
```bash
# Validar BAM
samtools quickcheck input.bam

# Validar VCF
vcf-validator input.vcf.gz

# Recriar índice se necessário
samtools index input.bam
tabix -p vcf input.vcf.gz
```

## Problemas com Índices

### Índices Ausentes ou Corrompidos
```bash
# Recriar índice FASTA
samtools faidx reference.fasta

# Recriar índice BAM
samtools index input.bam

# Recriar índice VCF
tabix -p vcf variants.vcf.gz

# Verificar integridade
md5sum reference.fasta.fai
```

### Incompatibilidade de Versões
```bash
# Verificar versão do formato
samtools view -H input.bam | grep @HD

# Atualizar se necessário
samtools view -h input.bam | samtools view -bS - > updated.bam
samtools index updated.bam
```

## Problemas com VCF

### Headers Inconsistentes
```bash
# Corrigir header VCF
bcftools view -h input.vcf.gz > header.txt
# Editar header.txt manualmente
bcftools reheader -h header.txt input.vcf.gz > corrected.vcf.gz
```

### Coordenadas Inválidas
```bash
# Validar coordenadas
bcftools view input.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'

# Filtrar coordenadas válidas
bcftools view -i 'POS>0 && INFO/END<=250000000' input.vcf.gz
```

### Encoding de Caracteres
```bash
# Verificar encoding
file -i input.vcf

# Converter se necessário
iconv -f ISO-8859-1 -t UTF-8 input.vcf > input_utf8.vcf
```

## Problemas com Recursos

### dbSNP/gnomAD Indisponíveis
```bash
# Verificar conectividade
wget --spider ftp://ftp.ncbi.nih.gov/snp/latest_release/

# Download local se necessário
wget ftp://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz

# Configurar caminho local
--dbsnp /local/path/dbsnp.vcf.gz
```

### Versões Incompatíveis do Genoma
```bash
# Verificar build do genoma
grep '^##reference' input.vcf.gz

# Converter coordenadas se necessário
crossmapvcf.py hg38ToHg19.over.chain input_hg38.vcf hg19.fa output_hg19.vcf
```

## Problemas com Anotações

### VEP/SnpEff Falhas
```bash
# Verificar cache VEP
vep --list

# Download cache se necessário
vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/cache

# Teste simples
echo "1\t100000\t.\tA\tG\t.\t.\t." | vep --format vcf
```

### Anotações Customizadas
```bash
# Verificar formato BED/GFF
head -n 5 annotations.bed
bedtools sort -i annotations.bed > annotations_sorted.bed

# Validar overlap
bedtools intersect -a variants.vcf -b annotations.bed -wa -wb
```

## VCFs Vazios

### Diagnóstico
```bash
# Contar variantes
bcftools view -H variants.vcf.gz | wc -l

# Verificar filtros aplicados
bcftools view -f PASS variants.vcf.gz | head

# Estatísticas gerais
bcftools stats variants.vcf.gz
```

### Possíveis Causas
1. **Filtros muito restritivos**
2. **Região de interesse muito pequena**
3. **Baixa cobertura nas amostras**
4. **Parâmetros de calling muito conservadores**

### Soluções
```bash
# Relaxar filtros temporariamente
bcftools view -i 'QUAL>10' variants.vcf.gz

# Verificar cobertura
samtools depth input.bam | awk '{sum+=$3; count++} END {print "Average coverage:", sum/count}'

# Ajustar parâmetros de calling
--min-base-quality 15  # ao invés de 20
--min-mapping-quality 10  # ao invés de 20
```

## Exemplos de Comandos

### Variant Calling Básico
```bash
# GATK HaplotypeCaller
gatk HaplotypeCaller \
  -R reference.fasta \
  -I input.bam \
  -O raw_variants.vcf.gz \
  --native-pair-hmm-threads 4
```

### Filtragem
```bash
# Filtros básicos de qualidade
bcftools view -i 'QUAL>20 && DP>10 && MQ>40' \
  raw_variants.vcf.gz -O z -o filtered_variants.vcf.gz

# Filtros específicos para multicase
bcftools view -i 'AC>1 && AN>10' filtered_variants.vcf.gz
```

### Anotação
```bash
# VEP annotation
vep --input_file variants.vcf.gz \
    --output_file annotated.vcf.gz \
    --format vcf --vcf --compress_output bgzip \
    --species homo_sapiens --assembly GRCh38 \
    --cache --offline
```

### Quality Control
```bash
# bcftools stats
bcftools stats variants.vcf.gz > variants.stats

# Plot stats
plot-vcfstats -p plots/ variants.stats

# Verificar Ti/Tv ratio
bcftools query -f '%TYPE\n' variants.vcf.gz | sort | uniq -c
```

## Análise de Logs e Debugging

### Localização de Logs
```bash
# Logs do sistema
tail -f /var/log/messages
tail -f /var/log/syslog

# Logs específicos do pipeline
tail -f logs/variant_calling.log
tail -f logs/filtering.log
tail -f logs/annotation.log
```

### Informações Importantes nos Logs
```bash
# Erros críticos
grep -i "error\|exception\|failed" logs/*.log

# Warnings importantes
grep -i "warning\|warn" logs/*.log

# Progresso e timing
grep -i "progress\|completed\|elapsed" logs/*.log

# Uso de recursos
grep -i "memory\|cpu\|disk" logs/*.log
```

### Debugging Avançado
```bash
# Ativar modo verbose
export JAVA_OPTS="$JAVA_OPTS -Dgatk.logging.level=DEBUG"

# Profiling de memória
jmap -dump:format=b,file=heap.hprof <PID>

# Monitoramento em tempo real
htop
iotop
df -h
```

### Logs Customizados
```bash
# Adicionar timestamps
exec > >(while read line; do echo "[$(date '+%Y-%m-%d %H:%M:%S')] $line"; done) 2>&1

# Log rotativo
logrotate -f /etc/logrotate.d/pipeline
```

## Links Internos do Pipeline

- [Configuração Principal](../config/)
- [Scripts de Variant Calling](./variant_calling_scripts.py)
- [Módulo de Filtragem](./filtering/)
- [Anotação Automática](./annotation/)
- [Quality Control](./qc/)
- [Utilitários](../utils/)
- [Documentação da API](../docs/api.md)
- [Exemplos de Configuração](../examples/)

## Contribuição

### Como Contribuir

1. **Reporte Issues**: Documente problemas encontrados
2. **Sugira Melhorias**: Propose soluções para problemas recorrentes
3. **Atualize Documentação**: Mantenha este guia atualizado
4. **Teste Casos Específicos**: Valide soluções em diferentes cenários

### Template para Novos Problemas

```markdown
### Novo Problema: [Título Descritivo]

**Sintoma**: Descrição do erro observado

**Contexto**:
- Versão do pipeline: X.Y.Z
- Sistema operacional: Ubuntu 20.04
- Recursos: 32GB RAM, 8 CPU cores
- Número de amostras: N

**Reprodução**:
1. Comando executado
2. Arquivos de entrada utilizados
3. Parâmetros específicos

**Solução Proposta**:
Descrição da solução encontrada

**Comandos de Teste**:
```bash
# Comandos para validar a solução
```

**Links Relacionados**:
- Issues relacionadas
- Documentação relevante
```

### Maintainers

- [@maintainer1](https://github.com/maintainer1) - Variant Calling
- [@maintainer2](https://github.com/maintainer2) - Quality Control
- [@maintainer3](https://github.com/maintainer3) - Annotation

### Histórico de Atualizações

- **v1.0** (2024-01): Versão inicial do guia
- **v1.1** (2024-03): Adicionado seção de VCFs vazios
- **v1.2** (2024-06): Melhorias na seção de debugging

---

> 💡 **Dica**: Mantenha sempre backups dos dados originais antes de aplicar correções!
> 
> 🔧 **Support**: Para problemas não cobertos neste guia, abra uma issue no repositório.
>
> 📚 **Documentação Adicional**: Consulte a [Wiki do projeto](../../wiki) para informações detalhadas.
