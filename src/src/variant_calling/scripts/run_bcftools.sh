#!/bin/bash
# SCRIPT: run_bcftools.sh
# AUTHOR: Genomic Data Analysis Pipeline Team
# VERSION: 1.0
# DATE: 2025-09-19
# DESCRIPTION: Pipeline robusto para chamada de variantes com bcftools.
# USO:
#   ./run_bcftools.sh -b input.bam -r reference.fa -o output.vcf.gz [ -t threads ] [ -c config.txt ]

set -euo pipefail

# Defaults
THREADS=8
CONFIG_FILE=""
LOG_DIR="./logs"
LOG_FILE="$LOG_DIR/bcftools_$(date +%Y%m%d_%H%M%S).log"

# Logging
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

usage() {
cat << EOF
USO: $0 -b <bam> -r <ref.fa> -o <out.vcf.gz> [ -t <threads> ] [ -c <config.txt> ]
ARGUMENTOS:
 -b, --bam      BAM de entrada (indexado)
 -r, --ref      Genoma referência FASTA
 -o, --out      VCF de saída (.vcf.gz recomendado)
 -t, --threads  Número de threads (default: 8)
 -c, --config   Arquivo extra de config com opções adicionais
 -h, --help     Mostra esta mensagem

EXEMPLO: $0 -b s.bam -r genome.fa -o s.vcf.gz -t 16
EOF
}

# Args
BAM=""; REF=""; OUT="";
while [[ $# -gt 0 ]]; do
  case $1 in
    -b|--bam) BAM="$2"; shift 2 ;;
    -r|--ref) REF="$2"; shift 2 ;;
    -o|--out) OUT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -c|--config) CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) log "Argumento desconhecido: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$BAM" || -z "$REF" || -z "$OUT" ]] && log "Campos obrigatórios ausentes!"; usage; exit 1
[[ ! -f "$BAM" ]] && log "Arquivo BAM não encontrado: $BAM"; exit 1
[[ ! -f "${BAM}.bai" && ! -f "${BAM%.*}.bai" ]] && log "Arquivo .bai obrigatório."; exit 1
[[ ! -f "$REF" ]] && log "Referência não encontrada: $REF"; exit 1

mkdir -p "$LOG_DIR"

EXTRA_ARGS=""
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
  while IFS= read -r line; do
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue
    EXTRA_ARGS="$EXTRA_ARGS $line"
  done < "$CONFIG_FILE"
fi

log "Iniciando bcftools: $(date)"
log "BAM: $BAM | REF: $REF | SAÍDA: $OUT | THREADS: $THREADS"
log "Config extra: $CONFIG_FILE"

CMD="bcftools mpileup -Ou -f $REF -a AD,DP,SP -q 20 -Q 20 -d 250 $BAM | \
     bcftools call -mv -Oz -o $OUT $EXTRA_ARGS"

log "Comando: $CMD"

if eval "$CMD" 2>>"$LOG_FILE"; then
  log "Chamada bcftools finalizada com sucesso. VCF: $OUT"
  bcftools index "$OUT"
  log "Index criado para $OUT"
else
  log "Erro ao executar bcftools. Consulte log: $LOG_FILE"
  exit 1
fi

log "Fim: $(date)"
exit 0
