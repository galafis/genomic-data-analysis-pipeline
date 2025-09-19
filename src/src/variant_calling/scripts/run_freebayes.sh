#!/bin/bash
# SCRIPT: run_freebayes.sh
# AUTHOR: Genomic Data Analysis Pipeline Team
# VERSION: 1.0
# DATE: 2025-09-19
# DESCRIPTION: Pipeline profissional para chamada de variantes germinativas com FreeBayes.
# USO:
#   ./run_freebayes.sh -b input.bam -r reference.fa -o output.vcf [ -t threads ] [ -c config.txt ]

set -euo pipefail

# Defaults
THREADS=8
CONFIG_FILE=""
LOG_DIR="./logs"
LOG_FILE="$LOG_DIR/freebayes_$(date +%Y%m%d_%H%M%S).log"

# Funções de logging
log()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

usage() {
cat << EOF
USAGE: $0 -b <bam> -r <ref.fa> -o <out.vcf> [ -t <threads> ] [ -c <config.txt> ]
ARGUMENTOS:
 -b, --bam         BAM de entrada (sorted & indexed)
 -r, --ref         genoma de referência FASTA
 -o, --out         VCF de saída (recomendado .vcf)
 -t, --threads     número de threads (default: 8)
 -c, --config      arquivo de config extra (opcional)
 -h, --help        mostra esta mensagem

EXEMPLO: $0 -b sample.bam -r hg38.fa -o variants.vcf -t 16
EOF
}

# Parse args
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

[[ -z "$BAM" || -z "$REF" || -z "$OUT" ]] && log "BAM, referência e saída são obrigatórios!" && usage && exit 1
[[ ! -f "$BAM" ]] && log "Arquivo BAM não encontrado: $BAM" && exit 1
[[ ! -f "${BAM}.bai" && ! -f "${BAM%.*}.bai" ]] && log "Arquivo BAM não indexado (.bai ausente)." && exit 1
[[ ! -f "$REF" ]] && log "Arquivo de referência inexistente: $REF" && exit 1

mkdir -p "$LOG_DIR"

EXTRA_ARGS=""
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
  while IFS= read -r line; do
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue
    EXTRA_ARGS="$EXTRA_ARGS $line"
  done < "$CONFIG_FILE"
fi

log "Iniciando FreeBayes: $(date)"
log "BAM: $BAM | REF: $REF | SAIDA: $OUT | THREADS: $THREADS | CONFIG: $CONFIG_FILE"

CMD="freebayes -f $REF --bam $BAM --vcf $OUT --use-best-n-alleles 4 --min-coverage 10 --min-alternate-fraction 0.2 --min-base-quality 20 --min-mapping-quality 20 --threads $THREADS $EXTRA_ARGS"

log "Comando: $CMD"

if eval "$CMD" 2>> "$LOG_FILE"; then
  log "Chamada FreeBayes concluída com sucesso! VCF em: $OUT"
else
  log "Erro na chamada FreeBayes. Consulte log: $LOG_FILE"
  exit 1
fi

log "Fim da execução: $(date)"
exit 0
