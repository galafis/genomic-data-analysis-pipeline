#!/bin/bash
# SCRIPT: run_variant_filter.sh
# AUTHOR: Genomic Data Analysis Pipeline Team
# VERSION: 1.0
# DATE: 2025-09-19
# DESCRIPTION: Pipeline robusto para filtragem de variantes com bcftools.
# USO:
#   ./run_variant_filter.sh -v input.vcf.gz -o output.filtered.vcf.gz -q 30 -d 10 -m 40 [ -c config.txt ]

set -euo pipefail

# Defaults
QUAL=30
DP=10
MQ=40
CONFIG_FILE=""
LOG_DIR="./logs"
LOG_FILE="$LOG_DIR/filter_$(date +%Y%m%d_%H%M%S).log"

# Logging
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

usage() {
cat << EOF
USO: $0 -v <input.vcf.gz> -o <output.vcf.gz> [ -q <QUAL> -d <DP> -m <MQ> ] [ -c <config.txt> ]
- -v, --vcf     VCF de entrada (.gz indexado)
- -o, --out     VCF de saída filtrado
- -q, --qual    Qualidade mínima (QUAL, default: 30)
- -d, --dp      Profundidade mínima (DP, default: 10)
- -m, --mq      Qualidade de mapeamento mínima (MQ, default: 40)
- -c, --config  Arquivo de config extra (opcional, bcftools filter)
- -h, --help    Ajuda

Exemplo:
  $0 -v sample.vcf.gz -o sample.filtered.vcf.gz -q 30 -d 10 -m 40
EOF
}

# Args
IN_VCF=""; OUT_VCF="";
while [[ $# -gt 0 ]]; do
  case $1 in
    -v|--vcf) IN_VCF="$2"; shift 2 ;;
    -o|--out) OUT_VCF="$2"; shift 2 ;;
    -q|--qual) QUAL="$2"; shift 2 ;;
    -d|--dp) DP="$2"; shift 2 ;;
    -m|--mq) MQ="$2"; shift 2 ;;
    -c|--config) CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) log "Argumento desconhecido: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$IN_VCF" || -z "$OUT_VCF" ]] && log "VCF de entrada e saída são obrigatórios!" && usage && exit 1
[[ ! -f "$IN_VCF" ]] && log "Arquivo VCF não encontrado: $IN_VCF" && exit 1
[[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]] && log "Index do VCF (.tbi/.csi) não encontrado!"; exit 1

mkdir -p "$LOG_DIR"

FILTER_EXPR="QUAL>=$QUAL && DP>=$DP && MQ>=$MQ"
EXTRA_ARGS=""
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
  while IFS= read -r line; do
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue
    EXTRA_ARGS="$EXTRA_ARGS $line"
  done < "$CONFIG_FILE"
fi

log "Filtrando VCF: $IN_VCF | QUAL>=$QUAL, DP>=$DP, MQ>=$MQ | Config: $CONFIG_FILE"
CMD="bcftools filter -i '$FILTER_EXPR' $EXTRA_ARGS -Oz -o $OUT_VCF $IN_VCF"
log "Comando: $CMD"

if eval "$CMD" 2>>"$LOG_FILE"; then
  bcftools index "$OUT_VCF"
  log "VCF filtrado criado: $OUT_VCF"
else
  log "Erro ao filtrar VCF. Verifique: $LOG_FILE"
  exit 1
fi

log "Execução finalizada: $(date)"
exit 0
