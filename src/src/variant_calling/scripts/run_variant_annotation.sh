#!/bin/bash
# SCRIPT: run_variant_annotation.sh
# AUTHOR: Genomic Data Analysis Pipeline Team
# VERSION: 1.0
# DATE: 2025-09-19
# DESCRIPTION: Pipeline para anotação de variantes via VEP ou SnpEff.
# USAGE:
#   ./run_variant_annotation.sh -v input.vcf.gz -o output.ann.vcf.gz -a vep|snpeff [ -r ref.fa|db ] [ -c config ]

set -euo pipefail

ANNOTATOR="vep"
CONFIG_FILE=""
LOG_DIR="./logs"
LOG_FILE="$LOG_DIR/annot_$(date +%Y%m%d_%H%M%S).log"

# Logging
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

usage() {
cat << EOF
USO: $0 -v <vcf.gz> -o <output.vcf.gz> -a vep|snpeff [ -r <ref/db> ] [ -c <config> ]
- -v, --vcf       VCF comprimido (.vcf.gz, indexado)
- -o, --out       Arquivo de saída anotado
- -a, --annot     Anotador ('vep' ou 'snpeff', default: vep)
- -r, --ref       FASTA de referência (vep) ou nome do DB (snpeff)
- -c, --config    Config extra para chamada do anotador
- -h, --help      Exibe esta ajuda

Exemplo: $0 -v sample.filtered.vcf.gz -o sample.ann.vcf.gz -a vep -r genome.fa
Exemplo: $0 -v sample.vcf.gz -o sample.ann.vcf.gz -a snpeff -r GRCh38.86
EOF
}

IN_VCF=""; OUT_VCF=""; REF="";
while [[ $# -gt 0 ]]; do
  case $1 in
    -v|--vcf) IN_VCF="$2"; shift 2 ;;
    -o|--out) OUT_VCF="$2"; shift 2 ;;
    -a|--annot) ANNOTATOR="$2"; shift 2 ;;
    -r|--ref) REF="$2"; shift 2 ;;
    -c|--config) CONFIG_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) log "Argumento desconhecido: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$IN_VCF" || -z "$OUT_VCF" ]] && log "VCF de entrada e saída são obrigatórios!"; usage; exit 1
[[ ! -f "$IN_VCF" ]] && log "VCF não encontrado: $IN_VCF"; exit 1
[[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]] && log "VCF deve estar indexado (.tbi /.csi)"; exit 1

mkdir -p "$LOG_DIR"

EXTRA_ARGS=""
if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
  while IFS= read -r line; do
    [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue
    EXTRA_ARGS="$EXTRA_ARGS $line"
  done < "$CONFIG_FILE"
fi

if [[ "$ANNOTATOR" == "vep" ]]; then
  [[ -z "$REF" ]] && log "FASTA de referência é obrigatório para VEP." && exit 1
  CMD="vep --input_file $IN_VCF --output_file $OUT_VCF --vcf --fasta $REF --cache --offline $EXTRA_ARGS"
elif [[ "$ANNOTATOR" == "snpeff" ]]; then
  [[ -z "$REF" ]] && log "Nome do DB SnpEff é obrigatório." && exit 1
  CMD="snpeff $EXTRA_ARGS $REF $IN_VCF > $OUT_VCF"
else
  log "Anotador não suportado: $ANNOTATOR"; exit 1
fi

log "Iniciando anotação: $ANNOTATOR | IN: $IN_VCF | OUT: $OUT_VCF | REF: $REF | Config: $CONFIG_FILE"
log "Comando: $CMD"

if eval "$CMD" 2>>"$LOG_FILE"; then
  log "Anotação finalizada com sucesso: $OUT_VCF"
else
  log "Falha na anotação. Consulte: $LOG_FILE"
  exit 1
fi

log "Execução concluída: $(date)"
exit 0
