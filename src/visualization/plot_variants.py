#!/usr/bin/env python3
"""
plot_variants.py — Visualização de variantes a partir de um arquivo VCF

Autor: Gabriel
Descrição:
  Lê um arquivo VCF e gera visualizações:
    - Manhattan plot (−log10 p-valor se presente em INFO/FORMAT; caso contrário, usa QUAL como proxy)
    - Histograma da pontuação QUAL
    - Gráfico de profundidade de leitura (DP) por variante

Saídas:
  Salva os gráficos em PNG e SVG no diretório de saída especificado.

Uso:
  python src/visualization/plot_variants.py \
    --vcf data/variants.vcf.gz \
    --outdir results/plots \
    --prefix sample1 \
    [--qual-col QUAL] [--dp-col DP] [--p-col PVAL]

Dependências:
  - Python 3.8+
  - pandas, numpy, matplotlib, seaborn
  - gzip (para VCF.gz) — leitura transparente via Python

Exemplos:
  python src/visualization/plot_variants.py --vcf data/example.vcf.gz --outdir results/plots --prefix ex1

Observações:
  - O script faz um parse leve do VCF (sem pyvcf/cyvcf2) para evitar dependências extras,
    lendo cabeçalho e colunas CHROM, POS, ID, REF, ALT, QUAL, INFO.
  - Para o Manhattan plot, se houver coluna/INFO com p-valor (p.ex. P, PV, PVAL), utilize --p-col.
  - Para DP (profundidade), o script tenta extrair INFO; caso não encontre, tenta FORMAT/DP da primeira amostra.
"""

from __future__ import annotations
import argparse
import gzip
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(context="notebook", style="whitegrid")

# ----------------------------
# Utilidades de parsing de VCF
# ----------------------------
VCF_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

def _open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def _extract_from_info(info: str, key: str) -> Optional[str]:
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            if k == key:
                return v
        else:
            if field == key:
                return "1"
    return None

def _guess_p_keys(info_header: List[str]) -> List[str]:
    # heurística: campos contendo 'P' comum em p-valores
    candidates = []
    for k in info_header:
        if re.search(r"(^|_)p(val)?($|_|\b)", k, re.IGNORECASE) or k.upper() in {"P", "PV", "PVALUE", "PVAL"}:
            candidates.append(k)
    return candidates

def read_vcf_light(vcf_path: str, dp_key: str = "DP", p_key: Optional[str] = None,
                   max_records: Optional[int] = None) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Lê VCF de forma leve retornando DataFrame com colunas principais e metadados.

    Retorna:
      df: colunas [CHROM, POS, ID, QUAL, INFO, (DP?), (PVAL?)]
      meta: dicionário com informações úteis (e.g., amostras)
    """
    chroms: List[str] = []
    poss: List[int] = []
    ids: List[str] = []
    quals: List[float] = []
    infos: List[str] = []
    dps: List[Optional[float]] = []
    pvals: List[Optional[float]] = []

    meta: Dict[str, str] = {}
    samples: List[str] = []
    info_keys_from_header: List[str] = []

    with _open_maybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith("##INFO="):
                # tenta capturar a KEY entre <ID=KEY,
                m = re.search(r"<ID=([^,>]+)", line)
                if m:
                    info_keys_from_header.append(m.group(1))
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                header_map = {c: i for i, c in enumerate(parts)}
                if len(parts) > 8:
                    samples = parts[9:]
                meta["samples"] = ",".join(samples)
                continue
            # linhas de variantes
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            vid = parts[2]
            qual_str = parts[5]
            try:
                qual = float(qual_str) if qual_str != "." else np.nan
            except ValueError:
                qual = np.nan
            info = parts[7]

            # DP: primeiro tenta INFO; depois tenta FORMAT/ amostra 0
            dp_val: Optional[float] = None
            dp_info = _extract_from_info(info, dp_key)
            if dp_info is not None:
                try:
                    dp_val = float(dp_info)
                except ValueError:
                    dp_val = np.nan
            elif len(parts) > 8 and len(samples) > 0:
                fmt = parts[8].split(":")
                sample0 = parts[9].split(":") if len(parts) > 9 else []
                if dp_key in fmt:
                    idx = fmt.index(dp_key)
                    if idx < len(sample0):
                        try:
                            dp_val = float(sample0[idx]) if sample0[idx] != "." else np.nan
                        except ValueError:
                            dp_val = np.nan

            # PVAL: se informado via p_key
            p_val: Optional[float] = None
            if p_key:
                p_info = _extract_from_info(info, p_key)
                if p_info is not None:
                    try:
                        p_val = float(p_info)
                    except ValueError:
                        p_val = np.nan

            chroms.append(chrom)
            poss.append(pos)
            ids.append(vid)
            quals.append(qual)
            infos.append(info)
            dps.append(dp_val)
            pvals.append(p_val)

            if max_records is not None and len(chroms) >= max_records:
                break

    df = pd.DataFrame({
        "CHROM": chroms,
        "POS": poss,
        "ID": ids,
        "QUAL": quals,
        "INFO": infos,
        "DP": dps,
        "PVAL": pvals,
    })

    # tenta preencher chrs como categoria ordenada por número quando possível
    def _chr_key(c: str):
        m = re.match(r"^chr?(\d+)$", str(c), re.IGNORECASE)
        if m:
            return (0, int(m.group(1)))
        if str(c).lower() in {"x", "chrx"}: return (1, 23)
        if str(c).lower() in {"y", "chry"}: return (1, 24)
        return (2, str(c))

    chrom_order = sorted(df["CHROM"].dropna().unique(), key=_chr_key)
    df["CHROM"] = pd.Categorical(df["CHROM"], categories=chrom_order, ordered=True)

    # Se p_key não fornecido, sugere possíveis
    meta["pval_candidates"] = ",".join(_guess_p_keys(info_keys_from_header))
    return df, meta

# ----------------------------
# Plotters
# ----------------------------

def manhattan_plot(df: pd.DataFrame, out_base: str, title: str, use_p: bool) -> None:
    tmp = df.dropna(subset=["POS"]).copy()
    tmp = tmp.sort_values(["CHROM", "POS"])  # ordena por cromossomo e posição

    # posição cumulativa por cromossomo
    chr_offsets: Dict[str, int] = {}
    cum_pos = 0
    xticks = []
    xticklabels = []
    for chrom in tmp["CHROM"].cat.categories:
        chr_df = tmp[tmp["CHROM"] == chrom]
        if chr_df.empty:
            continue
        chr_offsets[str(chrom)] = cum_pos
        mid = cum_pos + (chr_df["POS"].max() - chr_df["POS"].min()) / 2
        xticks.append(mid)
        xticklabels.append(str(chrom))
        cum_pos += chr_df["POS"].max()

    tmp["cum_pos"] = tmp.apply(lambda r: r["POS"] + chr_offsets.get(str(r["CHROM"]), 0), axis=1)

    if use_p and ("PVAL" in tmp.columns) and tmp["PVAL"].notna().any():
        y = -np.log10(tmp["PVAL"].astype(float).replace(0, np.nan))
        ylab = "-log10(P)"
    else:
        y = tmp["QUAL"].astype(float)
        ylab = "QUAL"

    plt.figure(figsize=(12, 5))
    sns.scatterplot(x="cum_pos", y=y, data=tmp, hue="CHROM", palette="tab20", legend=False, s=8, edgecolor=None)
    plt.xlabel("Genomic position (cumulative)")
    plt.ylabel(ylab)
    plt.title(title)
    plt.xticks(xticks, xticklabels, rotation=90)
    plt.tight_layout()
    plt.savefig(out_base + "_manhattan.png", dpi=200)
    plt.savefig(out_base + "_manhattan.svg")
    plt.close()


def qual_histogram(df: pd.DataFrame, out_base: str, title: str) -> None:
    plt.figure(figsize=(8, 5))
    sns.histplot(df["QUAL"].dropna().astype(float), bins=50, kde=True, color="#2a9d8f")
    plt.xlabel("QUAL")
    plt.ylabel("Count")
    plt.title(f"QUAL distribution — {title}")
    plt.tight_layout()
    plt.savefig(out_base + "_qual_hist.png", dpi=200)
    plt.savefig(out_base + "_qual_hist.svg")
    plt.close()


def depth_plot(df: pd.DataFrame, out_base: str, title: str) -> None:
    if "DP" not in df.columns or df["DP"].dropna().empty:
        return
    plt.figure(figsize=(10, 5))
    sns.scatterplot(x="POS", y="DP", data=df.dropna(subset=["DP"]).astype({"DP": float}), s=8, color="#264653")
    plt.xlabel("POS")
    plt.ylabel("Depth (DP)")
    plt.title(f"Depth by position — {title}")
    plt.tight_layout()
    plt.savefig(out_base + "_depth.png", dpi=200)
    plt.savefig(out_base + "_depth.svg")
    plt.close()

# ----------------------------
# CLI
# ----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Gera gráficos (Manhattan, QUAL, DP) a partir de VCF.")
    p.add_argument("--vcf", required=True, help="Caminho para o arquivo VCF (.vcf ou .vcf.gz)")
    p.add_argument("--outdir", required=True, help="Diretório de saída para salvar imagens")
    p.add_argument("--prefix", required=True, help="Prefixo do arquivo de saída")
    p.add_argument("--dp-col", default="DP", help="Nome do campo para profundidade (INFO/FORMAT), padrão: DP")
    p.add_argument("--p-col", default=None, help="Nome do campo INFO com p-valor (opcional)")
    p.add_argument("--title", default=None, help="Título dos gráficos (opcional)")
    p.add_argument("--max-records", type=int, default=None, help="Limitar nº de variantes lidas (debug)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out_base = os.path.join(args.outdir, args.prefix)

    df, meta = read_vcf_light(args.vcf, dp_key=args.dp_col, p_key=args.p_col, max_records=args.max_records)

    title = args.title or f"VCF: {os.path.basename(args.vcf)}"
    if args.p_col is None and meta.get("pval_candidates"):
        print(f"[Info] Possíveis chaves de p-valor no INFO: {meta['pval_candidates']}")

    # Plots
    manhattan_plot(df, out_base, title, use_p=args.p_col is not None)
    qual_histogram(df, out_base, title)
    depth_plot(df, out_base, title)

    print(f"[OK] Gráficos salvos em: {out_base}_*.png/svg")


if __name__ == "__main__":
    main()
