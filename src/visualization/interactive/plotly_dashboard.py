#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotly/Dash Variant Browser (plotly_dashboard.py)

Author: Gabriel

Descrição
---------
Dashboard web interativo para explorar variantes a partir de um arquivo VCF.
Recursos principais:
- Upload/entrada de arquivo VCF (ou via linha de comando)
- Tabela interativa de variantes com filtros
- Gráfico dinâmico de frequência alélica (AF) por cromossomo/posição
- Histograma de QUAL
- Download de dados filtrados

Requisitos
----------
- Python >= 3.9
- dash, dash-bootstrap-components, pandas, plotly, cyvcf2 (ou pysam/pyvcf como alternativa)

Instalação rápida
-----------------
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install dash dash-bootstrap-components plotly pandas cyvcf2

Uso
---
1) Via linha de comando (aponta um VCF existente):
   python src/visualization/interactive/plotly_dashboard.py --vcf data/variants.vcf.gz --host 0.0.0.0 --port 8050

2) Sem argumentos: use o upload no topo da página para carregar um VCF local.

Notas
-----
- Para .vcf.gz indexado, recomenda-se bgzip + tabix. cyvcf2 lida bem com VCF comprimido.
- Colunas calculadas: CHROM, POS, REF, ALT, QUAL, FILTER, INFO_AF (quando disponível), DP (se disponível), ID.
- Se AF não estiver em INFO, tenta-se derivar a partir de AD (quando presente) para variantes bialélicas.
"""

from __future__ import annotations
import base64
import io
import os
import argparse
from typing import List, Dict, Any, Optional

import pandas as pd
import plotly.express as px
from dash import Dash, dcc, html, Input, Output, State, dash_table, callback_context
import dash_bootstrap_components as dbc

# Leitura de VCF com cyvcf2
try:
    from cyvcf2 import VCF
except ImportError:  # fallback leve orientativo
    VCF = None


def parse_vcf_to_df(vcf_path: str, max_records: Optional[int] = None) -> pd.DataFrame:
    if VCF is None:
        raise RuntimeError("cyvcf2 não instalado. Instale com: pip install cyvcf2")
    v = VCF(vcf_path)
    rows = []
    for i, rec in enumerate(v):
        if max_records is not None and i >= max_records:
            break
        chrom = rec.CHROM
        pos = rec.POS
        ref = rec.REF
        # ALT pode ter múltiplos; criar linhas separadas
        for alt in rec.ALT or [None]:
            qual = rec.QUAL
            filt = ";".join(rec.FILTER) if isinstance(rec.FILTER, list) else (rec.FILTER or "PASS")
            vid = rec.ID
            # INFO fields
            info_af = None
            if 'AF' in rec.INFO:
                try:
                    # AF pode ser lista
                    af_val = rec.INFO.get('AF')
                    if isinstance(af_val, (list, tuple)):
                        info_af = float(af_val[0]) if af_val else None
                    else:
                        info_af = float(af_val)
                except Exception:
                    info_af = None
            # Depth
            dp = None
            if 'DP' in rec.INFO:
                try:
                    dp_val = rec.INFO.get('DP')
                    dp = int(dp_val) if dp_val is not None else None
                except Exception:
                    dp = None
            # Se AF ausente, tentar AD por amostra 0
            if info_af is None and rec.genotypes:
                try:
                    # cyvcf2: rec.format('AD') retorna matriz (amostras x alleles)
                    ad = rec.format('AD')  # profundidade por alelo
                    if ad is not None and len(ad) > 0:
                        ref_count = ad[0][0] if len(ad[0]) > 0 else None
                        alt_count = ad[0][1] if len(ad[0]) > 1 else None
                        if ref_count is not None and alt_count is not None:
                            total = float(ref_count + alt_count)
                            info_af = (alt_count / total) if total > 0 else None
                except Exception:
                    pass
            rows.append({
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'QUAL': qual,
                'FILTER': filt,
                'ID': vid,
                'DP': dp,
                'INFO_AF': info_af,
            })
    df = pd.DataFrame(rows)
    # Tipos e colunas auxiliares
    if not df.empty:
        df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')
        df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
        # Cromossomos como categoria para ordenação comum
        chr_order = [str(x) for x in list(range(1, 23))] + ['X', 'Y', 'MT', 'M']
        df['CHROM'] = df['CHROM'].astype(str)
        df['CHROM'] = pd.Categorical(df['CHROM'], categories=chr_order, ordered=True)
    return df


def decode_uploaded_vcf(content: str, filename: str) -> str:
    """Salva conteúdo enviado (base64) como arquivo temporário e retorna o caminho."""
    content_type, content_string = content.split(',')
    decoded = base64.b64decode(content_string)
    # manter extensão
    suffix = ".vcf.gz" if filename.endswith(".vcf.gz") else ".vcf"
    tmp_path = os.path.join("/tmp", f"uploaded_{os.getpid()}_{filename}")
    with open(tmp_path, 'wb') as f:
        f.write(decoded)
    return tmp_path


def serve_layout() -> html.Div:
    return dbc.Container([
        html.H3("Variant Browser - Dash/Plotly"),
        dbc.Alert("Carregue um VCF via upload ou informe --vcf na linha de comando.", color="info"),
        dbc.Row([
            dbc.Col([
                dcc.Upload(
                    id='upload-vcf',
                    children=html.Div(['Arraste e solte ou ', html.A('selecione um VCF')]),
                    style={
                        'width': '100%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                        'textAlign': 'center', 'margin': '10px'
                    },
                    multiple=False
                ),
                dbc.Row([
                    dbc.Col([
                        dbc.Label('Filtro por cromossomo'),
                        dcc.Dropdown(id='chrom-filter', multi=True, placeholder='Selecione cromossomos')
                    ], md=4),
                    dbc.Col([
                        dbc.Label('AF mínima'),
                        dcc.Input(id='af-min', type='number', placeholder='0.0', debounce=True)
                    ], md=2),
                    dbc.Col([
                        dbc.Label('AF máxima'),
                        dcc.Input(id='af-max', type='number', placeholder='1.0', debounce=True)
                    ], md=2),
                    dbc.Col([
                        dbc.Label('QUAL mínima'),
                        dcc.Input(id='qual-min', type='number', placeholder='0', debounce=True)
                    ], md=2),
                    dbc.Col([
                        dbc.Label('Texto (REF/ALT/FILTER/ID)'),
                        dcc.Input(id='text-filter', type='text', placeholder='e.g., PASS, rsid, A>G', debounce=True)
                    ], md=2),
                ]),
                html.Br(),
                dash_table.DataTable(
                    id='variants-table',
                    columns=[{"name": c, "id": c} for c in ['CHROM','POS','REF','ALT','QUAL','FILTER','INFO_AF','DP','ID']],
                    page_size=15,
                    filter_action='native',
                    sort_action='native',
                    sort_mode='multi',
                    style_table={'overflowX': 'auto'},
                    style_cell={'fontFamily': 'monospace', 'fontSize': 12},
                    export_format='csv'
                ),
                html.Br(),
                dbc.Button("Baixar dados filtrados (CSV)", id='download-btn', color='primary'),
                dcc.Download(id="download-data")
            ], md=6),
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="AF por cromossomo", tab_id="tab-af"),
                    dbc.Tab(label="Histograma QUAL", tab_id="tab-qual"),
                ], id='tabs', active_tab='tab-af'),
                dcc.Graph(id='af-plot'),
                dcc.Graph(id='qual-hist')
            ], md=6)
        ]),
        # Armazena DataFrame processado em memória do navegador
        dcc.Store(id='vcf-store')
    ], fluid=True)


def filter_dataframe(df: pd.DataFrame,
                     chroms: Optional[List[str]],
                     af_min: Optional[float],
                     af_max: Optional[float],
                     qual_min: Optional[float],
                     text: Optional[str]) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    m = pd.Series([True] * len(df))
    if chroms:
        m &= df['CHROM'].astype(str).isin(chroms)
    if af_min is not None:
        m &= (df['INFO_AF'].fillna(-1) >= float(af_min))
    if af_max is not None:
        m &= (df['INFO_AF'].fillna(-1) <= float(af_max))
    if qual_min is not None:
        m &= (df['QUAL'].fillna(-1) >= float(qual_min))
    if text:
        t = str(text).lower()
        subset = df[['REF','ALT','FILTER','ID']].astype(str).apply(lambda s: s.str.lower())
        m &= subset.apply(lambda row: any(t in v for v in row.values), axis=1)
    return df[m].copy()


def build_af_plot(df: pd.DataFrame) -> Any:
    if df is None or df.empty:
        return px.scatter(title="Sem dados para AF")
    # Usar posição no eixo x e AF no eixo y
    fig = px.scatter(df, x='POS', y='INFO_AF', color='CHROM',
                     hover_data=['REF','ALT','QUAL','FILTER','ID'],
                     title='Frequência alélica (AF) por posição')
    fig.update_layout(legend_title_text='CHROM')
    fig.update_yaxes(range=[0,1])
    return fig


def build_qual_hist(df: pd.DataFrame) -> Any:
    if df is None or df.empty:
        return px.histogram(title="Sem dados para QUAL")
    fig = px.histogram(df, x='QUAL', nbins=50, title='Histograma de QUAL')
    return fig


def create_app(initial_df: Optional[pd.DataFrame] = None) -> Dash:
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
    app.title = "Variant Browser"
    app.layout = serve_layout

    @app.callback(
        Output('vcf-store', 'data'),
        Output('chrom-filter', 'options'),
        Input('upload-vcf', 'contents'),
        State('upload-vcf', 'filename'),
        prevent_initial_call=True
    )
    def on_upload(contents, filename):
        if contents is None:
            raise dash.exceptions.PreventUpdate
        tmp = decode_uploaded_vcf(contents, filename or 'uploaded.vcf')
        df = parse_vcf_to_df(tmp)
        chroms = sorted(df['CHROM'].dropna().astype(str).unique().tolist()) if not df.empty else []
        return df.to_dict('records'), [{'label': c, 'value': c} for c in chroms]

    @app.callback(
        Output('variants-table', 'data'),
        Output('af-plot', 'figure'),
        Output('qual-hist', 'figure'),
        Input('vcf-store', 'data'),
        Input('chrom-filter', 'value'),
        Input('af-min', 'value'),
        Input('af-max', 'value'),
        Input('qual-min', 'value'),
        Input('text-filter', 'value'),
        prevent_initial_call=True
    )
    def on_filter(data_records, chroms, af_min, af_max, qual_min, text):
        if data_records is None:
            raise dash.exceptions.PreventUpdate
        df = pd.DataFrame(data_records)
        fdf = filter_dataframe(df, chroms, af_min, af_max, qual_min, text)
        return fdf.to_dict('records'), build_af_plot(fdf), build_qual_hist(fdf)

    @app.callback(
        Output('download-data', 'data'),
        Input('download-btn', 'n_clicks'),
        State('vcf-store', 'data'),
        prevent_initial_call=True
    )
    def on_download(n_clicks, data_records):
        if not n_clicks or data_records is None:
            raise dash.exceptions.PreventUpdate
        df = pd.DataFrame(data_records)
        csv_bytes = df.to_csv(index=False).encode('utf-8')
        return dict(content=csv_bytes.decode('utf-8'), filename='variants_filtered.csv')

    # Se houver DF inicial (via --vcf), pré-carregar no store através de clientside callback simples
    if initial_df is not None and not initial_df.empty:
        # hack simples: usar uma Div oculta para injetar dados na primeira renderização
        app.layout = lambda: dbc.Container([
            dcc.Store(id='vcf-store', data=initial_df.to_dict('records')),
            serve_layout()
        ], fluid=True)

    return app


def main():
    parser = argparse.ArgumentParser(description="Dash/Plotly Variant Browser")
    parser.add_argument('--vcf', type=str, default=None, help='Caminho para arquivo VCF (.vcf ou .vcf.gz)')
    parser.add_argument('--host', type=str, default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8050)
    parser.add_argument('--max-records', type=int, default=None, help='Limite opcional de variantes para carregamento')
    args = parser.parse_args()

    initial_df = None
    if args.vcf:
        initial_df = parse_vcf_to_df(args.vcf, max_records=args.max_records)
    app = create_app(initial_df=initial_df)
    app.run_server(host=args.host, port=args.port, debug=True)


if __name__ == '__main__':
    main()
