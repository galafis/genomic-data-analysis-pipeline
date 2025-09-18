"""
Módulo de Métricas de Qualidade - Pipeline de Análise de Dados Genômicos

Este módulo fornece funcionalidades para calcular métricas básicas de qualidade
para dados de sequenciamento genômico, incluindo análise de qualidade de bases,
conteúdo GC e detecção de sequências duplicadas.

Autor: Equipe de Desenvolvimento
Versão: 1.0.0
Data: Setembro 2025
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from collections import Counter


def calcular_q30(qualidades: List[int]) -> float:
    """
    Calcula a porcentagem de bases com qualidade Q30 ou superior.
    
    Q30 indica uma probabilidade de erro de 1 em 1000 (99,9% de precisão).
    
    Args:
        qualidades (List[int]): Lista de valores de qualidade das bases
    
    Returns:
        float: Porcentagem de bases com Q30+
    
    Exemplo:
        >>> qualidades = [35, 28, 40, 32, 25, 38]
        >>> calcular_q30(qualidades)
        66.67
    """
    if not qualidades:
        return 0.0
    
    q30_bases = sum(1 for q in qualidades if q >= 30)
    return (q30_bases / len(qualidades)) * 100


def calcular_conteudo_gc(sequencia: str) -> float:
    """
    Calcula o conteúdo GC de uma sequência de DNA.
    
    Args:
        sequencia (str): Sequência de DNA (A, T, G, C)
    
    Returns:
        float: Porcentagem de conteúdo GC
    
    Exemplo:
        >>> sequencia = "ATGCGCTA"
        >>> calcular_conteudo_gc(sequencia)
        50.0
    """
    if not sequencia:
        return 0.0
    
    sequencia_upper = sequencia.upper()
    gc_count = sequencia_upper.count('G') + sequencia_upper.count('C')
    return (gc_count / len(sequencia_upper)) * 100


def detectar_duplicatas(sequencias: List[str]) -> Dict[str, any]:
    """
    Detecta e calcula a taxa de sequências duplicadas.
    
    Args:
        sequencias (List[str]): Lista de sequências para análise
    
    Returns:
        Dict[str, any]: Dicionário contendo:
            - 'taxa_duplicatas': Porcentagem de sequências duplicadas
            - 'total_sequencias': Número total de sequências
            - 'sequencias_unicas': Número de sequências únicas
            - 'duplicatas_encontradas': Número de duplicatas
    
    Exemplo:
        >>> sequencias = ["ATGC", "GCTA", "ATGC", "TTAA"]
        >>> resultado = detectar_duplicatas(sequencias)
        >>> print(resultado['taxa_duplicatas'])
        25.0
    """
    if not sequencias:
        return {
            'taxa_duplicatas': 0.0,
            'total_sequencias': 0,
            'sequencias_unicas': 0,
            'duplicatas_encontradas': 0
        }
    
    contador = Counter(sequencias)
    total_sequencias = len(sequencias)
    sequencias_unicas = len(contador)
    duplicatas = total_sequencias - sequencias_unicas
    taxa_duplicatas = (duplicatas / total_sequencias) * 100
    
    return {
        'taxa_duplicatas': taxa_duplicatas,
        'total_sequencias': total_sequencias,
        'sequencias_unicas': sequencias_unicas,
        'duplicatas_encontradas': duplicatas
    }


def analisar_distribuicao_qualidade(qualidades: List[int]) -> Dict[str, float]:
    """
    Analisa a distribuição de qualidades e calcula estatísticas descritivas.
    
    Args:
        qualidades (List[int]): Lista de valores de qualidade
    
    Returns:
        Dict[str, float]: Estatísticas de qualidade (média, mediana, etc.)
    
    Exemplo:
        >>> qualidades = [30, 35, 32, 28, 40, 38]
        >>> stats = analisar_distribuicao_qualidade(qualidades)
        >>> print(f"Qualidade média: {stats['media']:.2f}")
        Qualidade média: 33.83
    """
    if not qualidades:
        return {
            'media': 0.0,
            'mediana': 0.0,
            'desvio_padrao': 0.0,
            'minimo': 0.0,
            'maximo': 0.0
        }
    
    array_qualidades = np.array(qualidades)
    
    return {
        'media': float(np.mean(array_qualidades)),
        'mediana': float(np.median(array_qualidades)),
        'desvio_padrao': float(np.std(array_qualidades)),
        'minimo': float(np.min(array_qualidades)),
        'maximo': float(np.max(array_qualidades))
    }


def gerar_relatorio_qualidade(dados_fastq: Dict[str, any]) -> str:
    """
    Gera um relatório resumido das métricas de qualidade.
    
    Args:
        dados_fastq (Dict[str, any]): Dados processados do arquivo FASTQ
    
    Returns:
        str: Relatório formatado das métricas de qualidade
    
    Exemplo de uso programático:
        >>> # Exemplo de integração com pipeline principal
        >>> from quality_metrics import *
        >>> 
        >>> # Simulação de dados de entrada
        >>> dados_exemplo = {
        ...     'sequencias': ['ATGCGCTA', 'GCTAGCTA', 'ATGCGCTA'],
        ...     'qualidades': [35, 32, 40, 38, 28, 36, 30, 34]
        ... }
        >>> 
        >>> # Calcular métricas individuais
        >>> q30_percent = calcular_q30(dados_exemplo['qualidades'])
        >>> gc_content = np.mean([calcular_conteudo_gc(seq) for seq in dados_exemplo['sequencias']])
        >>> duplicatas_info = detectar_duplicatas(dados_exemplo['sequencias'])
        >>> stats_qualidade = analisar_distribuicao_qualidade(dados_exemplo['qualidades'])
        >>> 
        >>> # Gerar relatório
        >>> relatorio = gerar_relatorio_qualidade({
        ...     'q30_percent': q30_percent,
        ...     'gc_content': gc_content,
        ...     'duplicatas': duplicatas_info,
        ...     'stats_qualidade': stats_qualidade
        ... })
        >>> print(relatorio)
    """
    relatorio = """
╔══════════════════════════════════════════════════════════════╗
║              RELATÓRIO DE QUALIDADE - DADOS GENÔMICOS       ║
╠══════════════════════════════════════════════════════════════╣
║ Métricas de Qualidade de Bases:                             ║
║   • Q30+: {q30:.2f}%                                        ║
║   • Qualidade Média: {media_qual:.2f}                       ║
║   • Qualidade Mediana: {mediana_qual:.2f}                   ║
║                                                              ║
║ Análise de Composição:                                       ║
║   • Conteúdo GC: {gc:.2f}%                                  ║
║                                                              ║
║ Detecção de Duplicatas:                                      ║
║   • Taxa de Duplicatas: {dup_rate:.2f}%                     ║
║   • Sequências Totais: {total_seq}                          ║
║   • Sequências Únicas: {unique_seq}                         ║
╚══════════════════════════════════════════════════════════════╝
    """.format(
        q30=dados_fastq.get('q30_percent', 0),
        media_qual=dados_fastq.get('stats_qualidade', {}).get('media', 0),
        mediana_qual=dados_fastq.get('stats_qualidade', {}).get('mediana', 0),
        gc=dados_fastq.get('gc_content', 0),
        dup_rate=dados_fastq.get('duplicatas', {}).get('taxa_duplicatas', 0),
        total_seq=dados_fastq.get('duplicatas', {}).get('total_sequencias', 0),
        unique_seq=dados_fastq.get('duplicatas', {}).get('sequencias_unicas', 0)
    )
    
    return relatorio


if __name__ == "__main__":
    """
    Exemplo de execução direta para testes.
    """
    # Dados de exemplo para demonstração
    sequencias_teste = [
        "ATGCGCTACGATCGATC",
        "GCTAGCTACGATCGATC", 
        "ATGCGCTACGATCGATC",  # Duplicata
        "TTAAGCTACGATCGATC"
    ]
    
    qualidades_teste = [35, 32, 40, 38, 28, 36, 30, 34, 42, 29, 31, 37, 33, 39, 26, 41]
    
    print("=== DEMONSTRAÇÃO DO MÓDULO QUALITY_METRICS ===")
    print(f"Q30+: {calcular_q30(qualidades_teste):.2f}%")
    print(f"Conteúdo GC médio: {np.mean([calcular_conteudo_gc(seq) for seq in sequencias_teste]):.2f}%")
    print(f"Duplicatas: {detectar_duplicatas(sequencias_teste)}")
    print(f"Estatísticas de qualidade: {analisar_distribuicao_qualidade(qualidades_teste)}")
