#!/usr/bin/env python3
"""
MultiQC Reporter - Agregador de Relatórios FastQC
Este módulo fornece uma interface Python para agregar relatórios de controle
de qualidade usando MultiQC em dados de sequenciamento NGS.

NOTA: Este módulo depende do MultiQC. Para instalar:
    pip install multiqc

Autor: Genomic Data Analysis Pipeline
Versão: 1.0
Data: 2025-09-18
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from typing import List, Optional, Union
import logging

# Configuração de logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MultiQCReporter:
    """
    Wrapper para MultiQC que agrega relatórios de controle de qualidade.
    
    Esta classe encapsula as funcionalidades do MultiQC, permitindo execução
    programática e agregação de múltiplos relatórios de ferramentas de QC.
    
    DEPENDÊNCIA: Requer instalação do MultiQC (pip install multiqc)
    """
    
    def __init__(self, multiqc_path: str = 'multiqc'):
        """
        Inicializa o wrapper do MultiQC.
        
        Args:
            multiqc_path: Caminho para o executável do MultiQC
        """
        self.multiqc_path = multiqc_path
        self._check_multiqc_installation()
    
    def _check_multiqc_installation(self) -> None:
        """
        Verifica se o MultiQC está instalado e acessível.
        """
        try:
            result = subprocess.run([self.multiqc_path, '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"MultiQC encontrado: {result.stdout.strip()}")
            else:
                raise FileNotFoundError("MultiQC não encontrado")
        except FileNotFoundError:
            logger.error("MultiQC não está instalado ou não está no PATH")
            logger.error("Para instalar: pip install multiqc")
            sys.exit(1)
    
    def generate_report(self, 
                       input_dir: str,
                       output_dir: str = './multiqc_report',
                       report_title: str = 'Quality Control Report',
                       force: bool = True,
                       zip_data_dir: bool = False) -> None:
        """
        Gera relatório agregado MultiQC.
        
        Args:
            input_dir: Diretório contendo relatórios FastQC e outras ferramentas
            output_dir: Diretório de saída para o relatório
            report_title: Título do relatório
            force: Sobrescrever relatórios existentes
            zip_data_dir: Compactar diretório de dados
        """
        # Criar diretório de saída
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Construir comando
        cmd = [self.multiqc_path]
        cmd.extend(['--outdir', output_dir])
        cmd.extend(['--title', report_title])
        
        if force:
            cmd.append('--force')
        
        if zip_data_dir:
            cmd.append('--zip-data-dir')
        
        cmd.append(input_dir)
        
        # Executar comando
        logger.info(f"Executando: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("MultiQC executado com sucesso")
            logger.info(f"Relatório salvo em: {output_dir}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Erro na execução do MultiQC: {e}")
            logger.error(f"Saída de erro: {e.stderr}")
            raise
    
    def aggregate_fastqc_reports(self,
                                fastqc_dir: str,
                                output_dir: str = './aggregated_qc_report') -> None:
        """
        Agrega especificamente relatórios FastQC.
        
        Args:
            fastqc_dir: Diretório contendo relatórios FastQC
            output_dir: Diretório de saída
        """
        logger.info(f"Agregando relatórios FastQC de: {fastqc_dir}")
        self.generate_report(
            input_dir=fastqc_dir,
            output_dir=output_dir,
            report_title='FastQC Aggregated Report'
        )

def main():
    """
    Função principal para execução via linha de comando.
    """
    parser = argparse.ArgumentParser(
        description='MultiQC Reporter - Agregador de Relatórios de Controle de Qualidade'
    )
    
    parser.add_argument('input_dir', 
                       help='Diretório contendo relatórios de QC')
    parser.add_argument('-o', '--output', default='./multiqc_report',
                       help='Diretório de saída (padrão: ./multiqc_report)')
    parser.add_argument('-t', '--title', default='Quality Control Report',
                       help='Título do relatório')
    parser.add_argument('--no-force', action='store_true',
                       help='Não sobrescrever relatórios existentes')
    parser.add_argument('--zip-data', action='store_true',
                       help='Compactar diretório de dados')
    
    args = parser.parse_args()
    
    # Inicializar reporter
    reporter = MultiQCReporter()
    
    try:
        reporter.generate_report(
            args.input_dir,
            args.output,
            args.title,
            not args.no_force,
            args.zip_data
        )
        
    except Exception as e:
        logger.error(f"Erro durante execução: {e}")
        sys.exit(1)

if __name__ == '__main__':
    # Exemplo de uso básico
    print("""\nExemplo de uso do MultiQC Reporter:\n
# Uso via linha de comando:
# python multiqc_reporter.py /caminho/para/relatorios/fastqc -o ./meu_relatorio

# Uso programático:
from multiqc_reporter import MultiQCReporter

# Inicializar
reporter = MultiQCReporter()

# Agregar relatórios FastQC
reporter.aggregate_fastqc_reports(
    fastqc_dir='/caminho/para/relatorios/fastqc',
    output_dir='./relatorio_agregado'
)
""")
    
    main()
