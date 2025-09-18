#!/usr/bin/env python3
"""
FastQC Wrapper - Wrapper para FastQC

Este módulo fornece uma interface Python para executar análises de controle 
de qualidade usando FastQC em dados de sequenciamento NGS.

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

class FastQCWrapper:
    """
    Wrapper para FastQC que simplifica a execução de análises de controle de qualidade.
    
    Esta classe encapsula as funcionalidades do FastQC, permitindo execução
    programática e configuração flexível dos parâmetros.
    """
    
    def __init__(self, fastqc_path: str = 'fastqc'):
        """
        Inicializa o wrapper do FastQC.
        
        Args:
            fastqc_path: Caminho para o executável do FastQC
        """
        self.fastqc_path = fastqc_path
        self._check_fastqc_installation()
    
    def _check_fastqc_installation(self) -> None:
        """
        Verifica se o FastQC está instalado e acessível.
        """
        try:
            result = subprocess.run([self.fastqc_path, '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"FastQC encontrado: {result.stdout.strip()}")
            else:
                raise FileNotFoundError("FastQC não encontrado")
        except FileNotFoundError:
            logger.error("FastQC não está instalado ou não está no PATH")
            sys.exit(1)
    
    def run_fastqc(self, 
                   input_files: Union[str, List[str]], 
                   output_dir: str = './fastqc_output',
                   threads: int = 1,
                   quiet: bool = False,
                   extract: bool = True,
                   format_type: str = 'auto') -> None:
        """
        Executa análise FastQC nos arquivos de entrada.
        
        Args:
            input_files: Arquivo(s) de entrada (FASTQ/SAM/BAM)
            output_dir: Diretório de saída para os relatórios
            threads: Número de threads a usar
            quiet: Execução silenciosa
            extract: Extrair arquivos ZIP dos relatórios
            format_type: Formato dos arquivos de entrada
        """
        # Criar diretório de saída
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Preparar lista de arquivos
        if isinstance(input_files, str):
            input_files = [input_files]
        
        # Construir comando
        cmd = [self.fastqc_path]
        cmd.extend(['--outdir', output_dir])
        cmd.extend(['--threads', str(threads)])
        cmd.extend(['--format', format_type])
        
        if quiet:
            cmd.append('--quiet')
        
        if extract:
            cmd.append('--extract')
        
        cmd.extend(input_files)
        
        # Executar comando
        logger.info(f"Executando: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info("FastQC executado com sucesso")
            logger.info(f"Relatórios salvos em: {output_dir}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Erro na execução do FastQC: {e}")
            logger.error(f"Saída de erro: {e.stderr}")
            raise
    
    def batch_process(self, 
                     input_dir: str, 
                     output_dir: str = './fastqc_batch_output',
                     file_pattern: str = '*.fastq*',
                     threads: int = 4) -> None:
        """
        Processa múltiplos arquivos em lote.
        
        Args:
            input_dir: Diretório com arquivos de entrada
            output_dir: Diretório de saída
            file_pattern: Padrão de arquivos a processar
            threads: Número de threads
        """
        input_path = Path(input_dir)
        files = list(input_path.glob(file_pattern))
        
        if not files:
            logger.warning(f"Nenhum arquivo encontrado em {input_dir} com padrão {file_pattern}")
            return
        
        logger.info(f"Processando {len(files)} arquivos")
        
        file_paths = [str(f) for f in files]
        self.run_fastqc(file_paths, output_dir, threads=threads)

def main():
    """
    Função principal para execução via linha de comando.
    """
    parser = argparse.ArgumentParser(
        description='Wrapper Python para FastQC - Controle de Qualidade de Sequências'
    )
    
    parser.add_argument('input', nargs='+', 
                       help='Arquivo(s) de entrada ou diretório')
    parser.add_argument('-o', '--output', default='./fastqc_output',
                       help='Diretório de saída (padrão: ./fastqc_output)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Número de threads (padrão: 1)')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Execução silenciosa')
    parser.add_argument('--no-extract', action='store_true',
                       help='Não extrair arquivos ZIP')
    parser.add_argument('-f', '--format', default='auto',
                       choices=['auto', 'fastq', 'sam', 'bam'],
                       help='Formato dos arquivos de entrada')
    parser.add_argument('--batch', action='store_true',
                       help='Modo processamento em lote (input deve ser diretório)')
    parser.add_argument('--pattern', default='*.fastq*',
                       help='Padrão de arquivos para modo lote')
    
    args = parser.parse_args()
    
    # Inicializar wrapper
    wrapper = FastQCWrapper()
    
    try:
        if args.batch:
            if len(args.input) != 1 or not os.path.isdir(args.input[0]):
                logger.error("Modo lote requer exatamente um diretório de entrada")
                sys.exit(1)
            
            wrapper.batch_process(
                args.input[0], 
                args.output,
                args.pattern,
                args.threads
            )
        else:
            wrapper.run_fastqc(
                args.input,
                args.output,
                args.threads,
                args.quiet,
                not args.no_extract,
                args.format
            )
            
    except Exception as e:
        logger.error(f"Erro durante execução: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
