#!/usr/bin/env python3
"""
Quality Control Module for NGS Data

This module provides comprehensive quality control functionality for
next-generation sequencing data including FastQC analysis, read statistics,
and quality metric aggregation.

Author: Gabriel Demetrios Lafis
Project: Genomic Data Analysis Pipeline
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Optional, Union
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class QualityController:
    """
    Main class for performing quality control on FASTQ files.
    
    This class provides methods for running FastQC, parsing results,
    and generating comprehensive QC reports.
    """
    
    def __init__(self, output_dir: str = "qc_results", threads: int = 4):
        """
        Initialize QualityController.
        
        Args:
            output_dir: Directory to store QC results
            threads: Number of threads to use for analysis
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        logger.info(f"QualityController initialized with output_dir={output_dir}")
    
    def run_fastqc(self, fastq_files: List[str], 
                   extract: bool = True) -> Dict[str, str]:
        """
        Run FastQC on list of FASTQ files.
        
        Args:
            fastq_files: List of FASTQ file paths
            extract: Whether to extract FastQC zip files
        
        Returns:
            Dictionary mapping file names to output paths
        """
        logger.info(f"Running FastQC on {len(fastq_files)} files")
        
        results = {}
        
        for fastq in fastq_files:
            if not Path(fastq).exists():
                logger.error(f"File not found: {fastq}")
                continue
            
            logger.info(f"Processing: {fastq}")
            
            cmd = [
                "fastqc",
                "-o", str(self.output_dir),
                "-t", str(self.threads),
                fastq
            ]
            
            if extract:
                cmd.append("--extract")
            
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )
                
                # Store output path
                base_name = Path(fastq).stem.replace(".fastq", "")
                results[fastq] = str(self.output_dir / f"{base_name}_fastqc")
                
                logger.info(f"FastQC completed for {fastq}")
                
            except subprocess.CalledProcessError as e:
                logger.error(f"FastQC failed for {fastq}: {e.stderr}")
                continue
        
        return results
    
    def parse_fastqc_data(self, fastqc_dir: str) -> Dict:
        """
        Parse FastQC data file and extract key metrics.
        
        Args:
            fastqc_dir: Path to extracted FastQC directory
        
        Returns:
            Dictionary containing parsed metrics
        """
        data_file = Path(fastqc_dir) / "fastqc_data.txt"
        
        if not data_file.exists():
            logger.error(f"FastQC data file not found: {data_file}")
            return {}
        
        metrics = {
            "filename": None,
            "total_sequences": 0,
            "sequence_length": None,
            "gc_content": 0,
            "modules": {}
        }
        
        with open(data_file, 'r') as f:
            current_module = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith("Filename"):
                    metrics["filename"] = line.split("\t")[1]
                elif line.startswith("Total Sequences"):
                    metrics["total_sequences"] = int(line.split("\t")[1])
                elif line.startswith("Sequence length"):
                    metrics["sequence_length"] = line.split("\t")[1]
                elif line.startswith("%GC"):
                    metrics["gc_content"] = float(line.split("\t")[1])
                elif line.startswith(">>") and line.endswith("END_MODULE"):
                    current_module = None
                elif line.startswith(">>"):
                    parts = line.split("\t")
                    module_name = parts[0].replace(">>", "").strip()
                    status = parts[1] if len(parts) > 1 else "unknown"
                    metrics["modules"][module_name] = status
                    current_module = module_name
        
        logger.info(f"Parsed FastQC data for {metrics.get('filename', 'unknown')}")
        return metrics
    
    def calculate_read_stats(self, fastq_file: str) -> Dict:
        """
        Calculate basic read statistics from FASTQ file.
        
        Args:
            fastq_file: Path to FASTQ file
        
        Returns:
            Dictionary with read statistics
        """
        logger.info(f"Calculating read stats for {fastq_file}")
        
        stats = {
            "total_reads": 0,
            "total_bases": 0,
            "read_lengths": []
        }
        
        try:
            with open(fastq_file, 'r') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    if line_num % 4 == 2:  # Sequence line
                        seq = line.strip()
                        stats["total_reads"] += 1
                        stats["total_bases"] += len(seq)
                        stats["read_lengths"].append(len(seq))
            
            if stats["read_lengths"]:
                stats["avg_read_length"] = sum(stats["read_lengths"]) / len(stats["read_lengths"])
                stats["min_read_length"] = min(stats["read_lengths"])
                stats["max_read_length"] = max(stats["read_lengths"])
            
            # Remove full list to save memory
            del stats["read_lengths"]
            
            logger.info(f"Read stats calculated: {stats['total_reads']} reads")
            
        except Exception as e:
            logger.error(f"Error calculating read stats: {e}")
        
        return stats
    
    def generate_summary_report(self, qc_results: Dict, 
                               output_file: str = "qc_summary.json") -> str:
        """
        Generate summary report from QC results.
        
        Args:
            qc_results: Dictionary of QC results
            output_file: Output file name
        
        Returns:
            Path to summary report
        """
        output_path = self.output_dir / output_file
        
        with open(output_path, 'w') as f:
            json.dump(qc_results, f, indent=2)
        
        logger.info(f"Summary report generated: {output_path}")
        return str(output_path)


def main():
    """
    Example usage of QualityController.
    """
    # Initialize controller
    qc = QualityController(output_dir="results/qc", threads=4)
    
    # Example FASTQ files (replace with actual paths)
    fastq_files = [
        "data/sample1_R1.fastq.gz",
        "data/sample1_R2.fastq.gz"
    ]
    
    # Run FastQC
    logger.info("Starting quality control analysis...")
    fastqc_results = qc.run_fastqc(fastq_files)
    
    # Parse results and generate report
    all_results = {}
    for fastq, result_dir in fastqc_results.items():
        metrics = qc.parse_fastqc_data(result_dir)
        all_results[fastq] = metrics
    
    # Generate summary
    summary_path = qc.generate_summary_report(all_results)
    logger.info(f"QC analysis complete. Summary: {summary_path}")


if __name__ == "__main__":
    main()
