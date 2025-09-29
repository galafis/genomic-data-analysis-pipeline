#!/usr/bin/env python3
"""
VCF Export Tools

This module provides functions to export VCF variants to various formats (CSV, Excel, JSON)
with optional filtering by QUAL, DP (depth), and chromosome.

Author: Gabriel
Date: September 29, 2025
"""

import pandas as pd
import json
from typing import Optional, List, Dict, Any
from pathlib import Path


class VCFExporter:
    """
    A class to handle VCF file parsing and exporting to multiple formats.
    
    Attributes:
        vcf_file (str): Path to the input VCF file
        variants (pd.DataFrame): DataFrame containing parsed VCF variants
    """
    
    def __init__(self, vcf_file: str):
        """
        Initialize the VCFExporter.
        
        Args:
            vcf_file (str): Path to the VCF file to be processed
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
        """
        self.vcf_file = vcf_file
        self.variants = self._parse_vcf()
    
    def _parse_vcf(self) -> pd.DataFrame:
        """
        Parse VCF file into a pandas DataFrame.
        
        Returns:
            pd.DataFrame: DataFrame with variant information
        """
        variants = []
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
                
                # Parse INFO field for DP
                dp = None
                for item in info.split(';'):
                    if item.startswith('DP='):
                        dp = int(item.split('=')[1])
                        break
                
                variant = {
                    'CHROM': chrom,
                    'POS': int(pos),
                    'ID': var_id,
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': float(qual) if qual != '.' else None,
                    'FILTER': filt,
                    'INFO': info,
                    'DP': dp
                }
                variants.append(variant)
        
        return pd.DataFrame(variants)
    
    def filter_variants(
        self,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Filter variants based on quality, depth, and chromosome.
        
        Args:
            min_qual (float, optional): Minimum QUAL score threshold
            min_dp (int, optional): Minimum depth (DP) threshold
            chromosomes (List[str], optional): List of chromosomes to include (e.g., ['chr1', 'chr2'])
            
        Returns:
            pd.DataFrame: Filtered DataFrame
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
            >>> filtered = exporter.filter_variants(min_qual=30, min_dp=10, chromosomes=['chr1'])
        """
        df = self.variants.copy()
        
        if min_qual is not None:
            df = df[df['QUAL'] >= min_qual]
        
        if min_dp is not None:
            df = df[df['DP'] >= min_dp]
        
        if chromosomes is not None:
            df = df[df['CHROM'].isin(chromosomes)]
        
        return df
    
    def export_to_csv(
        self,
        output_file: str,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None
    ) -> str:
        """
        Export variants to CSV format with optional filtering.
        
        Args:
            output_file (str): Path for the output CSV file
            min_qual (float, optional): Minimum QUAL score threshold
            min_dp (int, optional): Minimum depth (DP) threshold
            chromosomes (List[str], optional): List of chromosomes to include
            
        Returns:
            str: Path to the created CSV file
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
            >>> exporter.export_to_csv('variants.csv', min_qual=30, min_dp=10)
            'variants.csv'
        """
        df = self.filter_variants(min_qual, min_dp, chromosomes)
        df.to_csv(output_file, index=False)
        print(f"Exported {len(df)} variants to {output_file}")
        return output_file
    
    def export_to_excel(
        self,
        output_file: str,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None,
        sheet_name: str = 'Variants'
    ) -> str:
        """
        Export variants to Excel format with optional filtering.
        
        Args:
            output_file (str): Path for the output Excel file
            min_qual (float, optional): Minimum QUAL score threshold
            min_dp (int, optional): Minimum depth (DP) threshold
            chromosomes (List[str], optional): List of chromosomes to include
            sheet_name (str): Name of the Excel sheet (default: 'Variants')
            
        Returns:
            str: Path to the created Excel file
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
            >>> exporter.export_to_excel('variants.xlsx', min_qual=30, chromosomes=['chr1', 'chr2'])
            'variants.xlsx'
        """
        df = self.filter_variants(min_qual, min_dp, chromosomes)
        df.to_excel(output_file, index=False, sheet_name=sheet_name)
        print(f"Exported {len(df)} variants to {output_file}")
        return output_file
    
    def export_to_json(
        self,
        output_file: str,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None,
        orient: str = 'records'
    ) -> str:
        """
        Export variants to JSON format with optional filtering.
        
        Args:
            output_file (str): Path for the output JSON file
            min_qual (float, optional): Minimum QUAL score threshold
            min_dp (int, optional): Minimum depth (DP) threshold
            chromosomes (List[str], optional): List of chromosomes to include
            orient (str): JSON orientation ('records', 'index', 'columns', etc.)
            
        Returns:
            str: Path to the created JSON file
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
            >>> exporter.export_to_json('variants.json', min_dp=20, orient='records')
            'variants.json'
        """
        df = self.filter_variants(min_qual, min_dp, chromosomes)
        df.to_json(output_file, orient=orient, indent=2)
        print(f"Exported {len(df)} variants to {output_file}")
        return output_file
    
    def get_summary_stats(self) -> Dict[str, Any]:
        """
        Get summary statistics for the VCF data.
        
        Returns:
            dict: Dictionary containing summary statistics
            
        Example:
            >>> exporter = VCFExporter('sample.vcf')
            >>> stats = exporter.get_summary_stats()
            >>> print(f"Total variants: {stats['total_variants']}")
        """
        return {
            'total_variants': len(self.variants),
            'chromosomes': self.variants['CHROM'].unique().tolist(),
            'mean_qual': self.variants['QUAL'].mean(),
            'mean_dp': self.variants['DP'].mean(),
            'qual_range': (self.variants['QUAL'].min(), self.variants['QUAL'].max()),
            'dp_range': (self.variants['DP'].min(), self.variants['DP'].max())
        }


def batch_export(
    vcf_file: str,
    output_dir: str,
    base_name: str,
    formats: List[str] = ['csv', 'excel', 'json'],
    min_qual: Optional[float] = None,
    min_dp: Optional[int] = None,
    chromosomes: Optional[List[str]] = None
) -> Dict[str, str]:
    """
    Export VCF variants to multiple formats at once.
    
    Args:
        vcf_file (str): Path to the input VCF file
        output_dir (str): Directory for output files
        base_name (str): Base name for output files (without extension)
        formats (List[str]): List of formats to export ('csv', 'excel', 'json')
        min_qual (float, optional): Minimum QUAL score threshold
        min_dp (int, optional): Minimum depth (DP) threshold
        chromosomes (List[str], optional): List of chromosomes to include
        
    Returns:
        dict: Dictionary mapping format to output file path
        
    Example:
        >>> files = batch_export(
        ...     'sample.vcf',
        ...     'output',
        ...     'filtered_variants',
        ...     formats=['csv', 'json'],
        ...     min_qual=30,
        ...     min_dp=10
        ... )
        >>> print(files)
        {'csv': 'output/filtered_variants.csv', 'json': 'output/filtered_variants.json'}
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    exporter = VCFExporter(vcf_file)
    output_files = {}
    
    if 'csv' in formats:
        csv_file = output_path / f"{base_name}.csv"
        output_files['csv'] = exporter.export_to_csv(
            str(csv_file), min_qual, min_dp, chromosomes
        )
    
    if 'excel' in formats:
        excel_file = output_path / f"{base_name}.xlsx"
        output_files['excel'] = exporter.export_to_excel(
            str(excel_file), min_qual, min_dp, chromosomes
        )
    
    if 'json' in formats:
        json_file = output_path / f"{base_name}.json"
        output_files['json'] = exporter.export_to_json(
            str(json_file), min_qual, min_dp, chromosomes
        )
    
    return output_files


if __name__ == '__main__':
    """
    Example usage of the VCF export tools.
    """
    print("VCF Export Tools - Example Usage")
    print("=" * 50)
    
    # Example 1: Basic export
    print("\nExample 1: Basic CSV export")
    print("-" * 50)
    print("exporter = VCFExporter('sample.vcf')")
    print("exporter.export_to_csv('variants.csv')")
    
    # Example 2: Filtered export
    print("\nExample 2: Filtered Excel export")
    print("-" * 50)
    print("exporter = VCFExporter('sample.vcf')")
    print("exporter.export_to_excel(")
    print("    'high_quality_variants.xlsx',")
    print("    min_qual=30,")
    print("    min_dp=10,")
    print("    chromosomes=['chr1', 'chr2']")
    print(")")
    
    # Example 3: Batch export
    print("\nExample 3: Batch export to multiple formats")
    print("-" * 50)
    print("files = batch_export(")
    print("    'sample.vcf',")
    print("    'output',")
    print("    'filtered_variants',")
    print("    formats=['csv', 'excel', 'json'],")
    print("    min_qual=30,")
    print("    min_dp=10")
    print(")")
    print("print(files)")
    
    # Example 4: Get summary statistics
    print("\nExample 4: Get summary statistics")
    print("-" * 50)
    print("exporter = VCFExporter('sample.vcf')")
    print("stats = exporter.get_summary_stats()")
    print("print(f'Total variants: {stats[\"total_variants\"]}')") 
    print("print(f'Mean QUAL: {stats[\"mean_qual\"]:.2f}')")
    print("print(f'Mean DP: {stats[\"mean_dp\"]:.2f}')")
