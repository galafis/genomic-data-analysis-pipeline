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


def _json_number(value):
    return None if pd.isna(value) else float(value)


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
        """Read the eight fixed VCF columns; report malformed records by line."""
        import math

        columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "DP"]
        variants = []
        header_seen = False
        with open(self.vcf_file, encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("##") or not line.strip():
                    continue
                if line.startswith("#CHROM"):
                    if (
                        line.rstrip().split("\t")[:8]
                        != ["#" + columns[0]] + columns[1:8]
                    ):
                        raise ValueError(
                            f"Invalid VCF header at line {line_number} / cabeçalho inválido"
                        )
                    header_seen = True
                    continue
                if line.startswith("#"):
                    continue
                try:
                    if not header_seen:
                        raise ValueError("missing #CHROM header")
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) < 8:
                        raise ValueError("expected eight columns")
                    chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
                    position = int(pos)
                    quality = float(qual) if qual != "." else None
                    depth = None
                    for item in info.split(";"):
                        if item.startswith("DP="):
                            raw = item[3:]
                            depth = None if raw == "." else int(raw)
                            break
                    if not chrom or not ref or not alt or position < 1:
                        raise ValueError("invalid locus")
                    if quality is not None and (
                        not math.isfinite(quality) or quality < 0
                    ):
                        raise ValueError("invalid QUAL")
                    if depth is not None and depth < 0:
                        raise ValueError("invalid DP")
                    variants.append(
                        [chrom, position, var_id, ref, alt, quality, filt, info, depth]
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid VCF record at line {line_number} / registro inválido: {exc}"
                    ) from exc
        if not header_seen:
            raise ValueError("Missing #CHROM header / cabeçalho ausente")
        result = pd.DataFrame(variants, columns=columns)
        for column in ("POS", "QUAL", "DP"):
            result[column] = pd.to_numeric(result[column])
        return result

    def filter_variants(
        self,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None,
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
        import math

        for name, value in (("min_qual", min_qual), ("min_dp", min_dp)):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError(
                    f"{name} must be finite and nonnegative / deve ser finito e não negativo"
                )
        if min_dp is not None and int(min_dp) != min_dp:
            raise ValueError("min_dp must be an integer / deve ser inteiro")
        df = self.variants.copy()

        if min_qual is not None:
            df = df[df["QUAL"] >= min_qual]

        if min_dp is not None:
            df = df[df["DP"] >= min_dp]

        if chromosomes is not None:
            df = df[df["CHROM"].isin(chromosomes)]

        return df

    def export_to_csv(
        self,
        output_file: str,
        min_qual: Optional[float] = None,
        min_dp: Optional[int] = None,
        chromosomes: Optional[List[str]] = None,
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
        sheet_name: str = "Variants",
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
        orient: str = "records",
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
            "total_variants": len(self.variants),
            "chromosomes": self.variants["CHROM"].unique().tolist(),
            "mean_qual": _json_number(self.variants["QUAL"].mean()),
            "mean_dp": _json_number(self.variants["DP"].mean()),
            "qual_range": (
                _json_number(self.variants["QUAL"].min()),
                _json_number(self.variants["QUAL"].max()),
            ),
            "dp_range": (
                _json_number(self.variants["DP"].min()),
                _json_number(self.variants["DP"].max()),
            ),
        }


def batch_export(
    vcf_file: str,
    output_dir: str,
    base_name: str,
    formats: Optional[List[str]] = None,
    min_qual: Optional[float] = None,
    min_dp: Optional[int] = None,
    chromosomes: Optional[List[str]] = None,
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
    import re

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", base_name):
        raise ValueError("base_name must be a plain filename / nome simples de arquivo")
    formats = ["csv", "excel", "json"] if formats is None else formats
    if not formats or not set(formats) <= {"csv", "excel", "json"}:
        raise ValueError("Unsupported export format / formato não suportado")
    exporter = VCFExporter(vcf_file)
    exporter.filter_variants(min_qual, min_dp, chromosomes)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    output_files = {}

    if "csv" in formats:
        csv_file = output_path / f"{base_name}.csv"
        output_files["csv"] = exporter.export_to_csv(
            str(csv_file), min_qual, min_dp, chromosomes
        )

    if "excel" in formats:
        excel_file = output_path / f"{base_name}.xlsx"
        output_files["excel"] = exporter.export_to_excel(
            str(excel_file), min_qual, min_dp, chromosomes
        )

    if "json" in formats:
        json_file = output_path / f"{base_name}.json"
        output_files["json"] = exporter.export_to_json(
            str(json_file), min_qual, min_dp, chromosomes
        )

    return output_files


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Filter and export VCF / Filtrar e exportar VCF"
    )
    parser.add_argument("vcf")
    parser.add_argument("--output-dir", default="output")
    parser.add_argument("--name", default="variants")
    parser.add_argument("--min-qual", type=float)
    parser.add_argument("--min-dp", type=int)
    args = parser.parse_args()
    batch_export(
        args.vcf,
        args.output_dir,
        args.name,
        ["csv", "json"],
        args.min_qual,
        args.min_dp,
    )
