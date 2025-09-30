import json
import os
from pathlib import Path

import pandas as pd
import pytest

# Import the module under test
from src.visualization.interactive.vcf_export_tools import VCFExporter, batch_export


@pytest.fixture()
def tmp_vcf_file(tmp_path: Path) -> Path:
    """Create a minimal VCF file for tests and return its path.

    The VCF contains 5 variants across 2 chromosomes with assorted QUAL and DP values.
    """
    vcf_content = """##fileformat=VCFv4.2
##source=pytest
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs1\tA\tG\t50\tPASS\tDP=12;OTHER=foo
chr1\t150\trs2\tC\tT\t20\tPASS\tDP=8
chr1\t200\trs3\tG\tA\t.\tLowQual\tDP=30
chr2\t100\trs4\tT\tC\t60\tPASS\tDP=25
chr2\t300\trs5\tA\tAT\t10\tPASS\tDP=3
"""
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text(vcf_content)
    return vcf_path


@pytest.fixture()
def exporter(tmp_vcf_file: Path) -> VCFExporter:
    return VCFExporter(str(tmp_vcf_file))


def test_parse_vcf_loads_dataframe(exporter: VCFExporter):
    df = exporter.variants
    # basic shape and columns
    assert isinstance(df, pd.DataFrame)
    assert set(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "DP"]).issubset(df.columns)
    assert len(df) == 5
    # types/parsing
    assert df.loc[0, "CHROM"] == "chr1"
    assert isinstance(df.loc[0, "POS"], (int,))
    # QUAL '.' becomes None
    assert pd.isna(df.loc[2, "QUAL"]) or df.loc[2, "QUAL"] is None
    # DP parsed from INFO
    assert df.loc[0, "DP"] == 12


def test_filter_variants_by_min_qual(exporter: VCFExporter):
    filtered = exporter.filter_variants(min_qual=30)
    # rs2 has QUAL 20, rs5 has 10, rs3 is None -> all should be filtered out
    assert set(filtered["ID"]).issuperset({"rs1", "rs4"})
    assert {"rs2", "rs5"}.isdisjoint(set(filtered["ID"]))


def test_filter_variants_by_min_dp(exporter: VCFExporter):
    filtered = exporter.filter_variants(min_dp=10)
    # rs5 has DP=3 should be removed
    assert "rs5" not in set(filtered["ID"])   
    # others remain
    assert {"rs1", "rs2", "rs3", "rs4"}.issubset(set(filtered["ID"]))


def test_filter_variants_by_chromosomes(exporter: VCFExporter):
    filtered = exporter.filter_variants(chromosomes=["chr2"]) 
    assert set(filtered["CHROM"]) == {"chr2"}
    assert set(filtered["ID"]) == {"rs4", "rs5"}


def test_export_to_csv_creates_file_and_matches_dataframe(exporter: VCFExporter, tmp_path: Path):
    out_csv = tmp_path / "variants.csv"
    path = exporter.export_to_csv(str(out_csv), min_qual=15)
    assert path == str(out_csv)
    assert out_csv.exists()

    df_expected = exporter.filter_variants(min_qual=15)
    df_written = pd.read_csv(out_csv)
    # Compare core columns and lengths
    assert len(df_written) == len(df_expected)
    assert list(df_written.columns) == list(df_expected.columns)


def test_export_to_excel_creates_file_and_sheet(exporter: VCFExporter, tmp_path: Path):
    out_xlsx = tmp_path / "variants.xlsx"
    path = exporter.export_to_excel(str(out_xlsx), min_dp=10, sheet_name="Variants")
    assert path == str(out_xlsx)
    assert out_xlsx.exists()

    # Read back and check shape/columns
    df_expected = exporter.filter_variants(min_dp=10)
    df_written = pd.read_excel(out_xlsx, sheet_name="Variants")
    assert len(df_written) == len(df_expected)
    assert list(df_written.columns) == list(df_expected.columns)


def test_export_to_json_creates_file_and_content(exporter: VCFExporter, tmp_path: Path):
    out_json = tmp_path / "variants.json"
    path = exporter.export_to_json(str(out_json), chromosomes=["chr1"], orient="records")
    assert path == str(out_json)
    assert out_json.exists()

    df_expected = exporter.filter_variants(chromosomes=["chr1"]).reset_index(drop=True)
    data = json.loads(out_json.read_text())
    assert isinstance(data, list)
    assert len(data) == len(df_expected)
    # Basic integrity: IDs and CHROM match
    ids = {rec.get("ID") for rec in data}
    assert ids == set(df_expected["ID"])


def test_batch_export_multiple_formats(tmp_vcf_file: Path, tmp_path: Path):
    outdir = tmp_path / "output"
    files = batch_export(
        vcf_file=str(tmp_vcf_file),
        output_dir=str(outdir),
        base_name="filtered_variants",
        formats=["csv", "excel", "json"],
        min_qual=15,
        min_dp=5,
        chromosomes=["chr1", "chr2"],
    )

    # Returned mapping and files exist
    assert set(files.keys()) == {"csv", "excel", "json"}
    for fmt, fpath in files.items():
        assert Path(fpath).exists(), f"Missing exported file for {fmt}: {fpath}"


def test_get_summary_stats(exporter: VCFExporter):
    stats = exporter.get_summary_stats()
    assert "total_variants" in stats and stats["total_variants"] == 5
    assert set(stats["chromosomes"]) == {"chr1", "chr2"}
    # Means and ranges should be computable given numeric/None QUAL and integer DP
    assert "mean_dp" in stats and stats["mean_dp"] > 0
    assert "dp_range" in stats and isinstance(stats["dp_range"], tuple)
