import json
import pytest
from src.visualization.interactive.vcf_export_tools import VCFExporter, batch_export

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def make(tmp_path, row=""):
    path = tmp_path / "sample.vcf"
    path.write_text(HEADER + row)
    return path


def test_header_only_has_stable_schema(tmp_path):
    exporter = VCFExporter(make(tmp_path))
    assert exporter.filter_variants(min_dp=10).empty
    stats = exporter.get_summary_stats()
    assert stats["total_variants"] == 0 and stats["mean_dp"] is None
    json.dumps(stats, allow_nan=False)


def test_missing_depth(tmp_path):
    exporter = VCFExporter(make(tmp_path, "chr1\t1\t.\tA\tG\t.\tPASS\tDP=.\n"))
    assert exporter.filter_variants(min_dp=1).empty
    assert exporter.get_summary_stats()["mean_dp"] is None


@pytest.mark.parametrize(
    "row",
    [
        "bad\n",
        "chr1\t0\t.\tA\tG\t20\tPASS\tDP=2\n",
        "chr1\t1\t.\tA\tG\tNaN\tPASS\tDP=2\n",
    ],
)
def test_malformed_line_reports_number(tmp_path, row):
    with pytest.raises(ValueError, match="line 3"):
        VCFExporter(make(tmp_path, row))


def test_batch_rejects_path_escape_before_writing(tmp_path):
    with pytest.raises(ValueError):
        batch_export(str(make(tmp_path)), str(tmp_path / "out"), "../escape")
    assert not (tmp_path / "out").exists()


def test_unknown_format(tmp_path):
    with pytest.raises(ValueError):
        batch_export(str(make(tmp_path)), str(tmp_path / "out"), "sample", ["pdf"])
