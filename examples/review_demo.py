from pathlib import Path
from src.visualization.interactive.vcf_export_tools import VCFExporter


def run():
    exporter = VCFExporter(Path(__file__).with_name("sample.vcf"))
    selected = exporter.filter_variants(min_qual=30, min_dp=10)
    assert list(selected.ID) == ["demo-1", "demo-4"]
    return {
        "synthetic": True,
        "filters": {"min_qual": 30, "min_dp": 10},
        "selected_ids": list(selected.ID),
        "summary": exporter.get_summary_stats(),
    }


if __name__ == "__main__":
    import json

    print(json.dumps(run(), ensure_ascii=False, indent=2, allow_nan=False))
