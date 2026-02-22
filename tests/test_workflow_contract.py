from __future__ import annotations

from pathlib import Path


def test_snakefile_has_required_rules() -> None:
    text = Path("workflow/rules/pipeline.smk").read_text(encoding="utf-8")
    required = [
        "rule validate_inputs",
        "rule ingest_eqtl_colon_local",
        "rule ingest_eqtl_brain",
        "rule parse_eqtl_to_parquet",
        "rule harmonize_variants",
        "rule filter_significant_eqtl",
        "rule shared_specific_colon",
        "rule heterogeneity_colon",
        "rule gwas_harmonization",
        "rule coloc_colon_and_brain",
        "rule genetic_correlation",
        "rule bidirectional_mr",
        "rule pathway_enrichment",
        "rule scrna_support_colon",
        "rule assemble_results",
    ]
    for item in required:
        assert item in text
