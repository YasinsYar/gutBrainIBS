ACTIVE_CONFIG = config.get("config_file", "config/config.yaml")
GWAS_MANIFEST = config.get("gwas", {}).get("manifest", "config/gwas_manifest.yaml")

rule all:
    input:
        "data_processed/tables/shared_specific_colon.tsv",
        "data_processed/tables/heterogeneity_colon.tsv",
        "data_processed/tables/coloc_colon.tsv",
        "data_processed/tables/coloc_brain.tsv",
        "data_processed/tables/genetic_correlation.tsv",
        "data_processed/tables/mr_results.tsv",
        "data_processed/tables/pathway_enrichment.tsv",
        "data_processed/tables/scrna_support.tsv",
        "figures/colon_shared_specific_barplot.png",
        "figures/beta_concordance_scatter.png",
        "figures/heterogeneity_volcano.png",
        "figures/pp4_distribution.png",
        "figures/genetic_correlation_panel.png",
        "figures/mr_forest.png",
        "figures/pathway_enrichment_dotplot.png",
        "reports/run_qc.html",
        "data_processed/final/hypothesis_scorecard.json",
        "data_processed/final/evidence_summary.tsv",
        "data_processed/sensitivity/coloc_sensitivity_summary.tsv",
        "data_processed/sensitivity/coloc_sensitivity_brain_detail.tsv",
        "data_processed/sensitivity/mr_robustness.tsv",
        "data_processed/sensitivity/harmonization_strict_vs_nonstrict.tsv",
        "data_processed/sensitivity/coloc_brain_nonstrict.tsv",
        "data_processed/sensitivity/sensitivity_summary.tsv",
        "figures/presentation/01_colon_overlap_upset.png",
        "figures/presentation/02_colon_overlap_venn.png",
        "figures/presentation/03_beta_concordance_hexbin.png",
        "figures/presentation/04_heterogeneity_volcano_annotated.png",
        "figures/presentation/05_brain_coloc_heatmap.png",
        "figures/presentation/06_brain_coloc_locus_bubble.png",
        "figures/presentation/07_coloc_sensitivity_heatmap.png",
        "figures/presentation/08_mhc_sensitivity_bar.png",
        "figures/presentation/09_genetic_correlation_lollipop.png",
        "figures/presentation/10_mr_forest_robustness.png",
        "data_processed/tables/presentation/top_brain_coloc_hits.tsv",
        "data_processed/tables/presentation/mr_presentation.tsv",
        "data_processed/tables/presentation/sensitivity_key_metrics.tsv",
        "data_processed/tables/presentation/presentation_manifest.json",
        "reports/presentation_tables.html",
        "reports/presentation_tables.md"


rule stage_preprocess:
    input:
        "data_processed/qc/input_validation.json",
        "data_raw/qtl/colon_manifest.tsv",
        "data_raw/qtl/brain_manifest.tsv",
        "data_interim/eqtl_parsed.parquet",
        "data_interim/eqtl_harmonized.parquet",
        "data_interim/eqtl_significant.parquet",
        "data_interim/gwas_harmonized.parquet"


rule stage_association:
    input:
        "data_processed/tables/shared_specific_colon.tsv",
        "data_processed/tables/heterogeneity_colon.tsv",
        "data_processed/tables/coloc_colon.tsv",
        "data_processed/tables/coloc_brain.tsv",
        "data_processed/tables/genetic_correlation.tsv",
        "data_processed/tables/mr_results.tsv",
        "data_processed/tables/pathway_enrichment.tsv",
        "data_processed/tables/scrna_support.tsv"


rule stage_reporting:
    input:
        "figures/colon_shared_specific_barplot.png",
        "figures/beta_concordance_scatter.png",
        "figures/heterogeneity_volcano.png",
        "figures/pp4_distribution.png",
        "figures/genetic_correlation_panel.png",
        "figures/mr_forest.png",
        "figures/pathway_enrichment_dotplot.png",
        "reports/run_qc.html",
        "data_processed/final/hypothesis_scorecard.json",
        "data_processed/final/evidence_summary.tsv",
        "data_processed/sensitivity/sensitivity_summary.tsv",
        "reports/presentation_tables.html"


rule validate_inputs:
    output:
        inventory="data_processed/qc/input_inventory.tsv",
        validation_json="data_processed/qc/input_validation.json"
    shell:
        "python -m src.pipeline.validate_inputs --config {ACTIVE_CONFIG} --inventory {output.inventory} --summary-json {output.validation_json}"


rule ingest_eqtl_colon_local:
    input:
        "data_processed/qc/input_validation.json"
    output:
        sigmoid="data_raw/qtl/colon_sigmoid.all.tsv.gz",
        transverse="data_raw/qtl/colon_transverse.all.tsv.gz",
        manifest="data_raw/qtl/colon_manifest.tsv"
    shell:
        "python -m src.pipeline.ingest_eqtl_colon_local --config {ACTIVE_CONFIG} --out-sigmoid {output.sigmoid} --out-transverse {output.transverse} --manifest {output.manifest}"


rule ingest_eqtl_brain:
    input:
        "data_processed/qc/input_validation.json"
    output:
        manifest="data_raw/qtl/brain_manifest.tsv"
    shell:
        "python -m src.pipeline.ingest_eqtl_brain --config {ACTIVE_CONFIG} --manifest {output.manifest}"


rule parse_eqtl_to_parquet:
    input:
        colon_manifest="data_raw/qtl/colon_manifest.tsv",
        brain_manifest="data_raw/qtl/brain_manifest.tsv"
    output:
        parquet="data_interim/eqtl_parsed.parquet",
        parse_manifest="data_interim/eqtl_parse_manifest.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.parse_eqtl_to_parquet --config {ACTIVE_CONFIG} --colon-manifest {input.colon_manifest} --brain-manifest {input.brain_manifest} --out-parquet {output.parquet} --out-manifest {output.parse_manifest}"


rule harmonize_variants:
    input:
        "data_interim/eqtl_parsed.parquet"
    output:
        parquet="data_interim/eqtl_harmonized.parquet",
        qc="data_processed/qc/harmonization_qc.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.harmonize_variants --config {ACTIVE_CONFIG} --in-parquet {input} --out-parquet {output.parquet} --qc {output.qc}"


rule filter_significant_eqtl:
    input:
        "data_interim/eqtl_harmonized.parquet"
    output:
        "data_interim/eqtl_significant.parquet"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.filter_significant_eqtl --config {ACTIVE_CONFIG} --in-parquet {input} --out-parquet {output}"


rule shared_specific_colon:
    input:
        "data_interim/eqtl_significant.parquet"
    output:
        table="data_processed/tables/shared_specific_colon.tsv",
        metrics="data_processed/qc/shared_specific_metrics.json"
    shell:
        "python -m src.pipeline.shared_specific_colon --config {ACTIVE_CONFIG} --in-parquet {input} --out-table {output.table} --out-metrics {output.metrics}"


rule heterogeneity_colon:
    input:
        "data_interim/eqtl_significant.parquet"
    output:
        "data_processed/tables/heterogeneity_colon.tsv"
    shell:
        "python -m src.pipeline.heterogeneity_colon --config {ACTIVE_CONFIG} --in-parquet {input} --out-table {output}"


rule gwas_harmonization:
    input:
        "data_processed/qc/input_validation.json"
    output:
        parquet="data_interim/gwas_harmonized.parquet",
        qc="data_processed/qc/gwas_harmonization_qc.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.gwas_harmonization --config {ACTIVE_CONFIG} --manifest {GWAS_MANIFEST} --out-parquet {output.parquet} --qc {output.qc}"


rule coloc_colon_and_brain:
    input:
        eqtl="data_interim/eqtl_significant.parquet",
        gwas="data_interim/gwas_harmonized.parquet"
    output:
        colon="data_processed/tables/coloc_colon.tsv",
        brain="data_processed/tables/coloc_brain.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.coloc_colon_and_brain --config {ACTIVE_CONFIG} --eqtl {input.eqtl} --gwas {input.gwas} --out-colon {output.colon} --out-brain {output.brain}"


rule genetic_correlation:
    input:
        "data_interim/gwas_harmonized.parquet"
    output:
        "data_processed/tables/genetic_correlation.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.genetic_correlation --config {ACTIVE_CONFIG} --in-parquet {input} --out-table {output}"


rule bidirectional_mr:
    input:
        "data_interim/gwas_harmonized.parquet"
    output:
        "data_processed/tables/mr_results.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.bidirectional_mr --config {ACTIVE_CONFIG} --in-parquet {input} --out-table {output}"


rule pathway_enrichment:
    input:
        coloc_colon="data_processed/tables/coloc_colon.tsv",
        coloc_brain="data_processed/tables/coloc_brain.tsv",
        mr="data_processed/tables/mr_results.tsv"
    output:
        "data_processed/tables/pathway_enrichment.tsv"
    shell:
        "python -m src.pipeline.pathway_enrichment --config {ACTIVE_CONFIG} --coloc-colon {input.coloc_colon} --coloc-brain {input.coloc_brain} --mr {input.mr} --out-table {output}"


rule scrna_support_colon:
    input:
        "data_processed/tables/coloc_colon.tsv"
    output:
        "data_processed/tables/scrna_support.tsv"
    shell:
        "python -m src.pipeline.scrna_support_colon --config {ACTIVE_CONFIG} --coloc-table {input} --out-table {output}"


rule build_figures:
    input:
        shared="data_processed/tables/shared_specific_colon.tsv",
        het="data_processed/tables/heterogeneity_colon.tsv",
        coloc_colon="data_processed/tables/coloc_colon.tsv",
        coloc_brain="data_processed/tables/coloc_brain.tsv",
        rg="data_processed/tables/genetic_correlation.tsv",
        mr="data_processed/tables/mr_results.tsv",
        pathway="data_processed/tables/pathway_enrichment.tsv"
    output:
        "figures/colon_shared_specific_barplot.png",
        "figures/beta_concordance_scatter.png",
        "figures/heterogeneity_volcano.png",
        "figures/pp4_distribution.png",
        "figures/genetic_correlation_panel.png",
        "figures/mr_forest.png",
        "figures/pathway_enrichment_dotplot.png"
    shell:
        "python -m src.pipeline.build_figures --config {ACTIVE_CONFIG} --shared {input.shared} --het {input.het} --coloc-colon {input.coloc_colon} --coloc-brain {input.coloc_brain} --rg {input.rg} --mr {input.mr} --pathway {input.pathway} --fig-dir figures"


rule build_presentation_assets:
    input:
        shared="data_processed/tables/shared_specific_colon.tsv",
        het="data_processed/tables/heterogeneity_colon.tsv",
        coloc_brain="data_processed/tables/coloc_brain.tsv",
        rg="data_processed/tables/genetic_correlation.tsv",
        mr="data_processed/tables/mr_results.tsv",
        coloc_sens="data_processed/sensitivity/coloc_sensitivity_summary.tsv",
        mr_rob="data_processed/sensitivity/mr_robustness.tsv",
        harm_sens="data_processed/sensitivity/harmonization_strict_vs_nonstrict.tsv"
    output:
        "figures/presentation/01_colon_overlap_upset.png",
        "figures/presentation/02_colon_overlap_venn.png",
        "figures/presentation/03_beta_concordance_hexbin.png",
        "figures/presentation/04_heterogeneity_volcano_annotated.png",
        "figures/presentation/05_brain_coloc_heatmap.png",
        "figures/presentation/06_brain_coloc_locus_bubble.png",
        "figures/presentation/07_coloc_sensitivity_heatmap.png",
        "figures/presentation/08_mhc_sensitivity_bar.png",
        "figures/presentation/09_genetic_correlation_lollipop.png",
        "figures/presentation/10_mr_forest_robustness.png",
        "data_processed/tables/presentation/top_brain_coloc_hits.tsv",
        "data_processed/tables/presentation/mr_presentation.tsv",
        "data_processed/tables/presentation/sensitivity_key_metrics.tsv",
        "data_processed/tables/presentation/presentation_manifest.json",
        "reports/presentation_tables.html",
        "reports/presentation_tables.md"
    shell:
        "python -m src.pipeline.build_presentation_assets --shared {input.shared} --het {input.het} --coloc-brain {input.coloc_brain} --rg {input.rg} --mr {input.mr} --coloc-sensitivity {input.coloc_sens} --mr-robustness {input.mr_rob} --harmonization-sensitivity {input.harm_sens} --fig-dir figures/presentation --table-dir data_processed/tables/presentation --out-html reports/presentation_tables.html --out-md reports/presentation_tables.md"


rule assemble_results:
    input:
        validation="data_processed/qc/input_validation.json",
        shared_metrics="data_processed/qc/shared_specific_metrics.json",
        rg="data_processed/tables/genetic_correlation.tsv",
        mr="data_processed/tables/mr_results.tsv",
        coloc_brain="data_processed/tables/coloc_brain.tsv",
        figures_done="figures/pathway_enrichment_dotplot.png"
    output:
        scorecard="data_processed/final/hypothesis_scorecard.json",
        evidence="data_processed/final/evidence_summary.tsv",
        run_qc="reports/run_qc.html"
    shell:
        "python -m src.pipeline.assemble_results --config {ACTIVE_CONFIG} --validation-json {input.validation} --shared-metrics {input.shared_metrics} --rg {input.rg} --mr {input.mr} --coloc-brain {input.coloc_brain} --out-scorecard {output.scorecard} --out-evidence {output.evidence} --out-qc-html {output.run_qc}"


rule sensitivity_all:
    input:
        "data_processed/sensitivity/sensitivity_summary.tsv"


rule coloc_sensitivity:
    input:
        eqtl="data_interim/eqtl_significant.parquet",
        gwas="data_interim/gwas_harmonized.parquet"
    output:
        summary="data_processed/sensitivity/coloc_sensitivity_summary.tsv",
        brain_detail="data_processed/sensitivity/coloc_sensitivity_brain_detail.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.coloc_sensitivity --config {ACTIVE_CONFIG} --eqtl {input.eqtl} --gwas {input.gwas} --out-summary {output.summary} --out-brain-detail {output.brain_detail}"


rule mr_sensitivity:
    input:
        "data_interim/gwas_harmonized.parquet"
    output:
        "data_processed/sensitivity/mr_robustness.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.mr_sensitivity --config {ACTIVE_CONFIG} --in-parquet {input} --out-table {output}"


rule harmonization_sensitivity:
    input:
        eqtl_parsed="data_interim/eqtl_parsed.parquet",
        eqtl_significant="data_interim/eqtl_significant.parquet",
        gwas="data_interim/gwas_harmonized.parquet",
        strict_coloc_brain="data_processed/tables/coloc_brain.tsv",
        strict_scorecard="data_processed/final/hypothesis_scorecard.json"
    output:
        summary="data_processed/sensitivity/harmonization_strict_vs_nonstrict.tsv",
        nonstrict_coloc_brain="data_processed/sensitivity/coloc_brain_nonstrict.tsv"
    resources:
        io_slots=1
    shell:
        "python -m src.pipeline.harmonization_sensitivity --config {ACTIVE_CONFIG} --eqtl-parsed {input.eqtl_parsed} --eqtl-significant {input.eqtl_significant} --gwas {input.gwas} --strict-coloc-brain {input.strict_coloc_brain} --strict-scorecard {input.strict_scorecard} --out-table {output.summary} --out-coloc-brain-nonstrict {output.nonstrict_coloc_brain}"


rule run_sensitivity:
    input:
        coloc_summary="data_processed/sensitivity/coloc_sensitivity_summary.tsv",
        mr_robustness="data_processed/sensitivity/mr_robustness.tsv",
        harmonization_compare="data_processed/sensitivity/harmonization_strict_vs_nonstrict.tsv"
    output:
        "data_processed/sensitivity/sensitivity_summary.tsv"
    shell:
        "python -m src.pipeline.run_sensitivity --coloc-summary {input.coloc_summary} --mr-robustness {input.mr_robustness} --harmonization-compare {input.harmonization_compare} --out-table {output}"
