# Materials and Methods

## Study Objective and Design
This study evaluated whether integrated human genetics supports gut-brain involvement in irritable bowel syndrome (IBS) within a disorders-of-gut-brain-interaction framework. The analysis combined colon and brain eQTL signals with IBS and psychiatric GWAS summary statistics, followed by harmonized association, colocalization, correlation, and Mendelian randomization (MR) layers in a single reproducible workflow.

The analytical design was observational and summary-statistics based, with predefined decision criteria and guardrails to avoid binary psychosomatic interpretation.

## Data Sources

### 1. Colon eQTL datasets (mandatory)
Local eQTL Catalogue gene-expression datasets were reused as the primary colon input:
- `QTD000226.all.tsv.gz` (colon sigmoid)
- `QTD000231.all.tsv.gz` (colon transverse)

### 2. Brain eQTL datasets (mandatory)
Brain gene-expression eQTL datasets were ingested from eQTL Catalogue (GTEx-derived, ge quantification), using manifest-based selection. In this run, four datasets were used:
- `QTD000146`
- `QTD000151`
- `QTD000156`
- `QTD000161`

### 3. GWAS summary statistics (mandatory)
Traits were configured through `config/gwas_manifest.yaml` and included:
- IBS (primary trait)
- Major depressive disorder (MDD)
- Anxiety
- Neuroticism

All GWAS were harmonized to a common schema and used in downstream coloc, genetic-correlation proxy analysis, and bidirectional MR.

### 4. Supportive single-cell context (non-causal layer)
Colon scRNA-seq files were used only for expression context of prioritized genes, not for causal inference.

## Computational Environment and Orchestration
- Operating system target: Windows (PowerShell execution)
- Workflow engine: Snakemake
- Primary language: Python
- Dependency management: Conda/Mamba (`envs/environment.yml`)
- Main entrypoint:
  - `snakemake --profile profiles/windows-default all`
- Managed execution (heartbeat/status/stop/cleanup controls):
  - `python -m src.pipeline.run_managed --profile profiles/windows-default --configfile config/config.yaml --target all`

## Software Stack
Core libraries:
- `pandas`, `numpy`, `pyarrow`, `duckdb`
- `scipy`, `pyyaml`, `requests`
- `matplotlib`, `seaborn`
- `anndata`, `h5py`

Sensitivity/presentation visualization:
- `upsetplot`
- `adjustText`
- `matplotlib-venn`
- `tabulate`

## Pipeline Structure and Steps
The workflow was implemented in `workflow/rules/pipeline.smk` and executed as a reproducible DAG:

1. `validate_inputs`
- Verifies required files, readability, and basic schema expectations.

2. `ingest_eqtl_colon_local`
- Stages local colon eQTL files into pipeline input structure with manifest logging.

3. `ingest_eqtl_brain`
- Retrieves/validates selected brain eQTL datasets from configured catalogue sources.

4. `parse_eqtl_to_parquet`
- Converts raw eQTL tables to normalized parquet representation for scalable processing.

5. `harmonize_variants`
- Standardizes allele representation and direction conventions.
- Strict mode excludes palindromic SNPs by default.

6. `filter_significant_eqtl`
- Applies significance threshold (`p < 1e-6`) per configuration.

7. `shared_specific_colon`
- Computes shared vs tissue-specific colon eQTL definitions across sigmoid/transverse tissues.

8. `heterogeneity_colon`
- Performs between-tissue heterogeneity testing on shared variant-gene pairs.

9. `gwas_harmonization`
- Harmonizes IBS and psychiatric GWAS into a unified schema.

10. `coloc_colon_and_brain`
- Performs ABF-based colocalization for colon and brain windows around IBS loci.

11. `genetic_correlation`
- Computes IBS-psychiatric proxy genetic-correlation estimates.

12. `bidirectional_mr`
- Runs psych->IBS and IBS->psych MR directions.
- Primary estimator: IVW; sensitivity estimators: weighted median and MR-Egger.

13. `pathway_enrichment`
- Performs enrichment analysis for prioritized genes if pathway database is provided.

14. `scrna_support_colon`
- Provides supportive cell-type expression context for selected genes.

15. `assemble_results`
- Produces final scorecard and evidence tables with interpretation labels.

16. Sensitivity layer
- `coloc_sensitivity`
- `mr_sensitivity`
- `harmonization_sensitivity`
- `run_sensitivity` (aggregation)

17. Presentation layer
- `build_figures` (core article figures)
- `build_presentation_assets` (enhanced visual panels and presentation tables)

## Statistical Methods

### Shared vs tissue-specific eQTL
Variant-gene pairs were classified as:
- shared: significant in both colon tissues
- tissue-specific: significant in only one tissue

Additional metrics included effect-sign concordance and absolute effect-size difference (`|delta beta|`).

### Between-tissue heterogeneity
For shared pairs, heterogeneity was tested with a Z-statistic based on tissue-specific effect and SE values, followed by multiple-testing adjustment.

### Colocalization
Colocalization used Wakefield approximate Bayes factors (ABF) with coloc posterior probabilities (`PP0..PP4`) under configured priors:
- `p1 = 1e-4`
- `p2 = 1e-4`
- `p12 = 1e-5`

Default locus window: `+/-500 kb` around IBS lead loci.

### Genetic-correlation layer
IBS-psychiatric overlap was estimated using a proxy Z-correlation implementation on harmonized summary statistics.

### Bidirectional MR
MR directions were evaluated using:
- IVW (primary)
- weighted median (sensitivity)
- MR-Egger slope/intercept (sensitivity)

Directional support rule in scorecard required: significant IVW result, directional consistency, and no pleiotropy flag.

## Sensitivity Analyses
Minimal high-yield sensitivity package included:

1. Brain-coloc MHC sensitivity
- baseline vs MHC-excluded vs chr6-only scenarios.

2. Coloc parameter sensitivity
- window perturbation (`+/-250 kb`, `+/-500 kb`, `+/-1 Mb`)
- prior perturbation (`p12 = 1e-6`, `1e-5`, `1e-4`)

3. MR robustness
- Cochran's Q heterogeneity
- leave-one-out influence checks
- Steiger proxy support flag

4. Harmonization sensitivity
- strict vs non-strict palindromic-SNP handling and impact assessment on coloc outputs.

## Quality Control and Reproducibility
- Config-driven run behavior (`config/config.yaml`)
- Deterministic workflow graph and tracked outputs
- Managed run status telemetry (`run_managed`, `run_control`)
- Unit tests for core statistical contracts:
  - ABF sanity
  - coloc posterior normalization
  - palindromic SNP handling
  - workflow rule contract checks

## Output Artifacts
Main reproducible outputs are generated to:
- `data_processed/tables/`
- `data_processed/final/`
- `data_processed/sensitivity/`
- `figures/` and `figures/presentation/`
- `reports/run_qc.html`
- `reports/presentation_tables.html`

## Interpretation Framework
Conclusions are reported under predefined labels:
- `Supported gut-brain genetic contribution`
- `Partial support`
- `Insufficient support`

A hard guardrail is applied: results do not constitute proof of a purely psychosomatic origin of IBS.
