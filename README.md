# Gut-Brain IBS Genetics Pipeline

Windows-first Snakemake workflow for gut-brain IBS analysis.

## Quick start

1. Create environment:
   - `mamba env create -f envs/environment.yml`
2. Activate:
   - `conda activate gut-brain-ibs`
3. Run main workflow:
   - `snakemake --profile profiles/windows-default all`
4. Run sensitivity workflow:
   - `snakemake --profile profiles/windows-sensitivity sensitivity_all`

## Managed execution (recommended)

Start with persistent status tracking:

- `python -m src.pipeline.run_managed --profile profiles/windows-default --configfile config/config.yaml --target all --heartbeat-sec 15`

Read current status:

- `python -m src.pipeline.run_control status --status-file reports/run_status.json`

Watch live status with subprogress:

- `python -m src.pipeline.run_control watch --status-file reports/run_status.json --interval-sec 10`

Stop running pipeline and unlock workdir:

- `python -m src.pipeline.run_control stop --pid-file reports/run_status.pid --workspace . --unlock`

Kill stale orphan processes from previous interrupted runs:

- `python -m src.pipeline.run_control cleanup --workspace . --unlock`

## Required configs

- `config/config.yaml`
- `config/gwas_manifest.yaml`
- `config/traits.yaml`

## Mandatory external GWAS input

Add GWAS summary stats files matching `config/gwas_manifest.yaml` columns.

## Main outputs

- `data_processed/tables/`
- `data_processed/tables/presentation/`
- `figures/`
- `figures/presentation/`
- `reports/run_qc.html`
- `reports/presentation_tables.html`
- `data_processed/final/hypothesis_scorecard.json`
- `data_processed/final/evidence_summary.tsv`
- `data_processed/sensitivity/`
