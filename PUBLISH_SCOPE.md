FINAL PIPELINE PUBLISH BUNDLE

Scope:
- Includes only scripts/configs used by the final reproducible route (Snakemake target: all), including sensitivity and presentation layers.
- Excludes exploratory/legacy scripts and any data/artifact folders.

Excluded pipeline scripts from src/pipeline/:
- build_article_docx.py
- make_run_summary_and_sections.py
- run_safe_staged.py
- run_with_progress.py

Run entrypoints:
- snakemake --profile profiles/windows-default all
- python -m src.pipeline.run_managed --profile profiles/windows-default --configfile config/config.yaml --target all
