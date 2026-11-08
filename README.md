# PIC Neonatal FT4–Cystatin C

Reproducibility code for a PIC database study of free thyroxine (FT4) and cystatin C in 760 hospitalized neonates before and after creatinine adjustment. Python performs the statistical analyses; R generates Figures 2–3 and S1–S3. Figure 1 uses the cohort-flowchart generator supplied with this repository.

The manuscript uses release [v1.0.0](https://github.com/doctorwinson/pic-neonatal-ft4-cystatin-c/releases/tag/v1.0.0).

## Run

Requirements: Python 3.13.9 and `requirements.txt`; R 4.5.2 with ggplot2, grid, jsonlite, svglite, ragg and pdftools.

Use an authorized [PIC 1.1.0](https://physionet.org/content/picdb/1.1.0/) PostgreSQL import with schema `pic_icu`. Set `PGPASSWORD` privately. Optional connection variables are `PGHOST`, `PGPORT`, `PGDATABASE` and `PGUSER`; defaults are localhost, 5432, pic and postgres.

```sh
python -m pip install -r requirements.txt
python analysis.py --extract
Rscript --vanilla figures.R
```

`analysis.py` executes `extract.sql` in read-only transactions, prepares the cohort, fits the models and generates Figure 1. Subsequent runs can use `python analysis.py` to reuse local source snapshots. Individual records and generated outputs are written under `results/` and remain local.

Primary estimates are in `results/primary_analysis/primary_results.json`; supplementary estimates are under `supplementary_analysis/evidence`, `verification/evidence` and `selection_analysis/evidence`. Figures are written to `results/figures/`, `results/supplementary_analysis/figures/` and `results/selection_analysis/figures/`.

The primary model retains the post-result freeze of 15 October 2026. Supplementary spline, optimizer and selection checks retain their exploratory status. The code preserves the original Python calculations, including HC3/t inference and the 2,000-resample bootstrap.

On Windows, R may require `LC_ALL=English_United States.utf8` when paths contain Chinese characters.

## Data and license

Patient-level data are not distributed. Access to PIC requires the provider's training, credentialing and data use agreement. The repository's `.gitignore` excludes generated outputs and permits only the listed source and documentation files.

The source code is released under the [MIT License](LICENSE). This license does not apply to PIC data or replace its data use agreement. Software citation metadata are in `CITATION.cff`.
