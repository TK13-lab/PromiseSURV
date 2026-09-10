# PROMISE-Surv reviewer code

This folder contains the concise code package for reproducing revised manuscript Figures 1–6 and auditing the main model/statistical benchmarks.

The publication R files are deliberately standalone. They do not call `source()`, do not import a project-wide plotting setup, and state the final text sizes, geometry sizes, legend settings, panel dimensions, layout ratios, and 400-dpi export values inside the relevant figure/panel code. The Python files likewise keep all seeds, hyperparameters, bootstrap counts, test choices, cancer order, pathways, and marker lists local to `main()`; there are no module-level configuration assignments.

## Contents

```text
R/figures/       one independent script for each final Figure 1–6
R/analysis/      standalone enrichment and reviewer statistical audit
python/figure_inputs/  Figure 1 tidy-panel-data builder
python/benchmarks/     internal, TCGA, GEO, and biology benchmarks/tests
docs/            input contract and statistical-test specification
outputs/         generated files only; ignored by Git
```

## Environment

Create the supplied Conda environment:

```bash
conda env create -f environment.yml
conda activate promise-surv-reviewer-code
```

Alternatively, install Python packages with `pip install -r requirements.txt` and install the R/Bioconductor packages listed in `R_PACKAGE_VERSIONS.txt`. Arial must be installed or registered for an exact publication-facing render.

## Reproduce Figures 1–6

Run from this repository root. Replace `/path/to/frozen_analysis` with the data tree described in `docs/INPUT_DATA.md`.

```bash
Rscript R/figures/Figure_1_cohort_construction.R --input-root=/path/to/frozen_analysis
Rscript R/figures/Figure_2_PROMISE_Surv_framework.R
Rscript R/figures/Figure_3_transfer_and_tcga_biology.R --input-root=/path/to/frozen_analysis
Rscript R/figures/Figure_4_external_survival_and_GEO_evidence.R --input-root=/path/to/frozen_analysis
Rscript R/figures/Figure_5_shared_marker_discovery.R --input-root=/path/to/frozen_analysis
Rscript R/figures/Figure_6_GEO_external_validation.R --input-root=/path/to/frozen_analysis
```

Each script writes a 400-dpi PNG plus vector PDF/SVG files to `outputs/figures/` by default. Use `--output-dir=PATH` to select another directory.

If the Figure 1 panel table must first be regenerated:

```bash
python python/figure_inputs/Figure_1_build_panel_data.py --input-root /path/to/frozen_analysis
```

## Reproduce the manuscript statistical audit

The following R script is the canonical route for the revised manuscript cohort table and Tables 2–4. It reads the locked student-strengthening metrics, runs the declared 5,000-resample internal uncertainty analysis and 2,000-resample external C-index analysis, and performs the paired two-sided Wilcoxon tests.

```bash
Rscript R/analysis/reviewer_statistical_analysis.R \
  --input-root=/path/to/frozen_analysis \
  --tcga-cdr-file=/path/to/TCGA-CDR-SupplementalTableS1.xlsx
```

## Run transparent benchmark recomputations

```bash
python python/benchmarks/01_internal_survival_benchmarks.py --input-root /path/to/frozen_analysis
python python/benchmarks/02_external_TCGA_benchmarks.py --input-root /path/to/frozen_analysis
python python/benchmarks/03_GEO_survival_and_marker_tests.py --input-root /path/to/frozen_analysis
python python/benchmarks/04_Figure_3_biology_tests.py --input-root /path/to/frozen_analysis
```

The Python scripts refit or recompute the benchmark families from their documented model inputs. They write machine-readable CSV result tables and a JSON record of their explicit settings under `outputs/benchmarks/`. See `docs/STATISTICAL_TESTS.md` for the test families and multiplicity handling.

The full internal benchmark refits 80 models and can be computationally intensive. To verify an installation without changing the locked defaults used by a normal run, add `--quick-check`; this runs one split and 20 bootstrap samples for each uncertainty calculation and records `"quick_check": true` in its settings JSON.

## Data availability boundary

This GitHub folder intentionally excludes raw controlled-access data, downloaded cohort mirrors, trained artifacts, and large generated outputs. The scripts accept `--input-root` so an authorized reviewer can point them to the documented frozen data tree without editing code or introducing machine-specific paths.
