# PROMISE-Surv

Analysis and figure code for PROMISE-Surv. See the [analysis overview](docs/README.md) for the scope and reproducibility requirements. Raw data, fitted models and generated outputs are excluded.

## Environment

```bash
conda env create -f environment.yml
conda activate promise-surv-reviewer-code
```

Alternatively, use `requirements.txt` and `R_PACKAGE_VERSIONS.txt`. Publication figures use Arial. Run commands from the repository root, replacing `/path/to/frozen_analysis` with your input directory.

## Figures

| Result | Script |
| --- | --- |
| Figure 1: cohort construction | `R/figures/Figure_1_cohort_construction.R` |
| Figure 2: framework | `R/figures/Figure_2_PROMISE_Surv_framework.R` |
| Figure 3: biological transfer | `R/figures/Figure_3_transfer_and_tcga_biology.R` |
| Figure 4: external survival | `R/figures/Figure_4_external_survival_and_GEO_evidence.R` |
| Figure 5: shared markers | `R/figures/Figure_5_shared_marker_discovery.R` |
| Figure 6: GEO validation | `R/figures/Figure_6_GEO_external_validation.R` |

Each R figure script accepts `--input-root=PATH` and `--output-dir=PATH`, except Figure 2, which needs no input data. Outputs are 400-dpi PNG and vector PDF/SVG files.

```bash
Rscript R/figures/Figure_1_cohort_construction.R --input-root=/path/to/frozen_analysis
Rscript R/figures/Figure_2_PROMISE_Surv_framework.R
```

`R/analysis/shared_marker_enrichment.R --input-root=PATH` regenerates the enrichment intermediates used by Figures 5–6.

## Statistical results

```bash
Rscript R/analysis/reviewer_statistical_analysis.R \
  --input-root=/path/to/frozen_analysis \
  --tcga-cdr-file=/path/to/TCGA-CDR-SupplementalTableS1.xlsx
```

This script assembles Tables 1–3 from frozen inputs, including the declared 5,000 split-level and 2,000 patient-level bootstrap calculations. Its output named `table4_external_cox_associations.csv` contains continuous TCGA Cox associations and is separate from manuscript Table 4, which combines TCGA and GEO survival contrasts.

Additional computations accept `--input-root PATH`:

| Script under `python/benchmarks/` | Purpose |
| --- | --- |
| `01_internal_survival_benchmarks.py` | Fixed-penalty model refit; distinct from the training-selected blended analysis underlying Table 2 |
| `02_external_TCGA_benchmarks.py` | TCGA concordance, Cox and log-rank tests |
| `03_GEO_survival_and_marker_tests.py` | GEO OS/DFS survival and marker contrasts |
| `04_Figure_3_biology_tests.py` | Figure 3 correlation tests |

For an installation check, the internal Python benchmark offers `--quick-check` (one split and 20 bootstrap resamples). These settings are insufficient for manuscript estimates.

See the [analysis overview](docs/README.md) and the manuscript's Supplementary Methods for methodological context and reproducibility considerations.

