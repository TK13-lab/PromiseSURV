# PROMISE-Surv reviewer code

This folder contains the concise code package for reproducing revised manuscript Figures 1–6 and auditing the main model/statistical benchmarks.

The publication R files are deliberately standalone. 

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

