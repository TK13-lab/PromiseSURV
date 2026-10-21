# PROMISE-Surv

PROMISE-Surv studies the transfer of protein-informed survival risk to models using RNA and clinical measurements. This repository contains statistical analysis code and scripts for selected publication figures.

## Software

Analyses use R and Python. Software requirements and recorded package versions are included in the repository. The supplied Conda environment can be created with:

```bash
conda env create -f environment.yml
conda activate promise-surv-reviewer-code
```

## Figures

Figures 1 and 2 were prepared in Microsoft PowerPoint and have no corresponding R scripts in this repository.

R scripts are provided for Figures 3–6, covering biological transfer, external survival evaluation, shared-marker discovery and GEO validation. These scripts use Arial and export raster and vector figures.

## Statistical analyses

The code covers survival evaluation, bootstrap uncertainty, model comparisons, biological associations and pathway enrichment. Python benchmarks provide additional model refitting and statistical checks. The fixed-penalty refit is distinct from the original training-selected blended analysis.

## Reproducibility

The analyses require separately prepared data and consistent software environments. Raw data, fitted models and generated results are not distributed here. Dataset details and methodological specifications are described in the manuscript and Supplementary Methods.

See the [analysis overview](docs/README.md) for further context and reproducibility considerations.
