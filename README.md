# PROMISE-Surv code

This folder contains the concise code package for reproducing revised manuscript Figures 1–6 and auditing the main model/statistical benchmarks.

The publication R files are deliberately standalone. 


## Environment

Create the supplied Conda environment:

```bash
conda env create -f environment.yml
conda activate promise-surv-reviewer-code
```

Alternatively, install Python packages with `pip install -r requirements.txt` and install the R/Bioconductor packages listed in `R_PACKAGE_VERSIONS.txt`. Arial must be installed or registered for an exact publication-facing render.

