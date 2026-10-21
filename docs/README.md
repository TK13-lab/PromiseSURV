# Analysis overview

PROMISE-Surv investigates the transfer of protein-informed survival risk to models that use RNA and clinical measurements. The analyses cover internal model comparisons, external survival evaluation, biological associations, shared-marker discovery and publication figures.

## Data and software

The study uses CPTAC/PDC discovery data and TCGA and GEO validation resources. Dataset identifiers, retrieval dates, access conditions and methodological details are provided in the manuscript and Supplementary Methods. Data must be obtained and prepared separately under the source repositories' conditions.

The repository provides R and Python analysis code with software environment specifications. Statistical analyses include survival modeling, bootstrap uncertainty estimation, model comparisons, correlation tests and pathway enrichment.

## Reproducibility

The scripts use prepared cohort data and precomputed analysis results. Reproduction requires consistent cohort definitions, preprocessing, variable encoding, data ordering and software versions. Random seeds control the specified stochastic steps, but numerical results may vary across software libraries, hardware and execution settings.

The supplied fixed-penalty model refit is distinct from the original training-selected blended analysis. Use the manuscript and Supplementary Methods to distinguish these analyses and interpret their results. The repository supports statistical evaluation and figure generation; it does not contain the complete data-acquisition and original model-training workflow.
