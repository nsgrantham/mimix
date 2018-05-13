# MIMIX

**Requires Julia 0.4**

MIMIX (MIcrobiome MIXed model) is hierarchical Bayesian model for the analysis of next-generation sequencing OTU abundance data from designed experiments. 
It achieves four scientific objectives:

1. Global tests of whether experimental treatments affect microbiome composition, 
2. Local tests for treatment effects on individual taxa and estimation of theses effects if present,
3. Quantification of how different sources of variability contribute to microbiome heterogeneity, and
4. Characterization of latent structure in the microbiome, which may suggest ecological subcommunities.

For more information, please see our open access e-print: [MIMIX: a Bayesian Mixed-Effects Model for Microbiome Data from Designed Experiments](https://arxiv.org/abs/1703.07747).

## Description of files

`analyze.jl` runs the full-scale NutNet data analysis with MIMIX, MIMIX w/o Factors, and PERMANOVA w/ Bray-Curtis.

`data` contains the NutNet data stored as `Y.csv` (sequence counts), `X.csv` (experimental treatments), and `Z.csv` (blocking factors).

`demo.ipynb` demonstrates fitting MIMIX to a subset of the real data.

`models.jl` defines the model hierarchy and posterior sampling schemes for MIMIX and MIMIX w/o Factors.

`results` contains three subdirectories, `analyze`, `simulate`, and `validate`, which store the output from their respective scripts.

`simulate.jl` runs the simulation study comparing MIMIX, MIMIX w/o Factors, and PERMANOVA w/ Bray-Curtis.

`utils.jl` defines helper functions for the simulation study and cross-validation files.

`validate.jl` runs a five-fold cross-validation to compare the fits of MIMIX and MIMIX w/o Factors to the real data.
