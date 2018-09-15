# MIMIX

MIMIX (MIcrobiome MIXed model) is a hierarchical Bayesian model for the analysis of next-generation sequencing microbiome abundance data from designed experiments. It achieves four scientific objectives:

1. Global tests of whether experimental treatments affect microbiome composition, 
2. Local tests for treatment effects on individual taxa and estimation of theses effects if present,
3. Quantification of how different sources of variability contribute to microbiome heterogeneity, and
4. Characterization of latent structure in the microbiome, which may suggest ecological subcommunities.

For more information, read the paper: [MIMIX: a Bayesian Mixed-Effects Model for Microbiome Data from Designed Experiments](https://arxiv.org/abs/1703.07747).

## Installation

MIMIX requires Julia 0.6, which can be downloaded from [https://julialang.org/downloads/oldreleases.html](https://julialang.org/downloads/oldreleases.html). Platform specific installation instuctions are available at [https://julialang.org/downloads/platform.html](https://julialang.org/downloads/platform.html).

After Julia is successfully installed, clone this repository and symlink it to Julia's package directory.

```bash
ln -s /path/to/mimix ~/.julia/v0.6/MicrobiomeMixedModels.jl
```

The following Julia packages are required:
- StatsBase
- DataFrames
- Distributions
- Mamba
- YAML
- RCall
- CSV
- ArgParse

To install each Julia package, run `Pkg.add("PackageName")` within an active Julia session.

The analyses in the paper also require R (3.4 or greater) and several R packages:
- tidyverse
- vegan
- argparse
- xtable
- reshape2

To install each R package, run `install.packages("PackageName")` within an active R session.

Finally, the simulation study requires the GNU tool `parallel` (read more at [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)) which can be installed with `brew install parallel` on macOS and `sudo apt-get install parallel` on Ubuntu.

## Simulation study

Run `./simulation-study/run-simulation-study.sh -r 50 -o simulation-study/results` to reproduce the simulation study.

This depends primarily on `scripts/sim-mcmc.jl` which allows the user to simulate a model with their own configurations. Example configurations can be found in `simulation-study/configs`. A call to this script may look like the following:

```
julia scripts/sim-mcmc.jl \
    --data path/to/data.yml \
    --hyper path/to/hyper.yml \
    --monitor path/to/monitor.yml \
    --inits path/to/inits.yml \
    --seed 123 \
    --factors 20 \
    path/to/output_dir
```

## Data analysis

Run `./nutnet-analysis/run-nutnet-analysis.sh -d nutnet-analysis/full-data -o nutnet-analysis -f 500 -i 20000 -b 10000 -t 1 -c 1` to reproduce the full NutNet data analysis. To demo the analysis on a dataset of reduced dimensionality, replace `nutnet-analysis/full-data` with `nutnet-analysis/reduced-data` and reduce the number of factors to `-f 100` or fewer.

This depends primarily on `scripts/fit-mcmc.jl` which allows the user to fit a model with their own configurations. Example configurations can be found in `nutnet-analysis/configs`. A call to this script may look like the following:

```
julia scripts/fit-mcmc.jl \
    --hyper path/to/hyper.yml \
    --monitor path/to/monitor.yml \
    --inits path/to/inits.yml \
    --factors 20 \
    path/to/data_dir \
    path/to/output_dir
```

The data directory must contain three files:
- `X.csv`: treatment covariates in samples (rows) by covariates (columns)
- `Y.csv`: microbiome abundance data in samples (rows) by taxa (columns)
- `Z.csv`: block identifiers in samples (rows) by number of blocking factors (columns)
