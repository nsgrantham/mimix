# MIMIX

MIMIX (MIcrobiome MIXed model) is a hierarchical Bayesian model for the analysis of next-generation sequencing microbiome abundance data from designed experiments. It achieves four scientific objectives:

1. Global tests of whether experimental treatments affect microbiome composition, 
2. Local tests for treatment effects on individual taxa and estimation of theses effects if present,
3. Quantification of how different sources of variability contribute to microbiome heterogeneity, and
4. Characterization of latent structure in the microbiome, which may suggest ecological subcommunities.

For more information, read the paper: [MIMIX: a Bayesian Mixed-Effects Model for Microbiome Data from Designed Experiments](https://arxiv.org/abs/1703.07747).

## Installation

### Install Julia

MIMIX requires Julia 1.0, which can be downloaded from [https://julialang.org/downloads/](https://julialang.org/downloads/). Platform specific installation instuctions are available at [https://julialang.org/downloads/platform.html](https://julialang.org/downloads/platform.html).

Add the `julia` command to your path. One way to achieve this on macOS is to symlink the command to `/usr/local/bin/`:

```
ln -s /Applications/Julia-1.0.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

### Install MicrobiomeMixedModels.jl

```
git clone https://github.com/nsgrantham/mimix
cd mimix
```

Open Julia with `julia`, type `]` to enter Pkg mode, and then run the following commands:

```
develop MicrobiomeMixedModels.jl
add StatsBase DataFrames Distributions Mamba YAML RCall CSV ArgParse
```


### Install R

MIMIX requires R 3.4 or greater to reproduce the tables and graphics from the paper and to compare the performance of MIMIX and its variants with PERMANOVA, implemented by `vegan::adonis()`.

Download and install R from [https://www.r-project.org/](https://www.r-project.org/).

Once installed, open R from the terminal with `R` and run the following command:

```
install.packages(c("tidyverse", "vegan", "argparse", "xtable", "reshape2"))
```

### Install parallel

To reproduce the simulation study, you will need to install the GNU tool `parallel` (read more at [https://www.gnu.org/software/parallel/](https://www.gnu.org/software/parallel/)).

To do so, run `brew install parallel` on macOS or `sudo apt-get install parallel` on Ubuntu.


## Simulation study

Run `scripts/sim-mcmc.jl` to simulate a model on artificially-generated data. A call to this script may look like the following:

```
mkdir simulation-results
julia scripts/sim-mcmc.jl \
    --data simulation-study/configs/data.yml \
    --hyper simulation-study/configs/hyper.yml \
    --monitor simulation-study/configs/monitor.yml \
    --inits simulation-study/configs/inits.yml \
    --seed 123 \
    --factors 20 \
    simulation-results
```

Run `./simulation-study/run-simulation-study.sh -r 50 -o simulation-study/results` to reproduce the simulation study with 50 replicates. A single simulation replicate, which fits all three models (MIMIX, MIMIX w/o Factors, and PERMANOVA), takes around 20 minutes to run on a single core of a 2.3 GHz Intel Core i5 processor. Therefore, a full run of the simulation study with 50 replicates of each of 30 settings parallelized across 20 cores requires about 30 hours.


## Data analysis

Run `scripts/fit-mcmc.jl` to fit a model to microbiome data from a designed experiment. A call to this script may look like the following:

```
mkdir nutnet-analysis-results
julia scripts/fit-mcmc.jl \
    --hyper nutnet-analysis/configs/hyper.yml \
    --monitor nutnet-analysis/configs/monitor-mimix.yml \
    --inits nutnet-analysis/configs/inits.yml \
    --factors 20 \
    nutnet-analysis/reduced-data \
    nutnet-analysis-results
```

The data directory must contain three files:
- `X.csv`: treatment covariates in `n` samples (rows) by `p` covariates (columns)
- `Y.csv`: microbiome abundance data in `n` samples (rows) by `K` taxa (columns)
- `Z.csv`: block identifiers in `n` samples (rows) by `q` blocking factors (columns)

Run `./nutnet-analysis/run-nutnet-analysis.sh -d nutnet-analysis/full-data -o nutnet-analysis -f 166 -i 20000 -b 10000 -t 20 -c 1` to reproduce the full NutNet data analysis, which will take aabout 26 hours on a 2.3 GHz Intel Core i5 processor.

To demo the analysis on a small dataset, replace `nutnet-analysis/full-data` with `nutnet-analysis/reduced-data` and select the number of factors (`-f`) to be 100 or fewer.


