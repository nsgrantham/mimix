#! /usr/bin/env bash

set -euo pipefail

while getopts "d:o:f:i:b:t:c:" option; do
    case "${option}" in
        d) DATA_DIR=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        f) FACTORS=${OPTARG};;
        i) MCMC_ITERS=${OPTARG};;
        b) MCMC_BURNIN=${OPTARG};;
        t) MCMC_THIN=${OPTARG};;
        c) MCMC_CHAINS=${OPTARG};;
        \?) echo "Invalid option: ${OPTARG}";;
        :) echo "Invalid option: ${OPTARG} requires an argument" 1>&2;;
    esac
done

mkdir -p "$OUTPUT_DIR"/permanova
julia scripts/fit-mcmc.jl \
   --permanova \
   "$DATA_DIR" \
   "$OUTPUT_DIR"/permanova

mkdir -p $OUTPUT_DIR/mimix-no-factors
julia scripts/fit-mcmc.jl \
   --hyper nutnet-analysis/configs/hyper.yml \
   --inits nutnet-analysis/configs/inits.yml \
   --monitor nutnet-analysis/configs/monitor-mimix-no-factors.yml \
   --factors 0 \
   --iters $MCMC_ITERS \
   --burnin $MCMC_BURNIN \
   --thin $MCMC_THIN \
   --chains $MCMC_CHAINS \
   --post-pred-check \
   "$DATA_DIR" \
   "$OUTPUT_DIR"/mimix-no-factors

Rscript nutnet-analysis/summarize-nutnet-analysis.R "$OUTPUT_DIR"/mimix-no-factors "$OUTPUT_DIR"/mimix-no-factors "$DATA_DIR"

mkdir -p "$OUTPUT_DIR"/mimix-gaussian
julia scripts/fit-mcmc.jl \
   --hyper nutnet-analysis/configs/hyper.yml \
   --inits nutnet-analysis/configs/inits.yml \
   --monitor nutnet-analysis/configs/monitor-mimix.yml \
   --post-pred-check \
   --factors $FACTORS \
   --iters $MCMC_ITERS \
   --burnin $MCMC_BURNIN \
   --thin $MCMC_THIN \
   --chains $MCMC_CHAINS \
   --loadings 'G' \
   --post-pred-check \
   "$DATA_DIR" \
   "$OUTPUT_DIR"/mimix-gaussian

Rscript nutnet-analysis/summarize-nutnet-analysis.R "$OUTPUT_DIR"/mimix-gaussian "$OUTPUT_DIR"/mimix-gaussian "$DATA_DIR"

mkdir -p "$OUTPUT_DIR"/mimix
julia scripts/fit-mcmc.jl \
   --hyper nutnet-analysis/configs/hyper.yml \
   --inits nutnet-analysis/configs/inits.yml \
   --monitor nutnet-analysis/configs/monitor-mimix.yml \
   --post-pred-check \
   --factors $FACTORS \
   --iters $MCMC_ITERS \
   --burnin $MCMC_BURNIN \
   --thin $MCMC_THIN \
   --chains $MCMC_CHAINS \
   --loadings 'DL' \
   "$DATA_DIR" \
   "$OUTPUT_DIR"/mimix

Rscript nutnet-analysis/summarize-nutnet-analysis.R "$OUTPUT_DIR"/mimix "$OUTPUT_DIR"/mimix "$DATA_DIR"