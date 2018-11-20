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
    "$DATA_DIR" \
    "$OUTPUT_DIR"/mimix-no-factors

mkdir -p "$OUTPUT_DIR"/mimix
julia scripts/fit-mcmc.jl \
   --hyper nutnet-analysis/configs/hyper.yml nutnet-analysis-configs/nu-hyper-1.yml\
   --inits nutnet-analysis/configs/inits.yml \
   --monitor nutnet-analysis/configs/monitor-mimix.yml \
   --post-pred-check \
   --factors $FACTORS \
   --iters $MCMC_ITERS \
   --burnin $MCMC_BURNIN \
   --thin $MCMC_THIN \
   --chains $MCMC_CHAINS \
   "$DATA_DIR" \
   "$OUTPUT_DIR"/mimix

Rscript nutnet-analysis/summarize-nutnet-analysis.R "$OUTPUT_DIR" "$OUTPUT_DIR" "$DATA_DIR"