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

configs_dir=nutnet-analysis/configs
nu_hypers=(nu-hyper-1.yml nu-hyper-2.yml nu-hyper-3.yml)

mimix () {
    echo -e "Sensitivity analysis: running with nu hyperparameters defined in $configs_dir/$1"
    results_dir=$OUTPUT_DIR/$(basename "$1" .yml)
    mkdir -p $results_dir
    julia scripts/fit-mcmc.jl \
        --hyper nutnet-analysis/configs/hyper.yml $configs_dir/$1 \
        --inits nutnet-analysis/configs/inits.yml \
        --monitor nutnet-analysis/configs/monitor-mimix.yml \
        --post-pred-check \
        --factors $FACTORS \
        --iters $MCMC_ITERS \
        --burnin $MCMC_BURNIN \
        --thin $MCMC_THIN \
        --chains $MCMC_CHAINS \
        --post-pred-check \
        "$DATA_DIR" \
        $results_dir
    Rscript nutnet-analysis/summarize-nutnet-analysis.R \
        $results_dir \
        $results_dir \
        "$DATA_DIR"
}

for nu_hyper in "${nu_hypers[@]}"; do
    mimix $nu_hyper
done
