#! /usr/bin/env bash

set -eo pipefail

while getopts "r:o:" option; do
    case "${option}" in
        r) REPLICATES=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        \?) echo "Invalid option: ${OPTARG}";;
        :) echo "Invalid option: ${OPTARG} requires an argument" 1>&2;;
    esac
done
readonly REPLICATES
readonly OUTPUT_DIR

# Define simulation study settings

settings_dir=simulation-study/configs
setting01=(high-dense.yml       low-error-var.yml  low-block-var.yml)
setting02=(high-dense.yml       low-error-var.yml high-block-var.yml)
setting03=(high-dense.yml      high-error-var.yml  low-block-var.yml)
setting04=(high-dense.yml      high-error-var.yml high-block-var.yml)
setting05=(high-dense.yml very-high-error-var.yml  low-block-var.yml)
setting06=(high-dense.yml very-high-error-var.yml high-block-var.yml)
setting07=( low-dense.yml       low-error-var.yml  low-block-var.yml)
setting08=( low-dense.yml       low-error-var.yml high-block-var.yml)
setting09=( low-dense.yml      high-error-var.yml  low-block-var.yml)
setting10=( low-dense.yml      high-error-var.yml high-block-var.yml)
setting11=( low-dense.yml very-high-error-var.yml  low-block-var.yml)
setting12=( low-dense.yml very-high-error-var.yml high-block-var.yml)
setting13=(  no-dense.yml       low-error-var.yml  low-block-var.yml)
setting14=(  no-dense.yml       low-error-var.yml high-block-var.yml)
setting15=(  no-dense.yml      high-error-var.yml  low-block-var.yml)
setting16=(  no-dense.yml      high-error-var.yml high-block-var.yml)
setting17=(  no-dense.yml very-high-error-var.yml  low-block-var.yml)
setting18=(  no-dense.yml very-high-error-var.yml high-block-var.yml)

# Simulate mimix models

mimix () {
    settings=()
    for setting in ${@:3}; do
        settings+=("$settings_dir/$setting")
    done
    results_dir=$OUTPUT_DIR/mimix-$1-factors/$2
    mkdir -p $results_dir
    parallel --no-notice julia scripts/simulate-model.jl \
        --inits   simulation-study/configs/inits.yml \
        --data    simulation-study/configs/data.yml "${settings[@]}" \
        --hyper   simulation-study/configs/hyper.yml \
        --monitor simulation-study/configs/monitor.yml \
        --factors $1 \
        --seed {} \
        $results_dir/rep{} ::: $(seq 1 $REPLICATES)
}

# Array of factors to simulate
factors=(0 20)

for factor in "${factors[@]}"; do
    mimix $factor setting01 "${setting01[@]}"
    mimix $factor setting02 "${setting02[@]}"
    mimix $factor setting03 "${setting03[@]}"
    mimix $factor setting04 "${setting04[@]}"
    mimix $factor setting05 "${setting05[@]}"
    mimix $factor setting06 "${setting06[@]}"
    mimix $factor setting07 "${setting07[@]}"
    mimix $factor setting08 "${setting08[@]}"
    mimix $factor setting09 "${setting09[@]}"
    mimix $factor setting10 "${setting10[@]}"
    mimix $factor setting11 "${setting11[@]}"
    mimix $factor setting12 "${setting12[@]}"
    mimix $factor setting13 "${setting13[@]}"
    mimix $factor setting14 "${setting14[@]}"
    mimix $factor setting15 "${setting15[@]}"
    mimix $factor setting16 "${setting16[@]}"
    mimix $factor setting17 "${setting17[@]}"
    mimix $factor setting18 "${setting18[@]}"
done

# Run PERMANOVA

permanova () {
    settings=()
    for setting in ${@:2}; do
        settings+=("$settings_dir/$setting")
    done
    results_dir=$OUTPUT_DIR/permanova/$1
    mkdir -p $results_dir
    parallel --no-notice $(which julia) scripts/simulate-model.jl \
        --data simulation-study/configs/data.yml "${settings[@]}" \
        --permanova \
        --seed {} \
        $results_dir/rep{} ::: $(seq 1 $REPLICATES)    
}

permanova setting01 "${setting01[@]}"
permanova setting02 "${setting02[@]}"
permanova setting03 "${setting03[@]}"
permanova setting04 "${setting04[@]}"
permanova setting05 "${setting05[@]}"
permanova setting06 "${setting06[@]}"
permanova setting07 "${setting07[@]}"
permanova setting08 "${setting08[@]}"
permanova setting09 "${setting09[@]}"
permanova setting10 "${setting10[@]}"
permanova setting11 "${setting11[@]}"
permanova setting12 "${setting12[@]}"
permanova setting13 "${setting13[@]}"
permanova setting14 "${setting14[@]}"
permanova setting15 "${setting15[@]}"
permanova setting16 "${setting16[@]}"
permanova setting17 "${setting17[@]}"
permanova setting18 "${setting18[@]}"

# Aggregate and summarize results

$(which julia) simulation/aggregate-simulation-study-results.jl $OUTPUT_DIR $OUTPUT_DIR
$(which Rscript) simulation/summarize-simulation-results.R $OUTPUT_DIR $OUTPUT_DIR