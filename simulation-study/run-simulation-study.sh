#! /usr/bin/env bash

set -euo pipefail

while getopts "r:o:" option; do
    case "${option}" in
        r) REPLICATES=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        \?) echo "Invalid option: ${OPTARG}";;
        :) echo "Invalid option: ${OPTARG} requires an argument" 1>&2;;
    esac
done

# Define simulation study settings

settings_dir=simulation-study/configs
#setting01=(high-dense.yml       low-error-var.yml  low-block-var.yml grouped-form.yml)
#setting02=(high-dense.yml       low-error-var.yml high-block-var.yml grouped-form.yml)
#setting03=(high-dense.yml      high-error-var.yml  low-block-var.yml grouped-form.yml)
#setting04=(high-dense.yml      high-error-var.yml high-block-var.yml grouped-form.yml)
#setting05=(high-dense.yml very-high-error-var.yml  low-block-var.yml grouped-form.yml)
#setting06=(high-dense.yml very-high-error-var.yml high-block-var.yml grouped-form.yml)
#setting07=( low-dense.yml       low-error-var.yml  low-block-var.yml grouped-form.yml)
#setting08=( low-dense.yml       low-error-var.yml high-block-var.yml grouped-form.yml)
#setting09=( low-dense.yml      high-error-var.yml  low-block-var.yml grouped-form.yml)
#setting10=( low-dense.yml      high-error-var.yml high-block-var.yml grouped-form.yml)
#setting11=( low-dense.yml very-high-error-var.yml  low-block-var.yml grouped-form.yml)
#setting12=( low-dense.yml very-high-error-var.yml high-block-var.yml grouped-form.yml)
#setting13=(high-dense.yml       low-error-var.yml  low-block-var.yml  random-form.yml)
#setting14=(high-dense.yml       low-error-var.yml high-block-var.yml  random-form.yml)
#setting15=(high-dense.yml      high-error-var.yml  low-block-var.yml  random-form.yml)
#setting16=(high-dense.yml      high-error-var.yml high-block-var.yml  random-form.yml)
#setting17=(high-dense.yml very-high-error-var.yml  low-block-var.yml  random-form.yml)
#setting18=(high-dense.yml very-high-error-var.yml high-block-var.yml  random-form.yml)
#setting19=( low-dense.yml       low-error-var.yml  low-block-var.yml  random-form.yml)
#setting20=( low-dense.yml       low-error-var.yml high-block-var.yml  random-form.yml)
#setting21=( low-dense.yml      high-error-var.yml  low-block-var.yml  random-form.yml)
#setting22=( low-dense.yml      high-error-var.yml high-block-var.yml  random-form.yml)
#setting23=( low-dense.yml very-high-error-var.yml  low-block-var.yml  random-form.yml)
#setting24=( low-dense.yml very-high-error-var.yml high-block-var.yml  random-form.yml)
setting25=(  no-dense.yml       low-error-var.yml  low-block-var.yml grouped-form.yml)
setting26=(  no-dense.yml       low-error-var.yml high-block-var.yml grouped-form.yml)
setting27=(  no-dense.yml      high-error-var.yml  low-block-var.yml grouped-form.yml)
setting28=(  no-dense.yml      high-error-var.yml high-block-var.yml grouped-form.yml)
setting29=(  no-dense.yml very-high-error-var.yml  low-block-var.yml grouped-form.yml)
setting30=(  no-dense.yml very-high-error-var.yml high-block-var.yml grouped-form.yml)

# Run PERMANOVA

permanova () {
    settings=()
    for setting in ${@:2}; do
        settings+=("$settings_dir/$setting")
    done
    results_dir=$OUTPUT_DIR/permanova/$1
    mkdir -p $results_dir
    parallel --no-notice /Applications/Julia-1.0.app/Contents/Resources/julia/bin/julia scripts/sim-mcmc.jl \
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
permanova setting19 "${setting19[@]}"
permanova setting20 "${setting20[@]}"
permanova setting21 "${setting21[@]}"
permanova setting22 "${setting22[@]}"
permanova setting23 "${setting23[@]}"
permanova setting24 "${setting24[@]}"
permanova setting25 "${setting25[@]}"
permanova setting26 "${setting26[@]}"
permanova setting27 "${setting27[@]}"
permanova setting28 "${setting28[@]}"
permanova setting29 "${setting29[@]}"
permanova setting30 "${setting30[@]}"

# Run MIMIX models

mimix () {
    settings=()
    for setting in ${@:4}; do
        settings+=("$settings_dir/$setting")
    done
    results_dir=$OUTPUT_DIR/mimix-$1-factors-$2-priors/$3
    mkdir -p $results_dir
    parallel --no-notice /Applications/Julia-1.0.app/Contents/Resources/julia/bin/julia scripts/sim-mcmc.jl \
        --inits   simulation-study/configs/inits.yml \
        --data    simulation-study/configs/data.yml "${settings[@]}" \
        --hyper   simulation-study/configs/hyper.yml \
        --monitor simulation-study/configs/monitor.yml \
        --factors $1 \
        --loadings $2 \
        --seed {} \
        $results_dir/rep{} ::: $(seq 1 $REPLICATES)
}

# Array of factors to simulate
factors=(0 20)

for factor in "${factors[@]}"; do
    mimix $factor DL setting01 "${setting01[@]}"
    mimix $factor DL setting02 "${setting02[@]}"
    mimix $factor DL setting03 "${setting03[@]}"
    mimix $factor DL setting04 "${setting04[@]}"
    mimix $factor DL setting05 "${setting05[@]}"
    mimix $factor DL setting06 "${setting06[@]}"
    mimix $factor DL setting07 "${setting07[@]}"
    mimix $factor DL setting08 "${setting08[@]}"
    mimix $factor DL setting09 "${setting09[@]}"
    mimix $factor DL setting10 "${setting10[@]}"
    mimix $factor DL setting11 "${setting11[@]}"
    mimix $factor DL setting12 "${setting12[@]}"
    mimix $factor DL setting13 "${setting13[@]}"
    mimix $factor DL setting14 "${setting14[@]}"
    mimix $factor DL setting15 "${setting15[@]}"
    mimix $factor DL setting16 "${setting16[@]}"
    mimix $factor DL setting17 "${setting17[@]}"
    mimix $factor DL setting18 "${setting18[@]}"
    mimix $factor DL setting19 "${setting19[@]}"
    mimix $factor DL setting20 "${setting20[@]}"
    mimix $factor DL setting21 "${setting21[@]}"
    mimix $factor DL setting22 "${setting22[@]}"
    mimix $factor DL setting23 "${setting23[@]}"
    mimix $factor DL setting24 "${setting24[@]}"
    mimix $factor DL setting25 "${setting25[@]}"
    mimix $factor DL setting26 "${setting26[@]}"
    mimix $factor DL setting27 "${setting27[@]}"
    mimix $factor DL setting28 "${setting28[@]}"
    mimix $factor DL setting29 "${setting29[@]}"
    mimix $factor DL setting30 "${setting30[@]}"
done

# Now with Gaussian priors
mimix 20 G setting01 "${setting01[@]}"
mimix 20 G setting02 "${setting02[@]}"
mimix 20 G setting03 "${setting03[@]}"
mimix 20 G setting04 "${setting04[@]}"
mimix 20 G setting05 "${setting05[@]}"
mimix 20 G setting06 "${setting06[@]}"
mimix 20 G setting07 "${setting07[@]}"
mimix 20 G setting08 "${setting08[@]}"
mimix 20 G setting09 "${setting09[@]}"
mimix 20 G setting10 "${setting10[@]}"
mimix 20 G setting11 "${setting11[@]}"
mimix 20 G setting12 "${setting12[@]}"
mimix 20 G setting13 "${setting13[@]}"
mimix 20 G setting14 "${setting14[@]}"
mimix 20 G setting15 "${setting15[@]}"
mimix 20 G setting16 "${setting16[@]}"
mimix 20 G setting17 "${setting17[@]}"
mimix 20 G setting18 "${setting18[@]}"
mimix 20 G setting19 "${setting19[@]}"
mimix 20 G setting20 "${setting20[@]}"
mimix 20 G setting21 "${setting21[@]}"
mimix 20 G setting22 "${setting22[@]}"
mimix 20 G setting23 "${setting23[@]}"
mimix 20 G setting24 "${setting24[@]}"
mimix 20 G setting25 "${setting25[@]}"
mimix 20 G setting26 "${setting26[@]}"
mimix 20 G setting27 "${setting27[@]}"
mimix 20 G setting28 "${setting28[@]}"
mimix 20 G setting29 "${setting29[@]}"
mimix 20 G setting30 "${setting30[@]}"

# Aggregate and summarize results

julia simulation-study/aggregate-results.jl $OUTPUT_DIR $OUTPUT_DIR
Rscript simulation-study/summarize-results.R $OUTPUT_DIR $OUTPUT_DIR