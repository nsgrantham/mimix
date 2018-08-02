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

# simulation study settings
settings_dir=simulation/simulation-study-settings
setting1=( low-dense.yml  low-error-var.yml  low-block-var.yml)
setting2=( low-dense.yml  low-error-var.yml high-block-var.yml)
setting3=( low-dense.yml high-error-var.yml  low-block-var.yml)
setting4=( low-dense.yml high-error-var.yml high-block-var.yml)
setting5=(high-dense.yml  low-error-var.yml  low-block-var.yml)
setting6=(high-dense.yml  low-error-var.yml high-block-var.yml)
setting7=(high-dense.yml high-error-var.yml  low-block-var.yml)
setting8=(high-dense.yml high-error-var.yml high-block-var.yml)

# factors to simulate
factors=(0 20)

mimix () {
    settings=()
    for setting in ${@:3}; do
        settings+=("$settings_dir/$setting")
    done
    results_dir=$OUTPUT_DIR/$2/mimix-$1-factors
    mkdir -p $results_dir
    parallel --no-notice $(which julia) simulation/simulate-model.jl \
        --inits   simulation/inits.yml \
        --data    simulation/data.yml "${settings[@]}" \
        --hyper   simulation/hyper.yml \
        --monitor simulation/monitor.yml \
        --factors $1 \
        --seed {} \
        $results_dir/rep{}.tsv ::: $(seq 1 $REPLICATES)
}

for factor in "${factors[@]}"; do
    mimix $factor setting1 "${setting1[@]}"
    mimix $factor setting2 "${setting2[@]}"
    mimix $factor setting3 "${setting3[@]}"
    mimix $factor setting4 "${setting4[@]}"
    mimix $factor setting5 "${setting5[@]}"
    mimix $factor setting6 "${setting6[@]}"
    mimix $factor setting7 "${setting7[@]}"
    mimix $factor setting8 "${setting8[@]}"
done