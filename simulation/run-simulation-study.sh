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

# MIMIX w/o Factors
parallel --no-notice $(which julia) simulation/simulate-model.jl \
    --inits   simulation/inits.yml \
    --data    simulation/data.yml \
    --hyper   simulation/hyper.yml \
    --monitor simulation/monitor.yml \
    --no-factors \
    --seed {} \
    $OUTPUT_DIR/mimix-no-factors-rep-{}.tsv ::: $(seq 1 $REPLICATES)

# MIMIX
parallel --no-notice $(which julia) simulation/simulate-model.jl \
    --inits   simulation/inits.yml \
    --data    simulation/data.yml \
    --hyper   simulation/hyper.yml \
    --monitor simulation/monitor.yml \
    --factors 10 \
    --seed {} \
    $OUTPUT_DIR/mimix-10-factors-rep-{}.tsv ::: $(seq 1 $REPLICATES)