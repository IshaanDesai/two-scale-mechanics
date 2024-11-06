#!/bin/sh
set -e -u

echo "Configuring to couple with the micro-mechanics solver FANS"
cp config_for_fans.yml config.yml

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
$HOME/calculix-adapter/bin/ccx_preCICE -i one_element -precice-participant Meso-structure

rm config.yml
