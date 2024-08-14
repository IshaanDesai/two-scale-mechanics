#!/bin/sh
set -e -u

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
$HOME/calculix-adapter/bin/ccx_preCICE -i one_element -precice-participant Meso-structure
