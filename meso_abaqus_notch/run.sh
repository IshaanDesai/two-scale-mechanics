#!/bin/sh
set -e -u

#num_elem=$1

#echo "Cleaning old CalculiX input files (if there are any)"
#rm -rfv *_elements.inp

#echo "Creating CalculiX input file with ${num_elem}"
#python3 generate_calculix_input.py $num_elem

echo "Using the CalculiX adapter configuration file for the coupling with the micro-mechanics solver FANS"
cp config_for_fans.yml config.yml

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1

#$HOME/calculix-adapter/bin/ccx_preCICE -i ${num_elem}_elements -precice-participant Meso-structure

#$HOME/calculix-adapter/bin/ccx_preCICE -i AbqTest1 -precice-participant Meso-structure

$HOME/calculix-adapter/bin/ccx_preCICE -i one_element -precice-participant Meso-structure

#gdb -ex run --args $HOME/calculix-adapter/bin/ccx_preCICE -i one_element -precice-participant Meso-structure


rm -rfv config.yml

#echo "Cleaning CalculiX input file"
#rm -rfv *_elements.inp
