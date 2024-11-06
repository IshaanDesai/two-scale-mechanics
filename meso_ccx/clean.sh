#!/bin/sh
set -e -u

# ccx files
rm -rf *.12d *.cvg *.dat *.frd *.sta *.out

# preCICE files
rm -rf precice-profiling* ../precice-run*

# CalculiX config.yml which is left over
rm -rf config.yml
