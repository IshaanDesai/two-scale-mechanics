#!/bin/bash
#
# Job name:
#SBATCH -J multi-precice
#
# Error and Output files
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#
# Working directory:
#SBATCH -D ./
#
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=ishaan.desai@ipvs.uni-stuttgart.de
#
# Wall clock limit:
#SBATCH --time=01:00:00
#
# Compute resources
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5

echo "SLURM_NNODES"=$SLURM_NNODES
echo "working directory="$SLURM_SUBMIT_DIR

# load the modules you need
#module purge
module load ipvs-epyc/gcc/10.2 ipvs-epyc/openmpi/4.0.4-gcc-10.2 ipvs-epyc/python/3.8.5 ub2004/libxml2/2.9.10 ub2004/boost/1.75.0
#module list

echo "Launching meso participant"

cd meso_dummy/

mpiexec -n 1 --bind-to core python3 dummy.py &> log_meso.log &

echo "Launching Micro Manager"

cd ../micro_nasmat_dummy/

LD_PRELOAD=/home/desaiin/petsc-3.15.5/arch-linux-c-opt/lib/libpetsc.so.3.15.5 mpiexec -n 4 --bind-to core python3 run_micro_manager.py --config micro-manager-config.json &> log_micro.log

echo "Meso simulation and Micro Manager have been launched"
