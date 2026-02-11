#!/bin/bash
#
# Licenses
#SBATCH --licenses=scratch:0,work:0
#
# Job name:
#SBATCH -J test-coupled-notch
#
# Error and Output files
#SBATCH -o ../output/%x/%j.out
#SBATCH -e ../output/%x/%j.err
#
# Working directory:
#SBATCH -D ../work
#
#Notification and type
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alex.hocks@tum.de
#
# Wall clock limit:
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#
#Setup of execution environment
#SBATCH --get-user-env
#SBATCH --account=pn76so
#SBATCH --partition=test
#
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#
#Ensure exclusive access to compute nodes
#SBATCH --exclusive

# Important
echo "Loading the necessary modules"
module load slurm_setup
module load stack/24.5.0
module load intel-toolkit
module switch intel/2025.1.1 gcc/14.3.0
module load boost/1.84.0-gcc14-impi
module load eigen/3.4.0-gcc14
module load fftw/3.3.10-gcc14-impi-openmp
module load hdf5/1.14.5-gcc14-impi
module load hypre/2.32.0-gcc14-impi
#module load mumps/5.7.3-gcc14-impi-openmp-shared
#module load petsc/3.22.1-gcc14-impi-real
module load scotch/7.0.4-gcc14-impi-i64
module load parmetis/4.0.3-gcc14-impi-i64-r64
module load autoconf/
module load glib
module load pkg-config/
module load zlib
module load python/3.10.12-extended
export SLURM_EAR_LOAD_MPI_VERSION="intel"

source ~/venv/py310/bin/activate

echo "Starting meso sim"
touch ./output/${SLURM_JOB_NAME}/meso_${SLURM_JOB_ID}.log
(I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_PROCESSOR_LIST=40,41,42,43,44,45,46,47 mpiexec -n 8 python macro.py examples/coupled-notch/config-coupled-notch.json &> ./output/${SLURM_JOB_NAME}/meso_${SLURM_JOB_ID}.log)&

echo "Starting micro sim"
touch ./output/${SLURM_JOB_NAME}/micro_${SLURM_JOB_ID}.log
(cd ../micro-fans-notch-sphere && (I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_PROCESSOR_LIST=0,8,16,24 mpiexec -n 4 micro-manager-precice micro-manager-pyfans-config-stateless.json &> ../meso-fenics/output/${SLURM_JOB_NAME}/micro_${SLURM_JOB_ID}.log))

echo "Simulation with $SLURM_NNODES nodes launched."