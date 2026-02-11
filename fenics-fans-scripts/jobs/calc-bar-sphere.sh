#!/bin/bash
#
# Licenses
#SBATCH --licenses=scratch:0,work:0
#
# Job name:
#SBATCH -J calc-bar-sphere
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

source ../jobs/utils.sh
load_env

NUM_MM_RANKS=4
NUM_MM_WORKERS=4

# pwd = tsm/ffs/work/job-id
mkdir -p ./${SLURM_JOB_ID} && cd ./${SLURM_JOB_ID} || return
mkdir -p ../../output/${SLURM_JOB_NAME}/${SLURM_JOB_ID}
JOB_DIR=$(path::get_full .)
OUT_DIR=$(path::get_full ../../output/${SLURM_JOB_NAME}/${SLURM_JOB_ID})
TSM_PATH=$(path::get_tsm)
MESO_PATH=$(path::get_meso)
MICRO_PATH="${TSM_PATH}/micro-fans-bar-sphere"

cp $TSM_PATH/precice-config-fans-small-strain.xml ./precice-config.xml
cp $MESO_PATH/examples/coupled-bar/config-coupled-bar.json ./config-meso.json
cp $MICRO_PATH/PyFANS.so ./
cp $MICRO_PATH/input.json ./
cp $MICRO_PATH/micro-manager-pyfans-config-stateless.json ./micro-manager-config.json
cp $MICRO_PATH/sphere32.h5 ./

edit::meso_input ./config-meso.json "${OUT_DIR}/meso-geom" "${OUT_DIR}/meso-state" "${JOB_DIR}/precice-config.xml"
edit::precice_input ./precice-config.xml "${OUT_DIR}/precice.log" "${OUT_DIR}/precice-profiling" "${OUT_DIR}/precice-exports" "${JOB_DIR}"
edit::mm_input ./micro-manager-config.json "${JOB_DIR}/precice-config.xml" ${NUM_MM_RANKS} ${NUM_MM_WORKERS}
edit::fans_input ./input.json ${NUM_MM_WORKERS}

echo "Starting Meso Simulations"
touch "${OUT_DIR}/meso.log"
cd "${MESO_PATH}" || return
(I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_PROCESSOR_LIST=40,41,42,43,44,45,46,47 mpiexec -n 8 python macro.py "${JOB_DIR}/config-meso.json" &> "${OUT_DIR}/meso.log")&

echo "Starting Micro Simulations"
touch "${OUT_DIR}/micro.log"
cd "${JOB_DIR}" || return
I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_DOMAIN=${NUM_MM_WORKERS} mpiexec -n ${NUM_MM_RANKS} micro-manager-precice "${JOB_DIR}/micro-manager-config.json" &> "${OUT_DIR}/micro.log"

echo "Simulation with $SLURM_NNODES nodes launched."