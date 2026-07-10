#!/bin/bash
#
# Licenses
#SBATCH --licenses=scratch:0,work:0
#
# Job name:
#SBATCH -J calc-rivet1734-m32-loose-acc
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
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#
#Setup of execution environment
#SBATCH --get-user-env
#SBATCH --account=pn76so
#SBATCH --partition=general
#
#Number of nodes and MPI tasks per node:
#SBATCH --nodes=17
#SBATCH --ntasks-per-node=48
#
#Ensure exclusive access to compute nodes
#SBATCH --exclusive

source ../jobs/utils.sh
load_env

NUM_MESO_NODES=1
NUM_MM_NODES=16
NUM_NODES=${SLURM_NNODES}
NUM_MM_RANKS=768
NUM_MM_WORKERS=1
NUM_MM_PPN=48

# pwd = tsm/ffs/work/job-id
mkdir -p ./${SLURM_JOB_ID} && cd ./${SLURM_JOB_ID} || return
mkdir -p ../../output/${SLURM_JOB_NAME}/${SLURM_JOB_ID}
JOB_DIR=$(path::get_full .)
OUT_DIR=$(path::get_full ../../output/${SLURM_JOB_NAME}/${SLURM_JOB_ID})
TSM_PATH=$(path::get_tsm)
MESO_PATH=$(path::get_full "../../../meso-fenics-loose-acc")
MICRO_PATH="${TSM_PATH}/micro-fans-rivet"

cp $TSM_PATH/precice-config-fans-small-strain-acc.xml ./precice-config.xml
cp $MESO_PATH/examples/coupled-rivet/config-coupled-rivet.json ./config-meso.json
cp $MICRO_PATH/PyFANS.so ./
cp $MICRO_PATH/input.json ./
cp $MICRO_PATH/micro-manager-pyfans-config-stateless.json ./micro-manager-config.json
cp $MICRO_PATH/sphere3p32.h5 ./

edit::meso_input_checkpoint ./config-meso.json "${OUT_DIR}/meso-geom" "${OUT_DIR}/meso-state" "${JOB_DIR}/precice-config.xml" "${OUT_DIR}/checkpoint_rivet"
edit::precice_input ./precice-config.xml "${OUT_DIR}/precice.log" "${OUT_DIR}/precice-profiling" "${OUT_DIR}/precice-exports" "${JOB_DIR}" 8 "${NUM_MM_RANKS}"
edit::mm_input ./micro-manager-config.json "${JOB_DIR}/precice-config.xml" ${NUM_MM_RANKS} ${NUM_MM_WORKERS}
edit::fans_input ./input.json ${NUM_MM_WORKERS}
gen_host_files ${NUM_MESO_NODES} ${NUM_MM_NODES}

touch "${OUT_DIR}/meso.log"
touch "${OUT_DIR}/micro.log"
echo "Starting Simulations"
set -m
(
    (cd "${MESO_PATH}" && (I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_PROCESSOR_LIST=40,41,42,43,44,45,46,47 mpiexec -hostfile "${JOB_DIR}/hosts.meso" -n 8 -ppn 8 python macro.py "${JOB_DIR}/config-meso.json" &> "${OUT_DIR}/meso.log"))&
    (cd "${JOB_DIR}" && (I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_PIN_DOMAIN=${NUM_MM_WORKERS} I_MPI_PIN_ORDER=compact mpiexec -hostfile "${JOB_DIR}/hosts.micro" -n ${NUM_MM_RANKS} -ppn ${NUM_MM_PPN} micro-manager-precice "${JOB_DIR}/micro-manager-config.json" &> "${OUT_DIR}/micro.log"))&
    wait
)
echo "Simulations with $SLURM_NNODES nodes launched."

save_inputs "${OUT_DIR}"
cd ..
rm -rf ./${SLURM_JOB_ID}