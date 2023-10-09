# Run this file to install the Micro Manager from source on Great Lakes
# It is assumed that you have cloned the Micro Manager and have it in your home path
# If the Micro Manager is not cloned, clone it by running: `git clone https://github.com/precice/micro-manager.git` at /home/$USER

module purge

pip3 uninstall -y micro-manager-precice
pip3 uninstall -y pyprecice
pip3 uninstall -y mpi4py

module load use.own gcc/my_gcc10
module load openmpi/4.1.6

# module load use.own gcc/gcc_libs
# module load intel/2022.1.2 impi/2021.5.1

module list

cd ~/micro-manager

echo "Installing the Micro Manager"

pip3 install --user .
