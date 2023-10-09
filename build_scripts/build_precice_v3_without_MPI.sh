# Run this file to install preCICE from source on Great Lakes
# It is assumed that you have cloned preCICE and have it in your home path
# If preCICE is not cloned, clone it by running: git clone https://github.com/precice/precice.git at /home/$USER

module purge

# Use custom gcc to later use it alongside the intel module
module load use.own gcc/my_gcc10

# Modules may need to be updated depending on compatibility and availability on Great Lakes
module load cmake/3.26.3 eigen/3.4.0 boost/1.78.0

echo "Loaded modules:"
module list

echo "Libxml2 is installed on ARC HPC clusters, no need to load a module."

cd ~
rm -rfv precice-installation
mkdir precice-installation

cd precice

rm -rfv build
mkdir build && cd build

cmake -DCMAKE_INSTALL_PREFIX=~/precice-installation -DPRECICE_FEATURE_MPI_COMMUNICATION=OFF -DPRECICE_FEATURE_PETSC_MAPPING=OFF -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF ..

make -j 8

ctest

make install
