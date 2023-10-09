# Run this file to install preCICE from source on Great Lakes
# It is assumed that you have cloned preCICE and have it in your home path
# If preCICE is not cloned, clone it by running: git clone https://github.com/precice/precice.git at /home/$USER

module purge

# Use custom gcc to later use it alongside the intel module
module load use.own gcc/gcc_libs

module load intel/2022.1.2 impi/2021.5.1

module load cmake/3.26.3 eigen/3.4.0 boost/1.78.0

module list

echo "Libxml2 is installed on ARC HPC clusters, no need to load a module."

cd ~
rm -rfv precice-installation
mkdir precice-installation

cd precice

rm -rfv build
mkdir build && cd build

cmake -DCMAKE_INSTALL_PREFIX=~/precice-installation -DCMAKE_CXX_FLAGS=" -std=c++17" -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort -DPRECICE_FEATURE_MPI_COMMUNICATION=ON -DPRECICE_FEATURE_PETSC_MAPPING=OFF -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF ..

make -j 8

#ctest

make install
