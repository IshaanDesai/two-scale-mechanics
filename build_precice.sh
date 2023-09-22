# Run this file to install preCICE from source on Great Lakes
# It is assumed that you have cloned preCICE and have it in your home path
# If preCICE is not cloned, clone it by running: git clone https://github.com/precice/precice.git at /home/$USER

module purge

# Modules may need to be updated depending on compatibility and availability on Great Lakes
module load gcc/10.3.0 cmake/3.26.3 eigen/3.4.0 boost/1.78.0 openmpi/4.1.6
#module load intel/18.0.5

echo "Loaded modules:"
module list

echo "Libxml2 is installed on ARC HPC clusters, no need to load a module."

cd ~/precice

rm -rfv build
mkdir build
cd build

echo "Building and installing preCICE in $PWD"

cmake -DCMAKE_BUILD_TYPE=Debug -DPRECICE_FEATURE_PETSC_MAPPING=OFF -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF ..
make -j 8
ctest
#make install
