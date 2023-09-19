# Run this file to install the Micro Manager from source on Great Lakes
# It is assumed that you have cloned the Micro Manager and have it in your home path
# If the Micro Manager is not cloned, clone it by running: `git clone https://github.com/precice/micro-manager.git` at /home/$USER

module purge

# Modules may need to be updated depending on compatibility and availability on Great Lakes
module load gcc/10.3.0 openmpi/4.1.6

echo "Loaded modules:"
module list

cd ~/micro-manager

echo "Installing the Micro Manager in the user path used by pip3"

pip3 install --user .
