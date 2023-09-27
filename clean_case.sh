rm -fv slurm*
rm -rfv precice-profiling

echo "Cleaning macro participant"
cd macro_laminate
./clean.sh

cd ..

echo "Cleaning micro participant"
cd micro_ruc
./clean.sh
