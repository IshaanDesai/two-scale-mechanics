rm -fv slurm*
rm -rfv precice-profiling
rm -rfv precice-run

echo "Cleaning meso participant"
cd meso_laminate
./clean.sh

cd ..

echo "Cleaning micro participant"
cd micro_ruc
./clean.sh

cd ..

echo "Cleaning NASMAT micro participant"
cd micro_ruc_nasmat
./clean.sh
