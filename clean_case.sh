rm -fv slurm*
rm -rfv ./precice-profiling/
rm -rfv ./precice-run/
rm -rfv *.err
rm -rfv *.out

echo "Cleaning meso Abaqus participant"
cd meso_abaqus
./clean.sh

echo "Cleaning meso Abaqus participant"
cd ../meso_dummy
./clean.sh

echo "Cleaning micro Abaqus participant"
cd ../micro_abaqus
./clean.sh

echo "Cleaning micro NASMAT participant"
cd ../micro_nasmat
./clean.sh

echo "Cleaning NASMAT micro participant"
cd ../micro_nasmat_dummy
./clean.sh
