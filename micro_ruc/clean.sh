rm -rfv output/
rm -rfv ruc_*/
rm -fv *.log
rm -rfv __pycache__
rm -rfv precice-profiling/
rm -fv .nfs*

echo "Cleaning restart_method/"
cd restart_method/
rm -rfv ruc_*/
rm -rfv __pycache__/
rm -fv *.log
