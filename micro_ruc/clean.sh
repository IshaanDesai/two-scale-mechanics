rm -rfv output/
rm -fv *.log
rm -rfv __pycache__
rm -rfv precice-profiling/
rm -fv .nfs*

echo "Cleaning restart_method"
cd restart_method/
sh clean_ruc.sh
rm -rfv ruc_*/
rm -rfv __pycache__