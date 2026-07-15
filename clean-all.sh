rm -rfv ./precice-profiling/
rm -rfv ./precice-run/
rm -rfv *.err
rm -rfv *.out

for dir in */; do
    if [ -f "${dir}clean.sh" ]; then
        echo "Cleaning ${dir}"
        (cd "${dir}" && ./clean.sh)
    fi
done
