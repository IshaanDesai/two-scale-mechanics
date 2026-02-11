#!/bin/bash

load_env() {
    echo "Loading the necessary modules"
    module load slurm_setup
    module load stack/24.5.0
    module load intel-toolkit
    module switch intel/2025.1.1 gcc/14.3.0
    module load boost/1.84.0-gcc14-impi
    module load eigen/3.4.0-gcc14
    module load fftw/3.3.10-gcc14-impi-openmp
    module load hdf5/1.14.5-gcc14-impi
    module load hypre/2.32.0-gcc14-impi
    module load scotch/7.0.4-gcc14-impi-i64
    module load parmetis/4.0.3-gcc14-impi-i64-r64
    module load autoconf/
    module load glib
    module load pkg-config/
    module load zlib
    module load python/3.10.12-extended
    export SLURM_EAR_LOAD_MPI_VERSION="intel"
    source ~/venv/py310/bin/activate
}

path::get_full() (
    cd "$1" && pwd
)

path::get_tsm() {
    path::get_full "../../.."
}

path::get_meso() {
    path::get_full "../../../meso-fenics"
}

edit::meso_input() {
    local file="$1"
    local meso_out="$2"
    local meso_state_out="$3"
    local precice_path="$4"

    sed -Ei "s#\"output_path\".*\$#\"output_path\": \"${meso_out}\",#g" "$file"
    sed -Ei "s#\"write_state\".*\$#\"write_smeso-statetate\": \"${meso_state_out}\",#g" "$file"
    sed -Ei "s#\"precice_xml_path\".*\$#\"precice_xml_path\": \"${precice_path}\",#g" "$file"
    sed -Ei "s#\"slurm_id\".*\$#\"slurm_id\": \"${SLURM_JOB_ID}\"#g" "$file"
}

edit::precice_input() {
    local file="$1"
    local log_out="$2"
    local prof_dir="$3"
    local export_dir="$4"
    local exch_dir="$5"

    sed -Ei "s#<sink.*\$#<sink type=\"file\" output=\"${log_out}\" filter=\"%Severity% > debug\" enabled=\"true\" />#g" "$file"
    sed -Ei "s#<!--profiling.*\$#<profiling directory=\"${prof_dir}\" />#g" "$file"
    sed -Ei "s#<export:vtu.*\$#<export:vtu directory=\"${export_dir}\" />#g" "$file"
    sed -Ei "s#<m2n:sockets.*\$#<m2n:sockets acceptor=\"Meso-structure\" connector=\"Micro-Manager\" exchange-directory=\"${exch_dir}\" network=\"ib0\" use-two-level-initialization=\"true\" />#g" "$file"
}

edit::mm_input() {
    local file="$1"
    local precice_path="$2"
    local num_mm="$3"
    local num_worker="$4"

    sed -Ei "s#precice_config_file_name\":.*\$#precice_config_file_name\": \"${precice_path}\",#g" "$file"
    sed -Ei ":a;N;\$!ba;s#\"decomposition\":[[:space:]]*\[[^]]*]#\"decomposition\": [${num_mm}, 1, 1]#g" "$file"
    sed -Ei "s#num_workers\":.*\$#num_workers\": ${num_worker},#g" "$file"
}

edit::fans_input() {
    local file="$1"
    local num_worker="$2"

    if [[ $num_worker -eq 0 ]]; then
        sed -Ei "s#results\":.*\$#results\": [],\n\"no_mpi\": true#g" "$file"
    fi
}