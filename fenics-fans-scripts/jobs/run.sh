#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    echo "Usage: $0 <jobname>"
    exit 1
}

extract() {
    local option="$1"
    local file="$2"
    grep "^#SBATCH.*$option" $file | sed -E "s/^#SBATCH[[:space:]]+$option[[:space:]]*=?//"
}

main() {
    # Check number of arguments
    if [[ $# -ne 1 ]]; then
        usage
    fi

    local jobname="$1"
    local jobscript="../jobs/${jobname}.sh"

    # Check if job script exists
    if [[ ! -f "$jobscript" ]]; then
        echo "Error: Job script '$jobscript' does not exist." >&2
        exit 1
    fi

    # Create output directory
    mkdir -p "../output/${jobname}"

    # Submit job
    #sbatch "$jobscript"
    #source "$jobscript"

    # run with salloc
    # we use salloc as a workaround because sub-job spawning does not work with sbatch
    licenses=$(extract "--licenses" "$jobscript")
    name=$(extract "-J" "$jobscript")
    work_dir=$(extract "-D" "$jobscript")
    mail_t=$(extract "--mail-type" "$jobscript")
    mail_u=$(extract "--mail-user" "$jobscript")
    time_limit=$(extract "--time" "$jobscript")
    account=$(extract "--account" "$jobscript")
    partition=$(extract "--partition" "$jobscript")
    nodes=$(extract "--nodes" "$jobscript")

    part="${PARTITION:-$partition}"
    allowed=(test micro general large fat)
    valid=false
    for v in "${allowed[@]}"; do
        if [[ "$part" == "$v" ]]; then
            valid=true
            break
        fi
    done
    if ! "$valid"; then
        echo "Error: invalid partition '$part'"
        exit 1
    fi
    time="${TIME:-$time_limit}"

    salloc_args="--licenses=${licenses} -J ${name} -D ${work_dir} --mail-type=${mail_t} --mail-user=${mail_u} --time=${time} --account=${account} --partition=${part} --nodes=${nodes} --exclusive ${jobscript}"
    echo "${salloc_args}"
    salloc_args=(
        --licenses="${licenses}"
        -J "${name}"
        -D "${work_dir}"
        --mail-type="${mail_t}"
        --mail-user="${mail_u}"
        --time="${time}"
        --account="${account}"
        --partition="${part}"
        --nodes="${nodes}"
        --exclusive
    )
    (salloc "${salloc_args[@]}" bash "${jobscript}") &> /dev/null &
    disown %%
}

main "$@"
