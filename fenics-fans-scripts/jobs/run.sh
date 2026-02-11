#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

usage() {
    echo "Usage: $0 <jobname>"
    exit 1
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
    sbatch "$jobscript"
    #source "$jobscript"
}

main "$@"
