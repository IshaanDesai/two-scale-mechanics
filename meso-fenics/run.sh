#!/usr/bin/env bash
set -e -u

# Parse command line arguments
PROBLEM=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -problem)
            PROBLEM="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 -problem <pseudo|bar|notch>"
            exit 1
            ;;
    esac
done

# Validate problem argument
if [[ -z "$PROBLEM" ]]; then
    echo "Error: -problem argument is required"
    echo "Usage: $0 -problem <pseudo|bar|notch>"
    exit 1
fi

# Set configuration file based on problem type
case $PROBLEM in
    pseudo)
        CONFIG_FILE="examples/config-pseudo.json"
        ;;
    bar)
        CONFIG_FILE="examples/coupled-bar/config-coupled.json"
        ;;
    notch)
        CONFIG_FILE="examples/coupled-notch/config-coupled.json"
        ;;
    *)
        echo "Error: Invalid problem type '$PROBLEM'"
        echo "Valid options: pseudo, bar, notch"
        exit 1
        ;;
esac

python3 -m venv --system-site-package .venv
. .venv/bin/activate
pip install -r requirements.txt && pip freeze > pip-installed-packages.log

python3 macro.py -problem $PROBLEM