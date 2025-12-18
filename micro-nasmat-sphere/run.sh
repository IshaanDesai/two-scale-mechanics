#!/usr/bin/env bash
set -e -u

python3 -m venv --system-site-packages .venv
. .venv/bin/activate
pip install -r requirements.txt

cd pynasmat/
pip install .

cd ../

./run_micromanager_nasmat.sh
