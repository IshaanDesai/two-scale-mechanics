#!/usr/bin/env bash
set -e -u

usage() { echo "Usage: cmd [-s] [-p n]" 1>&2; exit 1; }

# Check if no input argument was provided
if [ -z "$*" ] ; then
  echo "No input argument provided. Micro Manager is launched in serial"
  micro-manager-precice micro-manager-config-one-element.json
fi

while getopts ":sp" opt; do
  case ${opt} in
  s)
    micro-manager-precice micro-manager-config.json
    ;;
  p)
    mpiexec -n "$2" micro-manager-precice micro-manager-config-one-element.json
    ;;
  *)
    usage
    ;;
  esac
done
