#!/bin/bash

for arg in "$@"; do
  if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
    Rscript /scripts/within-family-plots.R -h | sed "s/\.R//"
    exit 0
  fi
done

Rscript /scripts/within-family-plots.R $@ > within-family-plots.log 2>&1

if [ $? -ne 0 ]; then
  cat within-family-plots.log | sed "s/\.R//"
  exit 1
fi