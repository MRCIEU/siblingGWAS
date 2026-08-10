#!/bin/bash

for arg in "$@"; do
  if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
    Rscript /scripts/within-family.R -h | sed "s/\.R//"
    exit 0
  fi
done

Rscript /scripts/within-family.R $@ > within-family.log 2>&1

if [ $? -ne 0 ]; then
  cat within-family.log | sed "s/\.R//"
  exit 1
fi
