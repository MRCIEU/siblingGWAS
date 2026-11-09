#!/usr/bin/env bash
# Runs within-family.R and the Rust binary on a simulated dataset and compares the outputs.
# Usage: run_comparison.sh <data dir> <path to within-family.R> [extra args, e.g. --no-se]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
data="$1"; rscript="$2"; shift 2
bin="$here/../target/release/siblinggwas"
args=(--raw "$data/geno.raw" --bim "$data/geno.bim" --phenotypes "$data/pheno.txt"
      --covariates "$data/covariates.txt" --siblings "$data/siblings.txt")
suffix="${1:+_${1#--}}"

echo "== R"
start=$(date +%s.%N)
Rscript "$rscript" "${args[@]}" --out "$data/r$suffix" "$@" > "$data/r$suffix.log" 2>&1 || { cat "$data/r$suffix.log"; exit 1; }
echo "$(echo "$(date +%s.%N) - $start" | bc) s"
grep -v -E "^(Attaching|The following|    |$)" "$data/r$suffix.log" | sed 's/^/  /' | head -20 || true

echo "== Rust"
start=$(date +%s.%N)
"$bin" "${args[@]}" --out "$data/rs$suffix" "$@"
echo "$(echo "$(date +%s.%N) - $start" | bc) s"

echo "== Compare"
Rscript "$here/compare.R" "$data/r${suffix}_within-family-table.tsv" "$data/rs${suffix}_within-family-table.tsv"
