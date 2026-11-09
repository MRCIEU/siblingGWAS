#!/usr/bin/env bash
# Renders the --plot output of siblinggwas and, from the same results table, the plot made by
# validation/within-family-plots.R, with and without standard errors, as PDFs and as PNGs
# (via Ghostscript) for side-by-side comparison.
# Usage: compare_plots.sh <data dir from simulate.R ... clean> <path to within-family-plots.R> <output dir>
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
data="$(cd "$1" && pwd)"; rscript="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"
mkdir -p "$3"; out="$(cd "$3" && pwd)"
(cd "$here/.." && cargo build --release -q)
bin="$here/../target/release/siblinggwas"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

for flag in "" "--no-se"; do
  suffix="${flag:+_no-se}"
  "$bin" --raw "$data/geno.raw" --bim "$data/geno.bim" --phenotypes "$data/pheno.txt" \
    --covariates "$data/covariates.txt" --siblings "$data/siblings.txt" \
    --out "$tmp/bmi" --plot $flag 2> /dev/null
  mv "$tmp/bmi_pop_vs_family_beta.pdf" "$out/rust$suffix.pdf"
  cp "$tmp/bmi_within-family-table.tsv" "$out/table$suffix.tsv"
  (cd "$tmp" && Rscript "$rscript" --table "$out/table$suffix.tsv" --out bmi > /dev/null)
  mv "$tmp/bmi_pop_vs_family_beta.pdf" "$out/R$suffix.pdf"
  for impl in R rust; do
    gs -q -sDEVICE=png16m -r100 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
      -o "$out/$impl$suffix.png" "$out/$impl$suffix.pdf"
  done
done
echo "Wrote R/rust .pdf and .png (with and without SEs), and the tables they were drawn from, to $out"
