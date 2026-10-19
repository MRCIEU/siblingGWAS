#!/usr/bin/env Rscript
# Compares two *_within-family-table.tsv files cell by cell.
# Usage: Rscript compare.R <R output> <Rust output> [relative tolerance, default 1e-8]
# Exits non-zero if the columns, variants or NA patterns differ, or any value differs by more than the tolerance.

suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
a <- fread(args[1], colClasses = "character", na.strings = "")
b <- fread(args[2], colClasses = "character", na.strings = "")
tol <- if (length(args) >= 3) as.numeric(args[3]) else 1e-8

ok <- TRUE
fail <- function(...) { cat("FAIL:", sprintf(...), "\n"); ok <<- FALSE }

if (!identical(names(a), names(b))) fail("column names differ")
if (nrow(a) != nrow(b)) fail("row counts differ: %d vs %d", nrow(a), nrow(b))
for (col in c("CHR", "SNP", "BP", "A1", "A2", "N_REG")) {
  if (!identical(a[[col]], b[[col]])) fail("column %s differs", col)
}
num <- setdiff(names(a), c("CHR", "SNP", "BP", "A1", "A2", "N_REG"))
worst <- data.table(column = num, max_rel_diff = NA_real_, n_values = NA_integer_)
for (col in num) {
  x <- as.numeric(a[[col]]); y <- as.numeric(b[[col]])
  if (!identical(is.na(x), is.na(y))) {
    fail("NA pattern differs in %s for variants: %s", col, paste(a$SNP[is.na(x) != is.na(y)], collapse = ","))
    next
  }
  keep <- !is.na(x)
  rel <- abs(x[keep] - y[keep]) / pmax(abs(x[keep]), abs(y[keep]), .Machine$double.xmin)
  rel[x[keep] == y[keep]] <- 0
  worst[column == col, `:=`(max_rel_diff = if (any(keep)) max(rel) else 0, n_values = sum(keep))]
  if (any(rel > tol)) fail("%s differs by up to %.3g (relative) for %s", col, max(rel), paste(a$SNP[keep][rel > tol], collapse = ","))
}
print(worst, nrows = 100)
if (ok) cat(sprintf("OK: %d variants, all %d numeric columns agree to %.0e relative, identical NA patterns\n", nrow(a), length(num), tol))
quit(status = if (ok) 0 else 1)
