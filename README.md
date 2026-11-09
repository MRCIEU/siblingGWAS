# siblingGWAS-rs

A fast Rust port of `within-family.R` and `within-family-plots.R` from [siblingGWAS](https://github.com/MRCIEU/siblingGWAS). It estimates population, between-family and within-family effects for each variant, with standard errors clustered on family, and doesn't need R.

It takes the same inputs and writes the same `<prefix>_within-family-table.tsv` as the R script, and agrees with it to within about 1e-11. With `--plot` it also draws the plot that `within-family-plots.R` makes.

This is a work in progress, so the interface may change.

## Installation

With Cargo:

```
cargo install --git https://github.com/MRCIEU/siblingGWAS-rs
```

Or clone this repository and build it:

```
cargo build --release
./target/release/siblinggwas --help
```

The plot is drawn with [rs-graphics](https://github.com/MRCIEU/rs-graphics), which Cargo fetches over SSH, so you need access to that repository.

## Usage

```
siblinggwas --raw <raw file> --bim <bim file> --phenotypes <phenotype file> \
    --covariates <covariate file> --siblings <sibling file> --out <output prefix> \
    [--no-se] [--plot] [--threads N]
```

| Option | Description |
| --- | --- |
| `--raw` | Genotypes from `plink2 --recode A` |
| `--bim` | The matching plink `.bim` file |
| `--phenotypes` | FID, IID and phenotype |
| `--covariates` | FID, IID and one or more covariates |
| `--siblings` | IID and FID_FS, which identifies sibling families |
| `--out` | Output prefix |
| `--no-se` | Don't calculate standard errors and p-values |
| `--plot` | Also write `<prefix>_pop_vs_family_beta.pdf` |
| `--threads` | Number of threads (default: all cores) |

The arguments are the same as for `within-family.R`; see the [siblingGWAS README](https://github.com/MRCIEU/siblingGWAS) for details of the input files.

## Licence

GPL-3.0-or-later; see [LICENSE](LICENSE). Third-party notices are in [NOTICE](NOTICE). `siblinggwas -l` prints a summary and `siblinggwas --licences` the full text.
