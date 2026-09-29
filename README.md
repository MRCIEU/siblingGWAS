# SiblingGWAS

Use information from genotyped dizygotic siblings to assess the robustness of GWAS analyses, by estimating within-family (WF) and between-family (BF) effects of genetic variants on continuous traits. This work is derived from [https://github.com/LaurenceHowe/SiblingGWAS](https://github.com/LaurenceHowe/SiblingGWAS).

## Use

Needed as input:

* Genotype data for all samples, consisting of a plink `raw` file of genotypes + a plink `bim` file of variant information (we input LD-clumped, significantly associated variants only)

* A value-separated file containing information for a single phenotype, one row for each sample in the entire dataset

* A value-separated file with covariate information, one row for each sample in the entire dataset

* A value-separated file of sibling relationships, with one column consisting of sample IDs, and one column consisting of family IDs which define sibling relationships. One row for each sibling in the dataset

There is an optional `--no-se` flag, which will force scripts to not calculate standard errors or p-values for the estimated effects, which will significantly speed up the calculations.

### Rscript

You can run things using the R scripts (`within-family.R` and `within-family-plots.R`) provided:

```
❯ Rscript within-family.R -h
Usage: within-family.R --raw <raw file> --bim <bim file> --phenotypes <phenotype file> --covariates <covariate file> --siblings <sibling file> --out <output prefix> [--no-se]
...
```

If you provide the `--no-se` flag, the output table will contain fewer columns, but can still be parsed with `within-family-plots.R` to generate the plots (minus SEs).

#### Arguments

* `<raw file>` - The output of `plink2 --recode A` ("Sample-major additive (0/1/2) coding, suitable for loading from R...")

* `<bim file>` - [A plink extended variant information file](https://plink.readthedocs.io/en/latest/plink_fmt/#bim)

* `<phenotype file>` - A plain text (separated value) file consisting of EXACTLY three columns: `FID` `IID` `"phenotype"` where "phenotype" can be named sensibly

* `<covariate file>` - A plain text (separated value) file consisting of AT LEAST three columns. First: `FID` `IID` + at least one column containing covariates: `<"covariate 1">` `<"covariate 2">` `...`, etc.

* `<sibling file>` - A plain text (separated value) file containing sibling information with exactly two columns: `IID` `FID_FS` where `FID_FS` is a family identifier for dizygotic siblings (see [the original method](https://github.com/LaurenceHowe/SiblingGWAS/tree/master))

* `<output prefix>` - A string that will be prefixed to the output file according to: `paste0(<output prefix>, "_within-family-table.tsv")`

This will generate `<output prefix>_within-family-table.tsv`, which you can in turn feed to `within-family-plots.R`:

```
❯ Rscript within-family-plots.R -h
Usage: within-family-plots.R --table <input file> --out <output prefix>
...
```

where:

* `<input file>` - The output of `within-family.R`

* `<output prefix>` - A string that will be prefixed to the output file according to: `paste0(<output prefix>, "_pop_vs_family_beta.pdf")`

The plots look like:

<p align="left">
  <img src=".github/plot.png" alt="" width="738">
</p>

### Docker

The scripts above can also be run using the Docker image created from the `Dockerfile` in this repository, e.g.

```
docker build --platform linux/x86_64 --no-cache -t mrcieu/siblinggwas .
```

The Rscripts above are available in the `$PATH` as `within-family` and `within-family-plots` (don't include the `.R` extension):

```
❯ docker run mrcieu/siblinggwas:latest within-family -h

❯ docker run mrcieu/siblinggwas:latest within-family-plots -h
```
