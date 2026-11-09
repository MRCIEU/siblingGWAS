//! Rust port of `within-family.R` and `within-family-plots.R` from
//! https://github.com/MRCIEU/siblingGWAS, itself derived from https://github.com/LaurenceHowe/SiblingGWAS
//! (MIT-licensed; see NOTICE). Licensed under GPL-3.0-or-later; see LICENSE.
//!
//! For each variant, fits
//!
//! * model 1 (all individuals): `PHENOTYPE ~ GENOTYPE + covariates`, giving the population (TOTAL) effect;
//! * model 2 (siblings only): `PHENOTYPE ~ FAM_MEAN + CENTRED_GENOTYPE + covariates`, giving the
//!   between-family (BF) and within-family (WF) effects;
//!
//! with standard errors clustered on family, as `sandwich::vcovCL` + `lmtest::coeftest` do.

mod input;
mod model;
mod output;
mod plot;
mod tdist;

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use model::{Design, Layout, Terms};
use output::VariantResult;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser)]
#[command(
    name = "siblinggwas",
    version,
    disable_version_flag = true,
    about = "Estimate within-family (WF) and between-family (BF) effects of genetic variants using siblings"
)]
struct Args {
    /// The output of plink2 --recode A ("Sample-major additive (0/1/2) coding, suitable for loading from R...")
    #[arg(long)]
    raw: PathBuf,
    /// A plink extended variant information file
    #[arg(long)]
    bim: PathBuf,
    /// A plain text (separated value) file with three columns: FID IID "phenotype"
    #[arg(long)]
    phenotypes: PathBuf,
    /// A plain text (separated value) file of FID IID followed by one or more covariate columns
    #[arg(long)]
    covariates: PathBuf,
    /// A plain text (separated value) file with two columns: IID FID_FS, where FID_FS identifies dizygotic sibling families
    #[arg(long)]
    siblings: PathBuf,
    /// Output prefix; results are written to <OUT>_within-family-table.tsv
    #[arg(long)]
    out: String,
    /// Don't calculate p-values and standard errors
    #[arg(long = "no-se")]
    no_se: bool,
    /// Also plot population vs between- and within-family effects to <OUT>_pop_vs_family_beta.pdf
    #[arg(long)]
    plot: bool,
    /// Number of threads (default: all cores)
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// Print version
    #[arg(short = 'v', long, action = clap::ArgAction::Version)]
    version: Option<bool>,
    /// Print a summary of licence information
    #[arg(short = 'l', exclusive = true)]
    #[allow(dead_code)] // read from the matches in `run`, before `Args` is built
    licences_summary: bool,
    /// Print full licence information, including all third-party notices
    #[arg(long = "licences", alias = "licenses", exclusive = true)]
    #[allow(dead_code)] // read from the matches in `run`, before `Args` is built
    licences: bool,
}

/// The NOTICE file: this program's licence, and the notices of its third-party dependencies
const NOTICE: &str = include_str!("../NOTICE");

/// The opening statement of NOTICE, i.e. everything before the first third-party notice
fn get_notice_summary() -> &'static str {
    NOTICE
        .split_once("\n---")
        .map_or(NOTICE, |(summary, _)| summary)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {e:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    // The licence flags are checked before building `Args`, whose required fields they are exempt from.
    let matches = Args::command().get_matches();
    if matches.get_flag("licences") {
        print!("{NOTICE}");
        return Ok(());
    }
    if matches.get_flag("licences_summary") {
        println!("{}", get_notice_summary());
        println!("Use \"siblinggwas --licences\" to print the full text, including all third-party notices.");
        return Ok(());
    }
    let args = Args::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
    }
    let with_se = !args.no_se;
    let start = Instant::now();

    let mut ds = input::load(&input::Paths {
        raw: &args.raw,
        bim: &args.bim,
        phenotypes: &args.phenotypes,
        covariates: &args.covariates,
        siblings: &args.siblings,
    })?;
    let s = &ds.summary;
    eprintln!(
        "Read {} variants and {} individuals ({} with phenotype and covariates, {} complete) in {:.1}s",
        ds.bim.len(),
        s.raw_individuals,
        s.merged_individuals,
        ds.n,
        start.elapsed().as_secs_f64()
    );
    eprintln!("Phenotype: {}", ds.phenotype_name);
    eprintln!(
        "Covariates ({}): {}",
        ds.covariates.len(),
        describe_covariates(&ds)
    );
    eprintln!(
        "Siblings: {} in {} families",
        s.siblings, s.sibling_families
    );

    let layout = Layout::new(&ds);
    let all = Design::new(&ds, &layout, (0..ds.n as u32).collect());
    let sibs = Design::new(
        &ds,
        &layout,
        (0..ds.n as u32)
            .filter(|&i| ds.is_sibling[i as usize])
            .collect(),
    );
    eprintln!(
        "Model 1: {} individuals; model 2: {} individuals",
        all.n(),
        sibs.n()
    );
    ds.cov = Vec::new(); // now held, expanded and centred, by the designs

    let fit_start = Instant::now();
    let (results, skipped): (Vec<VariantResult>, Vec<Option<String>>) = (0..ds.bim.len())
        .into_par_iter()
        .map(|j| analyse_variant(ds.genotypes(j), &all, &sibs, with_se))
        .unzip();
    eprintln!(
        "Analysed {} variants in {:.1}s",
        ds.bim.len(),
        fit_start.elapsed().as_secs_f64()
    );
    report_skipped(&ds, &skipped);

    let path = PathBuf::from(format!("{}_within-family-table.tsv", args.out));
    output::write(&path, &ds.bim, &results, with_se)?;
    eprintln!("Wrote {}", path.display());
    if args.plot {
        let path = PathBuf::from(format!("{}_pop_vs_family_beta.pdf", args.out));
        plot::write(&path, &plot_title(&args.out), &results, with_se)?;
        eprintln!("Wrote {}", path.display());
    }
    eprintln!("Done in {:.1}s", start.elapsed().as_secs_f64());
    Ok(())
}

/// The plot's title: the last component of the output prefix (so `--out results/bmi` gives "bmi",
/// as running within-family-plots.R with `--out bmi` would).
fn plot_title(out: &str) -> String {
    Path::new(out)
        .file_name()
        .map_or_else(|| out.to_string(), |f| f.to_string_lossy().into_owned())
}

fn describe_covariates(ds: &input::Dataset) -> String {
    ds.covariates
        .iter()
        .map(|c| match &c.kind {
            input::CovKind::Numeric => c.name.clone(),
            input::CovKind::Logical => format!("{} (logical)", c.name),
            input::CovKind::Factor { levels } => {
                format!("{} (factor, {} levels)", c.name, levels.len())
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Prints one line per distinct reason a model was skipped, rather than one per variant.
fn report_skipped(ds: &input::Dataset, skipped: &[Option<String>]) {
    let mut by_reason: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for (b, reason) in ds.bim.iter().zip(skipped) {
        if let Some(reason) = reason {
            by_reason.entry(reason).or_default().push(&b.snp);
        }
    }
    for (reason, snps) in by_reason {
        let examples = snps.iter().take(5).copied().collect::<Vec<_>>().join(", ");
        let more = if snps.len() > 5 { ", ..." } else { "" };
        eprintln!(
            "Warning: {reason} for {} variant(s): {examples}{more}",
            snps.len()
        );
    }
}

/// The per-variant logic of `within-family.R`, including which outputs are left NA when.
/// Also returns why a model was skipped, if it was (`lm()` would have failed).
fn analyse_variant(
    geno: &[f64],
    all: &Design,
    sibs: &Design,
    with_se: bool,
) -> (VariantResult, Option<String>) {
    let mut r = VariantResult::default();

    let m1 = match all.fit(geno, Terms::Genotype, with_se) {
        Ok(fit) => fit,
        Err(e) => return (r, Some(format!("model 1 skipped ({e})"))),
    };
    r.beta[0] = m1.coef[0];
    r.beta[2] = m1.coef[1];
    if let Some(inf) = &m1.inference {
        let v = &inf.vcv;
        r.vcv1 = [v[0], v[1], v[3]];
        r.se[0] = inf.se[0];
        r.se[2] = inf.se[1];
        r.p[0] = inf.p[0];
        r.p[2] = inf.p[1];
    }

    let m2 = match sibs.fit(geno, Terms::FamilyMeanAndDeviation, with_se) {
        Ok(fit) => fit,
        Err(e) => return (r, Some(format!("model 2 skipped ({e})"))),
    };
    r.n_reg = Some(m2.n);
    r.beta[1] = m2.coef[0];
    r.beta[3] = m2.coef[1];
    r.beta[4] = m2.coef[2];
    if let Some(inf) = &m2.inference {
        let v = &inf.vcv;
        r.vcv2 = [v[0], v[1], v[2], v[4], v[5], v[8]];
        r.se[1] = inf.se[0];
        r.se[3] = inf.se[1];
        r.se[4] = inf.se[2];
        r.p[1] = inf.p[0];
        r.p[3] = inf.p[1];
        r.p[4] = inf.p[2];
    }
    (r, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_notice_summary_stops_before_first_third_party_notice() {
        let summary = get_notice_summary();
        assert!(
            summary.starts_with("siblingGWAS is licensed under the GNU General Public License")
        );
        assert!(!summary.contains("---"));
        assert!(summary.len() < NOTICE.len());
    }
}
