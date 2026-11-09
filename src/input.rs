//! Reading the input files and assembling the analysis sample.
//!
//! Mirrors the `fread`/`merge` steps of `within-family.R`: phenotypes and covariates are inner-joined
//! to the plink `.raw` file on (FID, IID), siblings are left-joined on IID, and individuals with a
//! missing phenotype or covariate are dropped (they would be removed by `na.omit` for every variant).

use anyhow::{bail, ensure, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Copy)]
enum Sep {
    Tab,
    Comma,
    Whitespace,
}

impl Sep {
    fn detect(header: &str) -> Sep {
        if header.contains('\t') {
            Sep::Tab
        } else if header.contains(',') {
            Sep::Comma
        } else {
            Sep::Whitespace
        }
    }

    fn split(self, line: &str) -> Vec<&str> {
        match self {
            Sep::Tab => line.split('\t').map(clean).collect(),
            Sep::Comma => line.split(',').map(clean).collect(),
            Sep::Whitespace => line.split_ascii_whitespace().map(clean).collect(),
        }
    }
}

fn clean(field: &str) -> &str {
    let f = field.trim();
    f.strip_prefix('"')
        .and_then(|f| f.strip_suffix('"'))
        .unwrap_or(f)
}

fn is_na(field: &str) -> bool {
    field.is_empty() || field == "NA"
}

fn open(path: &Path) -> Result<BufReader<File>> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    Ok(BufReader::with_capacity(1 << 20, file))
}

/// Reads a delimited file with a header line, calling `f(line_number, fields)` for every non-blank
/// data line. The delimiter (tab, comma or whitespace) is detected from the header. Returns the header.
fn scan<F>(path: &Path, mut f: F) -> Result<Vec<String>>
where
    F: FnMut(usize, &[&str]) -> Result<()>,
{
    let mut reader = open(path)?;
    let mut line = String::new();
    let mut lineno = 0;
    let header = loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            bail!("{} is empty", path.display());
        }
        lineno += 1;
        if !line.trim().is_empty() {
            break line.trim_end_matches(['\n', '\r']).to_string();
        }
    };
    let sep = Sep::detect(&header);
    let header: Vec<String> = sep.split(&header).into_iter().map(String::from).collect();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        lineno += 1;
        let l = line.trim_end_matches(['\n', '\r']);
        if l.trim().is_empty() {
            continue;
        }
        let fields = sep.split(l);
        ensure!(
            fields.len() == header.len(),
            "{} line {}: expected {} fields (as in the header) but found {}",
            path.display(),
            lineno,
            header.len(),
            fields.len()
        );
        f(lineno, &fields).with_context(|| format!("{} line {}", path.display(), lineno))?;
    }
    Ok(header)
}

/// `fread` reads integer-like ID columns as integers, so "007" and "7" refer to the same individual.
/// Canonicalise such columns the same way; columns containing any non-integer ID are left as-is.
fn normalise_ids(ids: &mut [String]) {
    let int_like = |s: &str| {
        let digits = s.strip_prefix(['-', '+']).unwrap_or(s);
        !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
    };
    if !ids.iter().all(|s| int_like(s)) {
        return;
    }
    for s in ids.iter_mut() {
        let neg = s.starts_with('-');
        let digits = s.trim_start_matches(['-', '+']).trim_start_matches('0');
        *s = match (digits.is_empty(), neg) {
            (true, _) => "0".to_string(),
            (false, true) => format!("-{digits}"),
            (false, false) => digits.to_string(),
        };
    }
}

pub struct BimRow {
    pub chr: String,
    pub snp: String,
    pub bp: String,
    pub a1: String,
    pub a2: String,
}

fn read_bim(path: &Path) -> Result<Vec<BimRow>> {
    let mut out = Vec::new();
    for (i, line) in open(path)?.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_ascii_whitespace().collect();
        ensure!(
            f.len() == 6,
            "{} line {}: expected 6 columns, found {}",
            path.display(),
            i + 1,
            f.len()
        );
        out.push(BimRow {
            chr: f[0].to_string(),
            snp: f[1].to_string(),
            bp: f[3].to_string(),
            a1: f[4].to_string(),
            a2: f[5].to_string(),
        });
    }
    Ok(out)
}

struct Phenotypes {
    fid: Vec<String>,
    iid: Vec<String>,
    value: Vec<f64>,
}

fn read_phenotypes(path: &Path) -> Result<Phenotypes> {
    let mut p = Phenotypes {
        fid: Vec::new(),
        iid: Vec::new(),
        value: Vec::new(),
    };
    let header = scan(path, |_, f| {
        ensure!(
            f.len() >= 3,
            "expected at least 3 columns: FID IID <phenotype>"
        );
        let v = if is_na(f[2]) {
            f64::NAN
        } else {
            f[2].parse()
                .with_context(|| format!("phenotype value '{}' is not numeric", f[2]))?
        };
        p.fid.push(f[0].to_string());
        p.iid.push(f[1].to_string());
        p.value.push(v);
        Ok(())
    })?;
    ensure!(
        header.len() >= 3,
        "{}: expected at least 3 columns: FID IID <phenotype>",
        path.display()
    );
    normalise_ids(&mut p.fid);
    normalise_ids(&mut p.iid);
    Ok(p)
}

#[derive(Debug, Clone)]
pub enum CovKind {
    /// Used as-is.
    Numeric,
    /// TRUE/FALSE (as recognised by `fread`): a single 0/1 column, like `lm` does for logicals.
    Logical,
    /// Anything else: a factor with levels in sorted order, expanded to treatment-contrast dummies.
    Factor { levels: Vec<String> },
}

pub struct Covariate {
    pub name: String,
    pub kind: CovKind,
}

struct Covariates {
    fid: Vec<String>,
    iid: Vec<String>,
    covariates: Vec<Covariate>,
    /// Row-major (row × covariate). Factor values are level indices. NaN is missing.
    values: Vec<f64>,
}

fn parse_logical(s: &str) -> Option<f64> {
    match s {
        "TRUE" | "True" | "true" => Some(1.0),
        "FALSE" | "False" | "false" => Some(0.0),
        _ => None,
    }
}

fn read_covariates(path: &Path) -> Result<Covariates> {
    // Pass 1: IDs, and the type of each covariate column.
    let mut fid = Vec::new();
    let mut iid = Vec::new();
    let mut numeric: Vec<bool> = Vec::new();
    let mut logical: Vec<bool> = Vec::new();
    let header = scan(path, |_, f| {
        if numeric.is_empty() {
            numeric = vec![true; f.len().saturating_sub(2)];
            logical = vec![true; f.len().saturating_sub(2)];
        }
        fid.push(f[0].to_string());
        iid.push(f[1].to_string());
        for (j, v) in f[2..].iter().enumerate() {
            if is_na(v) {
                continue;
            }
            if numeric[j] && v.parse::<f64>().is_err() {
                numeric[j] = false;
            }
            if logical[j] && parse_logical(v).is_none() {
                logical[j] = false;
            }
        }
        Ok(())
    })?;
    ensure!(
        header.len() >= 2,
        "{}: expected at least 2 columns: FID IID [covariates...]",
        path.display()
    );
    let ncov = header.len() - 2;
    numeric.resize(ncov, true);
    logical.resize(ncov, true);

    // Pass 2: values.
    let nrow = fid.len();
    let mut values = vec![f64::NAN; nrow * ncov];
    let mut level_ids: Vec<HashMap<String, usize>> = vec![HashMap::new(); ncov];
    let mut row = 0;
    scan(path, |_, f| {
        for (j, v) in f[2..].iter().enumerate() {
            if is_na(v) {
                continue;
            }
            values[row * ncov + j] = if numeric[j] {
                v.parse().unwrap()
            } else if logical[j] {
                parse_logical(v).unwrap()
            } else {
                let next = level_ids[j].len();
                *level_ids[j].entry(v.to_string()).or_insert(next) as f64
            };
        }
        row += 1;
        Ok(())
    })?;

    let mut covariates = Vec::with_capacity(ncov);
    for j in 0..ncov {
        let kind = if numeric[j] {
            CovKind::Numeric
        } else if logical[j] {
            CovKind::Logical
        } else {
            // Re-number the levels in sorted order so that the first level is the reference.
            let mut levels: Vec<(String, usize)> = level_ids[j].drain().collect();
            levels.sort();
            let mut remap = vec![0.0; levels.len()];
            for (sorted, (_, first_seen)) in levels.iter().enumerate() {
                remap[*first_seen] = sorted as f64;
            }
            for r in 0..nrow {
                let v = &mut values[r * ncov + j];
                if !v.is_nan() {
                    *v = remap[*v as usize];
                }
            }
            CovKind::Factor {
                levels: levels.into_iter().map(|(l, _)| l).collect(),
            }
        };
        covariates.push(Covariate {
            name: header[j + 2].clone(),
            kind,
        });
    }
    normalise_ids(&mut fid);
    normalise_ids(&mut iid);
    Ok(Covariates {
        fid,
        iid,
        covariates,
        values,
    })
}

/// IID -> FID_FS for individuals with a non-missing sibling family ID.
fn read_siblings(path: &Path) -> Result<HashMap<String, String>> {
    let mut iid = Vec::new();
    let mut fid_fs = Vec::new();
    let header = scan(path, |_, f| {
        ensure!(f.len() >= 2, "expected 2 columns: IID FID_FS");
        iid.push(f[0].to_string());
        fid_fs.push(f[1].to_string());
        Ok(())
    })?;
    ensure!(
        header.len() >= 2,
        "{}: expected 2 columns: IID FID_FS",
        path.display()
    );
    normalise_ids(&mut iid);
    let non_missing: Vec<usize> = (0..fid_fs.len()).filter(|&i| !is_na(&fid_fs[i])).collect();
    let mut fs: Vec<String> = non_missing.iter().map(|&i| fid_fs[i].clone()).collect();
    normalise_ids(&mut fs);
    let mut map = HashMap::with_capacity(fs.len());
    for (&i, f) in non_missing.iter().zip(fs) {
        if map.insert(iid[i].clone(), f).is_some() {
            bail!("{}: IID {} appears more than once", path.display(), iid[i]);
        }
    }
    Ok(map)
}

pub struct Paths<'a> {
    pub raw: &'a Path,
    pub bim: &'a Path,
    pub phenotypes: &'a Path,
    pub covariates: &'a Path,
    pub siblings: &'a Path,
}

pub struct Summary {
    pub raw_individuals: usize,
    pub merged_individuals: usize,
    pub siblings: usize,
    pub sibling_families: usize,
}

/// The analysis sample: individuals in the `.raw` file with phenotype, all covariates and a cluster ID.
pub struct Dataset {
    pub bim: Vec<BimRow>,
    pub phenotype_name: String,
    pub covariates: Vec<Covariate>,
    /// Number of individuals.
    pub n: usize,
    pub y: Vec<f64>,
    /// Row-major (individual × covariate); factor values are level indices.
    pub cov: Vec<f64>,
    /// Family cluster: FID_FS for siblings, FID otherwise.
    pub cluster: Vec<u32>,
    pub is_sibling: Vec<bool>,
    /// Column-major (variant × individual); NaN is a missing genotype.
    pub geno: Vec<f64>,
    pub summary: Summary,
}

impl Dataset {
    pub fn cov_row(&self, i: usize) -> &[f64] {
        let k = self.covariates.len();
        &self.cov[i * k..(i + 1) * k]
    }

    pub fn genotypes(&self, variant: usize) -> &[f64] {
        &self.geno[variant * self.n..(variant + 1) * self.n]
    }
}

fn parse_genotype(s: &str) -> Result<f64> {
    Ok(match s {
        "0" => 0.0,
        "1" => 1.0,
        "2" => 2.0,
        "NA" | "" => f64::NAN,
        _ => s
            .parse()
            .with_context(|| format!("genotype value '{s}' is not numeric"))?,
    })
}

pub fn load(paths: &Paths) -> Result<Dataset> {
    let bim = read_bim(paths.bim)?;
    let phen = read_phenotypes(paths.phenotypes)?;
    let phenotype_name = scan_header(paths.phenotypes)?
        .get(2)
        .cloned()
        .unwrap_or_default();
    let cov = read_covariates(paths.covariates)?;
    let sibs = read_siblings(paths.siblings)?;

    let phen_idx =
        index_by_key(&phen.fid, &phen.iid).context(paths.phenotypes.display().to_string())?;
    let cov_idx =
        index_by_key(&cov.fid, &cov.iid).context(paths.covariates.display().to_string())?;

    // Pass 1 over the .raw file: IDs only.
    let mut raw_fid = Vec::new();
    let mut raw_iid = Vec::new();
    let raw_header = scan(paths.raw, |_, f| {
        raw_fid.push(f[0].to_string());
        raw_iid.push(f[1].to_string());
        Ok(())
    })?;
    ensure!(
        raw_header.len() >= 6,
        "{}: expected FID IID PAT MAT SEX PHENOTYPE <variants...>",
        paths.raw.display()
    );
    let n_variants = raw_header.len() - 6;
    ensure!(
        n_variants == bim.len(),
        "{} has {} variant columns but {} has {} variants",
        paths.raw.display(),
        n_variants,
        paths.bim.display(),
        bim.len()
    );
    let mismatched = bim
        .iter()
        .zip(&raw_header[6..])
        .filter(|(b, col)| {
            !col.strip_prefix(b.snp.as_str())
                .is_some_and(|rest| rest.starts_with('_'))
        })
        .count();
    if mismatched > 0 {
        eprintln!(
            "Warning: {mismatched} variant column name(s) in the .raw file do not match the .bim SNP IDs \
             (expected <SNP>_<allele>); variants are matched by position"
        );
    }
    normalise_ids(&mut raw_fid);
    normalise_ids(&mut raw_iid);
    index_by_key(&raw_fid, &raw_iid).context(paths.raw.display().to_string())?;

    // Merge, and keep individuals with complete phenotype and covariate data.
    let ncov = cov.covariates.len();
    let mut keep: Vec<Option<usize>> = vec![None; raw_fid.len()];
    let mut y = Vec::new();
    let mut cov_values = Vec::new();
    let mut cluster_ids: HashMap<String, u32> = HashMap::new();
    let mut cluster = Vec::new();
    let mut is_sibling = Vec::new();
    let mut merged = 0;
    for r in 0..raw_fid.len() {
        let key = (raw_fid[r].as_str(), raw_iid[r].as_str());
        let (Some(&pi), Some(&ci)) = (phen_idx.get(&key), cov_idx.get(&key)) else {
            continue;
        };
        merged += 1;
        let row = &cov.values[ci * ncov..(ci + 1) * ncov];
        if phen.value[pi].is_nan() || row.iter().any(|v| v.is_nan()) {
            continue;
        }
        let (family, sib) = match sibs.get(&raw_iid[r]) {
            Some(fs) => (fs.as_str(), true),
            None => (raw_fid[r].as_str(), false),
        };
        let next = cluster_ids.len() as u32;
        cluster.push(*cluster_ids.entry(family.to_string()).or_insert(next));
        keep[r] = Some(y.len());
        y.push(phen.value[pi]);
        cov_values.extend_from_slice(row);
        is_sibling.push(sib);
    }
    let n = y.len();
    let sibling_families = {
        let mut fams: Vec<u32> = (0..n)
            .filter(|&i| is_sibling[i])
            .map(|i| cluster[i])
            .collect();
        fams.sort_unstable();
        fams.dedup();
        fams.len()
    };

    // Pass 2 over the .raw file: genotypes of the retained individuals, stored per variant. Rows are
    // parsed in blocks and transposed so that writes into the column-major matrix are contiguous.
    const BLOCK: usize = 256;
    let mut geno = vec![f64::NAN; n * n_variants];
    let mut block = Vec::with_capacity(BLOCK * n_variants);
    let mut block_start = 0;
    let flush = |block: &mut Vec<f64>, start: usize, geno: &mut Vec<f64>| {
        let rows = block.len() / n_variants.max(1);
        for j in 0..n_variants {
            let col = &mut geno[j * n + start..j * n + start + rows];
            for (b, g) in col.iter_mut().enumerate() {
                *g = block[b * n_variants + j];
            }
        }
        block.clear();
    };
    let mut r = 0;
    scan(paths.raw, |_, f| {
        if let Some(i) = keep[r] {
            if block.is_empty() {
                block_start = i;
            }
            for v in &f[6..] {
                block.push(parse_genotype(v)?);
            }
            if block.len() == BLOCK * n_variants {
                flush(&mut block, block_start, &mut geno);
            }
        }
        r += 1;
        Ok(())
    })?;
    if !block.is_empty() {
        flush(&mut block, block_start, &mut geno);
    }

    let siblings = is_sibling.iter().filter(|&&s| s).count();
    Ok(Dataset {
        bim,
        phenotype_name,
        covariates: cov.covariates,
        n,
        y,
        cov: cov_values,
        cluster,
        is_sibling,
        geno,
        summary: Summary {
            raw_individuals: raw_fid.len(),
            merged_individuals: merged,
            siblings,
            sibling_families,
        },
    })
}

fn scan_header(path: &Path) -> Result<Vec<String>> {
    let mut reader = open(path)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let l = line.trim_end_matches(['\n', '\r']);
        if !l.trim().is_empty() {
            return Ok(Sep::detect(l)
                .split(l)
                .into_iter()
                .map(String::from)
                .collect());
        }
        line.clear();
    }
    bail!("{} is empty", path.display())
}

fn index_by_key<'a>(
    fid: &'a [String],
    iid: &'a [String],
) -> Result<HashMap<(&'a str, &'a str), usize>> {
    let mut map = HashMap::with_capacity(fid.len());
    for (i, (f, d)) in fid.iter().zip(iid).enumerate() {
        if map.insert((f.as_str(), d.as_str()), i).is_some() {
            bail!("individual FID={f} IID={d} appears more than once");
        }
    }
    Ok(map)
}
