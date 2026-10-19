//! Per-variant OLS with family-clustered sandwich standard errors.
//!
//! Reproduces, for each variant,
//!
//! ```r
//! model <- lm(PHENOTYPE ~ <genotype terms> + <covariates>, data)
//! vcv   <- sandwich::vcovCL(model, cluster = data$FID)   # type "HC1", cadjust = TRUE
//! lmtest::coeftest(model, vcov. = vcv)
//! ```
//!
//! without refitting the covariates for every variant. The cross-product matrix of
//! `[1, covariates, phenotype]` is computed once per model; for each variant the rows with a missing
//! genotype are subtracted from it and the genotype columns are added, which costs O(n·p) instead of
//! the O(n·p²) of a QR decomposition. Columns are centred (and the intercept transformed back) to keep
//! the normal equations well conditioned.
//!
//! The clustered covariance is `V = c · (X'X)⁻¹ (Σ_g u_g u_g') (X'X)⁻¹` with `u_g = Σ_{i∈g} e_i x_i` and
//! `c = G/(G-1) · (n-1)/(n-k)`, exactly as in `sandwich::meatCL` + `sandwich::sandwich`. Only the rows of
//! `(X'X)⁻¹` belonging to the intercept and genotype terms are needed, so each observation contributes
//! a handful of dot products.

use crate::input::{CovKind, Dataset};
use crate::tdist;
use rayon::prelude::*;
use std::fmt;

/// `lm.fit`'s default tolerance for detecting linearly dependent columns (LINPACK `dqrdc2`).
const QR_TOL: f64 = 1e-7;

/// The genotype terms that follow the intercept in the model formula.
#[derive(Clone, Copy)]
pub enum Terms {
    /// `PHENOTYPE ~ GENOTYPE + covariates`
    Genotype,
    /// `PHENOTYPE ~ FAM_MEAN + CENTRED_GENOTYPE + covariates`, where `FAM_MEAN` is the mean genotype of
    /// the individual's family (among individuals in the regression) and `CENTRED_GENOTYPE` the deviation from it.
    FamilyMeanAndDeviation,
}

impl Terms {
    fn count(self) -> usize {
        match self {
            Terms::Genotype => 1,
            Terms::FamilyMeanAndDeviation => 2,
        }
    }
}

/// Errors under which `lm()` itself fails, and the R script skips the rest of the variant.
#[derive(Debug)]
pub enum FitError {
    NoObservations,
    SingleLevelFactor(String),
}

impl fmt::Display for FitError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FitError::NoObservations => write!(f, "0 (non-NA) cases"),
            FitError::SingleLevelFactor(name) => write!(
                f,
                "contrasts can be applied only to factors with 2 or more levels (covariate '{name}')"
            ),
        }
    }
}

pub struct Fit {
    pub n: usize,
    /// Intercept followed by the genotype terms; NaN for aliased coefficients.
    pub coef: Vec<f64>,
    /// Present when standard errors were requested and none of `coef` is aliased.
    pub inference: Option<Inference>,
}

pub struct Inference {
    /// Clustered covariance of `coef`, row-major.
    pub vcv: Vec<f64>,
    pub se: Vec<f64>,
    pub p: Vec<f64>,
}

/// A factor covariate's dummy columns: one per level, the reference level is chosen per fit.
struct FactorCols {
    name: String,
    start: usize,
    levels: usize,
}

/// Layout of the expanded covariate columns.
pub struct Layout {
    q: usize,
    starts: Vec<usize>,
    factors: Vec<FactorCols>,
}

impl Layout {
    pub fn new(ds: &Dataset) -> Layout {
        let mut q = 0;
        let mut starts = Vec::new();
        let mut factors = Vec::new();
        for c in &ds.covariates {
            starts.push(q);
            match &c.kind {
                CovKind::Numeric | CovKind::Logical => q += 1,
                CovKind::Factor { levels } => {
                    factors.push(FactorCols {
                        name: c.name.clone(),
                        start: q,
                        levels: levels.len(),
                    });
                    q += levels.len();
                }
            }
        }
        Layout { q, starts, factors }
    }

    fn fill(&self, ds: &Dataset, i: usize, out: &mut [f64]) {
        for (c, (&v, &s)) in ds.cov_row(i).iter().zip(&self.starts).enumerate() {
            match ds.covariates[c].kind {
                CovKind::Numeric | CovKind::Logical => out[s] = v,
                CovKind::Factor { ref levels } => {
                    out[s..s + levels.len()].fill(0.0);
                    out[s + v as usize] = 1.0;
                }
            }
        }
    }
}

/// The variant-independent part of one regression model.
pub struct Design<'a> {
    layout: &'a Layout,
    /// Individuals (indices into the dataset) in this model.
    rows: Vec<u32>,
    /// Covariates, row-major and centred by `zmean`.
    z: Vec<f64>,
    zmean: Vec<f64>,
    /// Phenotype, centred by `ymean`.
    y: Vec<f64>,
    ymean: f64,
    cluster: Vec<u32>,
    n_clusters: usize,
    /// Upper triangle of the cross-product matrix of `[1, z, y]`, (q+2)×(q+2) row-major.
    gram: Vec<f64>,
    /// Number of individuals with each factor level (indexed by covariate column; 0 elsewhere).
    level_counts: Vec<usize>,
}

impl<'a> Design<'a> {
    pub fn new(ds: &Dataset, layout: &'a Layout, rows: Vec<u32>) -> Design<'a> {
        let (n, q) = (rows.len(), layout.q);
        let mut z = vec![0.0; n * q];
        z.par_chunks_mut(q.max(1))
            .zip(&rows)
            .for_each(|(zi, &r)| layout.fill(ds, r as usize, zi));
        let mut y: Vec<f64> = rows.iter().map(|&r| ds.y[r as usize]).collect();

        let mut level_counts = vec![0; q];
        for f in &layout.factors {
            for i in 0..n {
                for l in 0..f.levels {
                    level_counts[f.start + l] += (z[i * q + f.start + l] == 1.0) as usize;
                }
            }
        }

        let mut zmean = vec![0.0; q];
        for zi in z.chunks_exact(q.max(1)).take(n) {
            for (m, v) in zmean.iter_mut().zip(zi) {
                *m += v;
            }
        }
        zmean.iter_mut().for_each(|m| *m /= n.max(1) as f64);
        z.par_chunks_mut(q.max(1))
            .for_each(|zi| zi.iter_mut().zip(&zmean).for_each(|(v, m)| *v -= m));
        let ymean = y.iter().sum::<f64>() / n.max(1) as f64;
        y.iter_mut().for_each(|v| *v -= ymean);

        let mut ids = vec![u32::MAX; ds.n];
        let mut n_clusters = 0;
        let cluster = rows
            .iter()
            .map(|&r| {
                let c = &mut ids[ds.cluster[r as usize] as usize];
                if *c == u32::MAX {
                    *c = n_clusters;
                    n_clusters += 1;
                }
                *c
            })
            .collect();

        let gram = gram(&z, &y, q);
        Design {
            layout,
            rows,
            z,
            zmean,
            y,
            ymean,
            cluster,
            n_clusters: n_clusters as usize,
            gram,
            level_counts,
        }
    }

    pub fn n(&self) -> usize {
        self.rows.len()
    }

    /// Fits the model for one variant. `geno` holds the genotypes of all individuals in the dataset.
    pub fn fit(&self, geno: &[f64], terms: Terms, with_se: bool) -> Result<Fit, FitError> {
        let (nr, q) = (self.rows.len(), self.layout.q);
        let ng = terms.count();
        let qw = q + 2; // [1, z, y]
        let g: Vec<f64> = self.rows.iter().map(|&r| geno[r as usize]).collect();
        let n = g.iter().filter(|v| !v.is_nan()).count();
        if n == 0 {
            return Err(FitError::NoObservations);
        }

        // Genotype terms, row-major nr×ng (rows with a missing genotype are ignored).
        let mut x = vec![f64::NAN; nr * ng];
        match terms {
            Terms::Genotype => x.copy_from_slice(&g),
            Terms::FamilyMeanAndDeviation => {
                let mut sum = vec![0.0; self.n_clusters];
                let mut count = vec![0u32; self.n_clusters];
                for (i, &gi) in g.iter().enumerate().filter(|(_, v)| !v.is_nan()) {
                    sum[self.cluster[i] as usize] += gi;
                    count[self.cluster[i] as usize] += 1;
                }
                for (i, &gi) in g.iter().enumerate().filter(|(_, v)| !v.is_nan()) {
                    let c = self.cluster[i] as usize;
                    let fam_mean = sum[c] / count[c] as f64;
                    x[2 * i] = fam_mean;
                    x[2 * i + 1] = gi - fam_mean;
                }
            }
        }
        let mut xmean = vec![0.0; ng];
        for (i, xi) in x.chunks_exact(ng).enumerate() {
            if !g[i].is_nan() {
                xmean.iter_mut().zip(xi).for_each(|(m, v)| *m += v);
            }
        }
        xmean.iter_mut().for_each(|m| *m /= n as f64);
        x.iter_mut()
            .zip(xmean.iter().cycle())
            .for_each(|(v, m)| *v -= m);

        // Cross products: drop rows with a missing genotype from the base matrix, add genotype terms.
        let mut base = self.gram.clone();
        let mut level_counts = self.level_counts.clone();
        let mut cross = vec![0.0; ng * qw];
        let mut xx = vec![0.0; ng * ng];
        let mut w = vec![1.0; qw];
        for i in 0..nr {
            let zi = &self.z[i * q..(i + 1) * q];
            if g[i].is_nan() {
                w[1..=q].copy_from_slice(zi);
                w[q + 1] = self.y[i];
                for a in 0..qw {
                    let wa = w[a];
                    base[a * qw + a..(a + 1) * qw]
                        .iter_mut()
                        .zip(&w[a..])
                        .for_each(|(b, wb)| *b -= wa * wb);
                }
                for f in &self.layout.factors {
                    for l in f.start..f.start + f.levels {
                        level_counts[l] -= (zi[l] + self.zmean[l] > 0.5) as usize;
                    }
                }
                continue;
            }
            let xi = &x[i * ng..(i + 1) * ng];
            for (c, &v) in xi.iter().enumerate() {
                let row = &mut cross[c * qw..(c + 1) * qw];
                row[0] += v;
                row[1..=q]
                    .iter_mut()
                    .zip(zi)
                    .for_each(|(r, zv)| *r += v * zv);
                row[q + 1] += v * self.y[i];
                xx[c * ng..(c + 1) * ng]
                    .iter_mut()
                    .zip(xi)
                    .for_each(|(r, xv)| *r += v * xv);
            }
        }

        // Full cross-product matrix of [1, genotype terms, z, y].
        let d = ng + qw;
        let p = d - 1;
        let pos = |a: usize| if a == 0 { 0 } else { a + ng }; // [1, z, y] index -> full index
        let mut m = vec![0.0; d * d];
        let mut set = |a: usize, b: usize, v: f64| {
            m[a * d + b] = v;
            m[b * d + a] = v;
        };
        for a in 0..qw {
            for b in a..qw {
                set(pos(a), pos(b), base[a * qw + b]);
            }
        }
        for c in 0..ng {
            for a in 0..qw {
                set(1 + c, pos(a), cross[c * qw + a]);
            }
            for c2 in c..ng {
                set(1 + c, 1 + c2, xx[c * ng + c2]);
            }
        }

        // Column means (of the uncentred columns) and the factor levels absent from this fit or used
        // as reference level, which R's model matrix does not contain.
        let mean = |j: usize| match j {
            0 => 0.0,
            j if j <= ng => xmean[j - 1],
            j => self.zmean[j - 1 - ng],
        };
        let mut excluded = vec![false; p];
        for f in &self.layout.factors {
            let present: Vec<usize> = (f.start..f.start + f.levels)
                .filter(|&l| level_counts[l] > 0)
                .collect();
            if present.len() < 2 {
                return Err(FitError::SingleLevelFactor(f.name.clone()));
            }
            for l in f.start..f.start + f.levels {
                excluded[1 + ng + l] = level_counts[l] == 0 || l == present[0];
            }
        }

        // Cholesky factorisation in column order, dropping columns that are linearly dependent on the
        // preceding ones by the criterion of R's dqrdc2: residual norm < tol × original column norm.
        let mut kept: Vec<usize> = Vec::with_capacity(p);
        let mut l = vec![0.0; p * p]; // row-major, in the order of `kept`
        for j in 0..p {
            if excluded[j] {
                continue;
            }
            let r = kept.len();
            for s in 0..r {
                let dot: f64 = (0..s).map(|t| l[r * p + t] * l[s * p + t]).sum();
                l[r * p + s] = (m[j * d + kept[s]] - dot) / l[s * p + s];
            }
            let resid = m[j * d + j] - (0..r).map(|s| l[r * p + s].powi(2)).sum::<f64>();
            let mj = mean(j);
            let norm2 = m[j * d + j] + 2.0 * mj * m[j] + n as f64 * mj * mj; // uncentred
            let norm2 = if norm2 == 0.0 { 1.0 } else { norm2 };
            if resid >= QR_TOL * QR_TOL * norm2 {
                l[r * p + r] = resid.sqrt();
                kept.push(j);
            } else {
                l[r * p..r * p + r].fill(0.0);
            }
        }
        let rank = kept.len();
        let forward = |rhs: &mut [f64]| {
            for s in 0..rank {
                let dot: f64 = (0..s).map(|t| l[s * p + t] * rhs[t]).sum();
                rhs[s] = (rhs[s] - dot) / l[s * p + s];
            }
        };
        let backward = |rhs: &mut [f64]| {
            for s in (0..rank).rev() {
                let dot: f64 = (s + 1..rank).map(|t| l[t * p + s] * rhs[t]).sum();
                rhs[s] = (rhs[s] - dot) / l[s * p + s];
            }
        };
        let to_full = |v: &[f64]| {
            let mut full = vec![0.0; p];
            kept.iter().zip(v).for_each(|(&j, &x)| full[j] = x);
            full
        };

        let mut beta: Vec<f64> = kept.iter().map(|&j| m[(d - 1) * d + j]).collect();
        forward(&mut beta);
        backward(&mut beta);
        let beta = to_full(&beta); // aliased coefficients are 0 here, NA in the output

        let intercept = beta[0] + self.ymean - (1..p).map(|j| mean(j) * beta[j]).sum::<f64>();
        let is_kept = |j: usize| kept.contains(&j);
        let coef: Vec<f64> = (0..=ng)
            .map(|j| match j {
                0 => intercept,
                j if is_kept(j) => beta[j],
                _ => f64::NAN,
            })
            .collect();
        if !with_se || coef.iter().any(|c| c.is_nan()) {
            return Ok(Fit {
                n,
                coef,
                inference: None,
            });
        }

        // Rows of (X'X)⁻¹ for the reported coefficients, in the original (uncentred) parametrisation:
        // the intercept is beta_0 - Σ mean_j beta_j.
        let k = ng + 1;
        let h: Vec<Vec<f64>> = (0..k)
            .map(|t| {
                let mut a: Vec<f64> = kept
                    .iter()
                    .map(|&j| match (t, j) {
                        (0, 0) => 1.0,
                        (0, j) => -mean(j),
                        (t, j) => f64::from(u8::from(t == j)),
                    })
                    .collect();
                forward(&mut a);
                backward(&mut a);
                to_full(&a)
            })
            .collect();

        // Per-cluster sums of residual × influence.
        let mut u = vec![0.0; self.n_clusters * k];
        let mut seen = vec![false; self.n_clusters];
        for i in (0..nr).filter(|&i| !g[i].is_nan()) {
            let zi = &self.z[i * q..(i + 1) * q];
            let xi = &x[i * ng..(i + 1) * ng];
            let linear = |b: &[f64]| b[0] + dot(xi, &b[1..=ng]) + dot(zi, &b[ng + 1..]);
            let e = self.y[i] - linear(&beta);
            let c = self.cluster[i] as usize;
            seen[c] = true;
            for (t, ht) in h.iter().enumerate() {
                u[c * k + t] += e * linear(ht);
            }
        }
        let clusters = seen.iter().filter(|&&s| s).count() as f64;
        let adj = clusters / (clusters - 1.0) * (n as f64 - 1.0) / (n as f64 - rank as f64);
        let mut vcv = vec![0.0; k * k];
        for uc in u.chunks_exact(k) {
            for a in 0..k {
                for b in 0..k {
                    vcv[a * k + b] += uc[a] * uc[b];
                }
            }
        }
        vcv.iter_mut().for_each(|v| *v *= adj);

        let se: Vec<f64> = (0..k).map(|t| vcv[t * k + t].sqrt()).collect();
        let df = n as f64 - rank as f64;
        let p = coef
            .iter()
            .zip(&se)
            .map(|(b, s)| {
                let t = b / s;
                if df > 0.0 {
                    tdist::t_two_sided(t, df)
                } else {
                    tdist::normal_two_sided(t)
                }
            })
            .collect();
        Ok(Fit {
            n,
            coef,
            inference: Some(Inference { vcv, se, p }),
        })
    }
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = [0.0; 4];
    let (ca, cb) = (a.chunks_exact(4), b.chunks_exact(4));
    let tail: f64 = ca
        .remainder()
        .iter()
        .zip(cb.remainder())
        .map(|(x, y)| x * y)
        .sum();
    for (x, y) in ca.zip(cb) {
        for i in 0..4 {
            acc[i] += x[i] * y[i];
        }
    }
    acc[0] + acc[1] + acc[2] + acc[3] + tail
}

/// Upper triangle of `W'W` for `W = [1, z, y]`, computed in parallel over blocks of rows.
fn gram(z: &[f64], y: &[f64], q: usize) -> Vec<f64> {
    const BLOCK: usize = 1024;
    let qw = q + 2;
    let n = y.len();
    (0..n.div_ceil(BLOCK))
        .into_par_iter()
        .map(|b| {
            let mut acc = vec![0.0; qw * qw];
            let mut w = vec![1.0; qw];
            for i in b * BLOCK..n.min((b + 1) * BLOCK) {
                w[1..=q].copy_from_slice(&z[i * q..(i + 1) * q]);
                w[q + 1] = y[i];
                for a in 0..qw {
                    let wa = w[a];
                    acc[a * qw + a..(a + 1) * qw]
                        .iter_mut()
                        .zip(&w[a..])
                        .for_each(|(s, wb)| *s += wa * wb);
                }
            }
            acc
        })
        .reduce(
            || vec![0.0; qw * qw],
            |mut a, b| {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                a
            },
        )
}
