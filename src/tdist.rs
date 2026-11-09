//! Two-sided p-values as computed by `lmtest::coeftest`: `2 * pt(|t|, df, lower.tail = FALSE)`,
//! or `2 * pnorm(|t|, lower.tail = FALSE)` when there are no residual degrees of freedom.

const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;

/// `2 * pt(|t|, df, lower.tail = FALSE)`, i.e. the regularised incomplete beta
/// `I_x(df/2, 1/2)` with `x = df / (df + t^2)`.
pub fn t_two_sided(t: f64, df: f64) -> f64 {
    if t.is_nan() || df.is_nan() {
        return f64::NAN;
    }
    let t = t.abs();
    if t.is_infinite() {
        return 0.0;
    }
    let t2 = t * t;
    let nx = 1.0 + t2 / df;
    if nx > 1e100 {
        // Same asymptotic branch as R's pt() for extreme t.
        let lval = -0.5 * df * (2.0 * t.ln() - df.ln()) - lbeta(0.5 * df, 0.5) - (0.5 * df).ln();
        return lval.exp();
    }
    let denom = df + t2;
    inc_beta(0.5 * df, 0.5, df / denom, t2 / denom)
}

/// `2 * pnorm(|z|, lower.tail = FALSE)`.
pub fn normal_two_sided(z: f64) -> f64 {
    libm::erfc(z.abs() / std::f64::consts::SQRT_2)
}

/// Regularised incomplete beta `I_x(a, b)`, with `y = 1 - x` supplied separately so that
/// neither tail loses precision to cancellation.
fn inc_beta(a: f64, b: f64, x: f64, y: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if y <= 0.0 {
        return 1.0;
    }
    if x < (a + 1.0) / (a + b + 2.0) {
        front(a, b, x, y) * beta_cf(a, b, x) / a
    } else {
        1.0 - front(b, a, y, x) * beta_cf(b, a, y) / b
    }
}

/// `x^a y^b / B(a, b)`.
fn front(a: f64, b: f64, x: f64, y: f64) -> f64 {
    let ln_x = if x > 0.5 { (-y).ln_1p() } else { x.ln() };
    let ln_y = if y > 0.5 { (-x).ln_1p() } else { y.ln() };
    (a * ln_x + b * ln_y - lbeta(a, b)).exp()
}

/// Continued fraction for the incomplete beta function (modified Lentz).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITER: usize = 100_000;
    const EPS: f64 = 1e-15;
    const TINY: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAX_ITER {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

/// `log(B(a, b))`, following R's `lbeta()` to avoid cancellation when an argument is large.
fn lbeta(a: f64, b: f64) -> f64 {
    let p = a.min(b);
    let q = a.max(b);
    if p >= 10.0 {
        let corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
        -0.5 * q.ln()
            + LN_SQRT_2PI
            + corr
            + (p - 0.5) * (p / (p + q)).ln()
            + q * (-p / (p + q)).ln_1p()
    } else if q >= 10.0 {
        let corr = lgammacor(q) - lgammacor(p + q);
        libm::lgamma(p) + corr + p - p * (p + q).ln() + (q - 0.5) * (-p / (p + q)).ln_1p()
    } else {
        libm::lgamma(p) + libm::lgamma(q) - libm::lgamma(p + q)
    }
}

/// Stirling-series remainder `lgamma(x) - ((x - 1/2) log(x) - x + log(sqrt(2 pi)))` for `x >= 10`.
fn lgammacor(x: f64) -> f64 {
    let x2 = 1.0 / (x * x);
    (1.0 / 12.0
        - x2 * (1.0 / 360.0
            - x2 * (1.0 / 1260.0
                - x2 * (1.0 / 1680.0
                    - x2 * (1.0 / 1188.0 - x2 * (691.0 / 360_360.0 - x2 / 156.0))))))
        / x
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs()
    }

    #[test]
    fn matches_r() {
        // Reference values: sprintf("%.17g", 2 * pt(t, df, lower.tail = FALSE)) in R 4.5.
        // R uses TOMS 708 for pbeta(); the continued fraction here agrees to ~1e-11 for very large df.
        let cases = [
            (0.5, 10.0, 0.62789360574297293),
            (2.0, 30.0, 0.054625044962983108),
            (2.0, 4e5, 0.045500938784241239),
            (6.0, 4e5, 1.9748619681590899e-09),
            (40.0, 4e5, 0.0),
            (12.0, 5.0, 7.0894925171615278e-05),
            (0.01, 3.0, 0.99264911141285195),
            (1e-8, 5e4, 0.99999999202119427),
            (9.0, 123_456.0, 2.2880987613016037e-19),
            (3.3, 1e6, 0.00096688207934000285),
        ];
        for (t, df, want) in cases {
            let got = t_two_sided(t, df);
            assert!(
                close(got, want, 1e-10),
                "t={t} df={df}: got {got:e}, want {want:e}"
            );
        }
    }

    #[test]
    fn normal() {
        assert!(close(normal_two_sided(1.959963984540054), 0.05, 1e-14));
    }
}
