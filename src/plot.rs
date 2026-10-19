//! The `<prefix>_pop_vs_family_beta.pdf` plot (`--plot`): a port of siblingGWAS's
//! `within-family-plots.R`, drawn with rs-graphics. Plots the absolute population effect of each
//! variant against its between-family (left) and within-family (right) effects, with ±1 SE bars
//! when standard errors were calculated.

use crate::output::VariantResult;
use anyhow::{bail, Context, Result};
use rs_graphics::backend::pdf::PdfBackend;
use rs_graphics::prelude::*;
use std::path::Path;

/// Effects, oriented so that the population effect is positive.
struct Effects {
    total: Vec<f64>,
    bf: Vec<f64>,
    wf: Vec<f64>,
    /// SEs of total, BF and WF, unless standard errors weren't calculated.
    se: Option<[Vec<f64>; 3]>,
}

/// Writes the plot to `path`, with `title` above both panels.
pub fn write(path: &Path, title: &str, results: &[VariantResult], with_se: bool) -> Result<()> {
    // VariantResult's estimate order: MODEL1_0, MODEL2_0, TOTAL, BF, WF.
    let column = |values: fn(&VariantResult) -> &[f64; 5], k: usize| -> Vec<f64> {
        results.iter().map(|r| values(r)[k]).collect()
    };
    let mut fx = Effects {
        total: column(|r| &r.beta, 2),
        bf: column(|r| &r.beta, 3),
        wf: column(|r| &r.beta, 4),
        se: with_se.then(|| {
            [
                column(|r| &r.se, 2),
                column(|r| &r.se, 3),
                column(|r| &r.se, 4),
            ]
        }),
    };
    // The R script errors on a missing BETA_TOTAL here; such variants are just left out of the plot.
    for i in 0..fx.total.len() {
        if fx.total[i] < 0.0 {
            fx.total[i] = -fx.total[i];
            fx.bf[i] = -fx.bf[i];
            fx.wf[i] = -fx.wf[i];
        }
    }

    let limits = axis_limits(&fx)?;
    let mut dev = Device::new(PdfBackend::new(14.0, 7.0), Par::default());
    dev.mfrow(1, 2);
    panel(&mut dev, &fx, limits, Family::Between);
    // mtext(phenname, side = 3, line = -2, outer = TRUE, cex = 1.4, font = 2)
    dev.mtext(
        title,
        3,
        &MtextOptions::default()
            .line(-2.0)
            .outer(true)
            .cex(1.4)
            .font(2),
    );
    panel(&mut dev, &fx, limits, Family::Within);
    std::fs::write(path, dev.into_backend().finish())
        .with_context(|| format!("cannot write {}", path.display()))
}

/// The common range of all effects (± their SEs), used for both axes of both panels.
/// R's `range()` here returns NA (and the script fails) if any value is missing; missing values
/// are ignored instead.
fn axis_limits(fx: &Effects) -> Result<(f64, f64)> {
    let mut values: Vec<f64> = Vec::new();
    match &fx.se {
        None => {
            for v in [&fx.total, &fx.bf, &fx.wf] {
                values.extend(v.iter());
            }
        }
        Some(se) => {
            for (v, s) in [(&fx.total, &se[0]), (&fx.bf, &se[1]), (&fx.wf, &se[2])] {
                values.extend(v.iter().zip(s).flat_map(|(v, s)| [v - s, v + s]));
            }
        }
    }
    let finite = values.iter().copied().filter(|v| v.is_finite());
    let lo = finite.clone().fold(f64::INFINITY, f64::min);
    let hi = finite.fold(f64::NEG_INFINITY, f64::max);
    if lo > hi {
        bail!("no variant has finite effect estimates to plot");
    }
    Ok((lo, hi))
}

#[derive(Clone, Copy)]
enum Family {
    Between,
    Within,
}

fn panel(dev: &mut Device<PdfBackend>, fx: &Effects, limits: (f64, f64), family: Family) {
    let (y, y_se, name) = match family {
        Family::Between => (&fx.bf, fx.se.as_ref().map(|s| &s[1]), "between-family"),
        Family::Within => (&fx.wf, fx.se.as_ref().map(|s| &s[2]), "within-family"),
    };
    let x = &fx.total;
    let ylab = format!("absolute {name} beta ± 1SE");
    dev.plot(
        x,
        y,
        &PlotOptions::default()
            .type_(PlotType::None)
            .xlim(limits.0, limits.1)
            .ylim(limits.0, limits.1)
            .xlab("absolute population beta ± 1SE")
            .ylab(ylab),
    );
    // abline(0, 1, col = "darkgrey", lty = "dashed")
    dev.abline(
        &AblineOptions::default()
            .coef(0.0, 1.0)
            .col("darkgrey")
            .lty(LineType::DASHED),
    );

    if let (Some(x_se), Some(y_se)) = (fx.se.as_ref().map(|s| &s[0]), y_se) {
        let gray = ArrowsOptions::default()
            .length(0.05)
            .angle(90.0)
            .code(3)
            .col("gray");
        let lo = |v: &[f64], s: &[f64]| v.iter().zip(s).map(|(v, s)| v - s).collect::<Vec<_>>();
        let hi = |v: &[f64], s: &[f64]| v.iter().zip(s).map(|(v, s)| v + s).collect::<Vec<_>>();
        // arrows(x0, x1, y0, length = 0.05, angle = 90, code = 3, col = "gray"), horizontal then vertical
        dev.arrows(&lo(x, x_se), y, &hi(x, x_se), y, &gray);
        dev.arrows(x, &lo(y, y_se), x, &hi(y, y_se), &gray);
    }
    // points(x, y, pch = 21, bg = "white")
    dev.points(x, y, &XyOptions::default().pch(21).col(1).bg("white"));

    let red = ColorSpec::from("red");
    let (intercept, slope) = simple_lm(x, y);
    if intercept.is_finite() && slope.is_finite() {
        dev.abline(&AblineOptions::default().coef(intercept, slope).col(red));
    }
    // paste("slope: ", round(slope, 2))
    let label = format!("slope:  {}", format_round2(slope));
    dev.legend(
        LegendPos::Keyword(LegendKeyword::BottomRight),
        &[LegendEntry::point(label, red, '_' as i32)],
        &LegendOptions::default().bty(LegendBty::None),
    );
    // mtext(name, side = 3, line = 0.5, cex = 1.1)
    dev.mtext(name, 3, &MtextOptions::default().line(0.5).cex(1.1));
}

/// Intercept and slope of `lm(y ~ x)` on the complete pairs. The slope is NaN when it is not
/// estimable (as `lm` reports NA for an aliased coefficient).
fn simple_lm(x: &[f64], y: &[f64]) -> (f64, f64) {
    let pairs: Vec<(f64, f64)> = x
        .iter()
        .zip(y)
        .map(|(&x, &y)| (x, y))
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .collect();
    let n = pairs.len() as f64;
    if pairs.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let mx = pairs.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pairs.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pairs.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let sxy: f64 = pairs.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let sx_raw: f64 = pairs.iter().map(|p| p.0 * p.0).sum();
    // lm's (dqrdc2) test: the centred column norm must be at least 1e-7 of the original norm.
    if sxx < 1e-14 * sx_raw || sxx == 0.0 {
        return (my, f64::NAN);
    }
    let slope = sxy / sxx;
    (my - slope * mx, slope)
}

/// `as.character(round(v, 2))`.
fn format_round2(v: f64) -> String {
    if !v.is_finite() {
        return "NA".to_string();
    }
    let r: f64 = format!("{v:.2}").parse().unwrap();
    if r == 0.0 {
        "0".to_string()
    } else {
        format!("{r}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round2() {
        assert_eq!(format_round2(0.5349), "0.53");
        assert_eq!(format_round2(1.0), "1");
        assert_eq!(format_round2(1.199), "1.2");
        assert_eq!(format_round2(-0.001), "0");
        assert_eq!(format_round2(f64::NAN), "NA");
    }

    #[test]
    fn lm_slope() {
        let (a, b) = simple_lm(&[1.0, 2.0, 3.0, f64::NAN], &[3.0, 5.0, 7.0, 1.0]);
        assert!((a - 1.0).abs() < 1e-12 && (b - 2.0).abs() < 1e-12);
        assert!(simple_lm(&[2.0, 2.0], &[1.0, 3.0]).1.is_nan());
    }
}
