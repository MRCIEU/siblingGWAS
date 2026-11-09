//! The `<prefix>_within-family-table.tsv` output, formatted like `data.table::fwrite`.

use crate::input::BimRow;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Results for one variant. NaN means NA.
pub struct VariantResult {
    pub n_reg: Option<usize>,
    /// MODEL1_0, MODEL2_0, TOTAL, BF, WF
    pub beta: [f64; 5],
    pub se: [f64; 5],
    pub p: [f64; 5],
    /// MODEL1_0, MODEL1_0_TOTAL, MODEL1_TOTAL
    pub vcv1: [f64; 3],
    /// MODEL2_0, MODEL2_0_BF, MODEL2_0_WF, MODEL2_BF, MODEL2_BF_WF, MODEL2_WF
    pub vcv2: [f64; 6],
}

impl Default for VariantResult {
    fn default() -> Self {
        VariantResult {
            n_reg: None,
            beta: [f64::NAN; 5],
            se: [f64::NAN; 5],
            p: [f64::NAN; 5],
            vcv1: [f64::NAN; 3],
            vcv2: [f64::NAN; 6],
        }
    }
}

const ID_COLUMNS: [&str; 6] = ["CHR", "SNP", "BP", "A1", "A2", "N_REG"];
const ESTIMATES: [&str; 5] = ["MODEL1_0", "MODEL2_0", "TOTAL", "BF", "WF"];
const VCV_COLUMNS: [&str; 9] = [
    "VCV_MODEL1_0",
    "VCV_MODEL1_0_TOTAL",
    "VCV_MODEL1_TOTAL",
    "VCV_MODEL2_0",
    "VCV_MODEL2_0_BF",
    "VCV_MODEL2_0_WF",
    "VCV_MODEL2_BF",
    "VCV_MODEL2_BF_WF",
    "VCV_MODEL2_WF",
];

pub fn write(path: &Path, bim: &[BimRow], results: &[VariantResult], with_se: bool) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    let mut out = BufWriter::new(file);
    let mut header: Vec<String> = ID_COLUMNS.iter().map(|s| s.to_string()).collect();
    header.extend(ESTIMATES.iter().map(|e| format!("BETA_{e}")));
    if with_se {
        header.extend(ESTIMATES.iter().map(|e| format!("SE_BETA_{e}")));
        header.extend(ESTIMATES.iter().map(|e| format!("P_BETA_{e}")));
        header.extend(VCV_COLUMNS.iter().map(|s| s.to_string()));
    }
    writeln!(out, "{}", header.join("\t"))?;

    for (b, r) in bim.iter().zip(results) {
        let n_reg = r.n_reg.map(|n| n.to_string()).unwrap_or_default();
        let mut fields = vec![
            b.chr.clone(),
            b.snp.clone(),
            b.bp.clone(),
            b.a1.clone(),
            b.a2.clone(),
            n_reg,
        ];
        let mut numbers: Vec<f64> = r.beta.to_vec();
        if with_se {
            numbers.extend(r.se);
            numbers.extend(r.p);
            numbers.extend(r.vcv1);
            numbers.extend(r.vcv2);
        }
        fields.extend(numbers.into_iter().map(format_number));
        writeln!(out, "{}", fields.join("\t"))?;
    }
    out.flush()?;
    Ok(())
}

/// Formats like `fwrite`: 15 significant digits, scientific notation only when it is shorter,
/// NA/NaN as an empty field.
pub fn format_number(v: f64) -> String {
    if v.is_nan() {
        return String::new();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let s = format!("{:.14e}", v.abs());
    let (mantissa, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits: String = mantissa.chars().filter(|c| c.is_ascii_digit()).collect();
    let digits = digits.trim_end_matches('0');
    let sign = if v < 0.0 { "-" } else { "" };

    let mut sci = format!("{sign}{}", &digits[..1]);
    if digits.len() > 1 {
        sci.push('.');
        sci.push_str(&digits[1..]);
    }
    sci.push_str(&format!(
        "e{}{:02}",
        if exp < 0 { '-' } else { '+' },
        exp.abs()
    ));

    let fixed = if exp >= 0 {
        let int_len = exp as usize + 1;
        if digits.len() <= int_len {
            format!("{sign}{digits}{}", "0".repeat(int_len - digits.len()))
        } else {
            format!("{sign}{}.{}", &digits[..int_len], &digits[int_len..])
        }
    } else {
        format!("{sign}0.{}{digits}", "0".repeat((-exp - 1) as usize))
    };
    if fixed.len() <= sci.len() {
        fixed
    } else {
        sci
    }
}

#[cfg(test)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::format_number;

    #[test]
    fn matches_fwrite() {
        // Expected strings from data.table::fwrite (1.18).
        let cases = [
            (1.234567890123456789, "1.23456789012346"),
            (1e-300, "1e-300"),
            (123456789.123, "123456789.123"),
            (0.0001234, "0.0001234"),
            (1.0 / 3.0, "0.333333333333333"),
            (2e15, "2e+15"),
            (1e5, "1e+05"),
            (100000.5, "100000.5"),
            (-2.5e-7, "-2.5e-07"),
            (-0.5, "-0.5"),
            (0.0, "0"),
            (f64::INFINITY, "Inf"),
            (f64::NAN, ""),
        ];
        for (v, want) in cases {
            assert_eq!(format_number(v), want, "{v}");
        }
    }
}
