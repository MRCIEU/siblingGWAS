#!/usr/bin/env Rscript

VERSIONSTRING <- "0.2.0"

library(optparse, quietly = TRUE, warn.conflicts = FALSE)
library(data.table, quietly = TRUE, warn.conflicts = FALSE)

parser <- OptionParser(
  formatter = IndentedHelpFormatter,
  description = "",
  usage = "Usage: within-family-plots.R --table <input file> --out <output prefix>"
)

parser <- add_option(
  parser, c("-v", "--version"), help = "Print version", action = "callback",
  callback = function(option, flag, value, parser) {
    cat(paste("siblingGWAS version: ", VERSIONSTRING, "\n"), sep = "")
    quit(save = "no", status = 0)
  }
)

parser <- add_option(parser, "--table", help = "The output of within-family.R", type = "character", required = T)
parser <- add_option(parser, "--out", help = "A string that will be prefixed to the output file according to: paste0(<output prefix>, \"_pop_vs_family_beta.pdf\")", type = "character", required = T)

args <- parse_args(parser)

ifile <- args$table
phenname <- args$out

expected_full <- c("CHR", "SNP", "BP", "A1", "A2", "N_REG", "BETA_MODEL1_0", "BETA_MODEL2_0", "BETA_TOTAL", "BETA_BF", "BETA_WF", "SE_BETA_MODEL1_0", "SE_BETA_MODEL2_0", "SE_BETA_TOTAL", "SE_BETA_BF", "SE_BETA_WF", "P_BETA_MODEL1_0", "P_BETA_MODEL2_0", "P_BETA_TOTAL", "P_BETA_BF", "P_BETA_WF", "VCV_MODEL1_0", "VCV_MODEL1_0_TOTAL", "VCV_MODEL1_TOTAL", "VCV_MODEL2_0", "VCV_MODEL2_0_BF", "VCV_MODEL2_0_WF", "VCV_MODEL2_BF", "VCV_MODEL2_BF_WF", "VCV_MODEL2_WF")
expected_sub <- c("CHR", "SNP", "BP", "A1", "A2", "N_REG", "BETA_MODEL1_0", "BETA_MODEL2_0", "BETA_TOTAL", "BETA_BF", "BETA_WF")

df <- fread(ifile)

if(all(expected_full %in% colnames(df))){
  no_se <- FALSE
} else if(all(expected_sub %in% colnames(df))){
  no_se <- TRUE
} else {
  cat("Expected column names not found\n")
  quit(save = "no", status = 1)
}

for(i in 1:nrow(df)) {
  if(df[i, "BETA_TOTAL"] < 0) {
    df[i, "BETA_TOTAL"] <- df[i, "BETA_TOTAL"] * -1
    df[i, "BETA_BF"] <- df[i, "BETA_BF"] * -1
    df[i, "BETA_WF"] <- df[i, "BETA_WF"] * -1
  }
}

###

if(no_se) {
  df_abs <- subset(df, select = c(
      BETA_TOTAL,
      BETA_BF,
      BETA_WF
    )
  )

  axis_limits <- range(
    df_abs$BETA_TOTAL,
    df_abs$BETA_BF,
    df_abs$BETA_WF
  )

} else {
  df_abs <- subset(df, select = c(
    BETA_TOTAL,
    SE_BETA_TOTAL,
    BETA_BF,
    SE_BETA_BF,
    BETA_WF,
    SE_BETA_WF
  ))

  axis_limits <- range(
    df_abs$BETA_TOTAL - df_abs$SE_BETA_TOTAL,
    df_abs$BETA_TOTAL + df_abs$SE_BETA_TOTAL,
    df_abs$BETA_BF - df_abs$SE_BETA_BF,
    df_abs$BETA_BF + df_abs$SE_BETA_BF,
    df_abs$BETA_WF - df_abs$SE_BETA_WF,
    df_abs$BETA_WF + df_abs$SE_BETA_WF
  )
}

###


pdf(
  file = paste0(phenname, "_pop_vs_family_beta.pdf"),
  height = 7,
  width = 14
)

par(mfrow=c(1,2))

regression_pop_bf <- lm(df_abs$BETA_BF ~ df_abs$BETA_TOTAL)
plot(
  df_abs$BETA_TOTAL,
  df_abs$BETA_BF,
  type = "n",
  xlim = axis_limits,
  ylim = axis_limits,
  xlab = "absolute population beta ± 1SE",
  ylab = "absolute between-family beta ± 1SE",
  main = ""
)
abline(
  a = 0,
  b = 1,
  col = "darkgrey",
  lty = "dashed"
)
if(!no_se) {
  arrows(
    x0 = df_abs$BETA_TOTAL - df_abs$SE_BETA_TOTAL,
    x1 = df_abs$BETA_TOTAL + df_abs$SE_BETA_TOTAL,
    y0 = df_abs$BETA_BF,
    length = 0.05,
    angle = 90,
    code = 3,
    col = "gray"
  )
  arrows(
    y0 = df_abs$BETA_BF - df_abs$SE_BETA_BF,
    y1 = df_abs$BETA_BF + df_abs$SE_BETA_BF,
    x0 = df_abs$BETA_TOTAL,
    length = 0.05,
    angle = 90,
    code = 3,
    col = "gray"
  )
}
points(
  df_abs$BETA_TOTAL,
  df_abs$BETA_BF,
  pch = 21,
  bg = "white"
)
abline(
  regression_pop_bf,
  col = "red"
)
legend(
  "bottomright",
  legend = c(
    paste("slope: ", round(regression_pop_bf$coefficients[[2]], 2))
  ),
  bty = "n",
  pch = c("_"),
  col = "red"
)
mtext(
  "between-family",
  side = 3,
  line = 0.5,
  cex = 1.1
)

mtext(
  phenname,
  side = 3,
  line = -2,
  outer = T,
  cex = 1.4,
  font = 2
)

###

regression_pop_wf <- lm(df_abs$BETA_WF ~ df_abs$BETA_TOTAL)
plot(
  df_abs$BETA_TOTAL,
  df_abs$BETA_WF,
  type = "n",
  xlim = axis_limits,
  ylim = axis_limits,
  xlab = "absolute population beta ± 1SE",
  ylab = "absolute within-family beta ± 1SE",
  main = ""
)
abline(
  a = 0,
  b = 1,
  col = "darkgrey",
  lty = "dashed"
)
if(!no_se) {
  arrows(
    x0 = df_abs$BETA_TOTAL - df_abs$SE_BETA_TOTAL,
    x1 = df_abs$BETA_TOTAL + df_abs$SE_BETA_TOTAL,
    y0 = df_abs$BETA_WF,
    length = 0.05,
    angle = 90,
    code = 3,
    col = "gray"
  )
  arrows(
    y0 = df_abs$BETA_WF - df_abs$SE_BETA_WF,
    y1 = df_abs$BETA_WF + df_abs$SE_BETA_WF,
    x0 = df_abs$BETA_TOTAL,
    length = 0.05,
    angle = 90,
    code = 3,
    col = "gray"
  )
}
points(
  df_abs$BETA_TOTAL,
  df_abs$BETA_WF,
  pch = 21,
  bg = "white"
)
abline(
  regression_pop_wf,
  col = "red"
)
legend(
  "bottomright",
  legend = c(
    paste("slope: ", round(regression_pop_wf$coefficients[[2]], 2))
  ),
  bty = "n",
  pch = c("_"),
  col = "red"
)
mtext(
  "within-family",
  side = 3,
  line = 0.5,
  cex = 1.1
)

dev.off()
