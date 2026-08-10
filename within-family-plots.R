
USAGESTRING <- "Usage: within-family-plots.R <input file> <phenotype name>"

arguments <- commandArgs(trailingOnly = T)

for(argument in arguments) {
  if(argument == "-h" || argument == "--help") {
    cat(USAGESTRING)
    quit(save = "no", status = 0)
  }
}

if(length(arguments) != 2) {
  cat("Error: Wrong number of arguments provided.\n")
  cat(USAGESTRING)
  quit(save = "no", status = 1)
}

library(data.table, quietly = TRUE, warn.conflicts = FALSE)

ifile <- arguments[1]
phenname <- arguments[2]

df <- fread(ifile)

###

df_abs <- subset(df, select = c(
  BETA_TOTAL,
  SE_BETA_TOTAL,
  BETA_BF,
  SE_BETA_BF,
  BETA_WF,
  SE_BETA_WF
))

for(i in 1:nrow(df_abs)) {
  if(df_abs[i, "BETA_TOTAL"] < 0) {
    df_abs[i, "BETA_TOTAL"] <- df_abs[i, "BETA_TOTAL"] * -1
    df_abs[i, "BETA_BF"] <- df_abs[i, "BETA_BF"] * -1
    df_abs[i, "BETA_WF"] <- df_abs[i, "BETA_WF"] * -1
  }
}

###

axis_limits <- range(
  df_abs$BETA_TOTAL - df_abs$SE_BETA_TOTAL, 
  df_abs$BETA_TOTAL + df_abs$SE_BETA_TOTAL, 
  df_abs$BETA_BF - df_abs$SE_BETA_BF,
  df_abs$BETA_BF + df_abs$SE_BETA_BF,
  df_abs$BETA_WF - df_abs$SE_BETA_WF,
  df_abs$BETA_WF + df_abs$SE_BETA_WF
)

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
