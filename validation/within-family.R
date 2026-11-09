#!/usr/bin/env Rscript

VERSIONSTRING <- "0.2.0"

# The code below is modified from: https://github.com/LaurenceHowe/SiblingGWAS/
# which is distributed under licence:
#
# MIT License
#
# Copyright (c) 2020 LaurenceHowe
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library(data.table, quietly = TRUE, warn.conflicts = FALSE)
library(sandwich, quietly = TRUE, warn.conflicts = FALSE)
library(lmtest, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

parser <- OptionParser(
  formatter = IndentedHelpFormatter,
  description = "",
  usage = "Usage: within-family.R --raw <raw file> --bim <bim file> --phenotypes <phenotype file> --covariates <covariate file> --siblings <sibling file> --out <output prefix> [--no-se]"
)

parser <- add_option(
  parser, c("-v", "--version"), help = "Print version", action = "callback",
  callback = function(option, flag, value, parser) {
    cat(paste("siblingGWAS version: ", VERSIONSTRING, "\n"), sep = "")
    quit(save = "no", status = 0)
  }
)

parser <- add_option(parser, "--raw", help = "The output of plink2 --recode A (\"Sample-major additive (0/1/2) coding, suitable for loading from R...\")", type = "character", required = T)
parser <- add_option(parser, "--bim", help = " A plink extended variant information file", type = "character", required = T)
parser <- add_option(parser, "--phenotypes", help = "A plain text (separated value) file consisting of EXACTLY three columns: FID IID \"phenotype\" whe\" where \"phenotype\" can be named sensibly", type = "character", required = T)
parser <- add_option(parser, "--covariates", help = "A plain text (separated value) file consisting of AT LEAST three columns. First: FID IID + at lease one column containing covariates: <\"covariate 1\"> <\"covariate 2\"> ..., etc", type = "character", required = T)
parser <- add_option(parser, "--siblings", help = "A plain text (separated value) file containing sibling information with exactly two columns: IID FID_FS where FID_FS is a family identifier for dizygotic siblings", type = "character", required = T)
parser <- add_option(parser, "--out", help = "A string that will be prefixed to the output file according to: paste0(<output prefix>, \"_within-family-table.tsv\")", type = "character", required = T)

parser <- add_option(parser, "--no-se", help = "Don't calculate p-values and standard errors", type = "logical", action = "store_true", default = FALSE, dest = "no_se")

args <- parse_args(parser)


rawfile <- args$raw
bimfile <- args$bim
phenfile <- args$phenotypes
covfile <- args$covariates
sibfile <- args$siblings
outprefix <- args$out


raw <- fread(rawfile) # FID	IID	PAT	MAT	SEX	PHENOTYPE	rs35665085_G ... <other snps (appended with sample major ref allele - i.e. the allele that has been counted)>
bim <- fread(bimfile) # No header: CHR RSID Map_position BP_coordinate Allele1 Allele2
phenos <- fread(phenfile) # FID	IID	participant.p50_i0
cov <- fread(covfile) # FID	IID	... <covariates?>

covariate_names <- colnames(cov)[3:ncol(cov)] # CHANGE THIS?

sibs <- fread(sibfile) # IID	FID_FS

phenos <- phenos[,1:3] # CHANGE THIS

original_phenotype_column_name <- colnames(phenos)[3]
colnames(phenos) <- c("FID", "IID", "Phenotype")

phencov <- merge(phenos, cov, by=c("FID", "IID"))
rawphencov <- merge(raw, phencov, by=c("FID", "IID"))
df_input <- merge(rawphencov, sibs, by = c("IID"), all.x = TRUE, all.y = FALSE)

df_input$IS_FS <- apply(df_input, 1, FUN = function(rw) {
  if(is.na(rw["FID_FS"])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})

df_input$FID_FS_FILLED <- apply(df_input, 1, FUN = function(rw) {
  if(is.na(rw["FID_FS"])) {
    return(rw["FID"])
  } else {
    return(rw["FID_FS"])
  }
})

output <- data.frame(CHR=bim$V1, SNP=bim$V2, BP=bim$V4, A1=bim$V5, A2=bim$V6, N_REG=NA,
                     BETA_MODEL1_0=NA, BETA_MODEL2_0=NA, BETA_TOTAL=NA, BETA_BF=NA, BETA_WF=NA,
                     SE_BETA_MODEL1_0=NA, SE_BETA_MODEL2_0=NA, SE_BETA_TOTAL=NA, SE_BETA_BF=NA, SE_BETA_WF=NA,
                     P_BETA_MODEL1_0=NA, P_BETA_MODEL2_0=NA, P_BETA_TOTAL=NA, P_BETA_BF=NA, P_BETA_WF=NA,
                     VCV_MODEL1_0=NA, VCV_MODEL1_0_TOTAL=NA, VCV_MODEL1_TOTAL=NA,
                     VCV_MODEL2_0=NA, VCV_MODEL2_0_BF=NA, VCV_MODEL2_0_WF=NA, VCV_MODEL2_BF=NA, VCV_MODEL2_BF_WF=NA, VCV_MODEL2_WF=NA)


for(i in 1:nrow(bim)) {
  snp_ind <- i+6

  df_all <- data.table(FID=df_input$FID_FS_FILLED,
                       IS_FS=df_input$IS_FS,
                       PHENOTYPE=df_input$Phenotype,
                       GENOTYPE=as.numeric(unlist(df_input[,snp_ind, with=F])))
  for(covariate in covariate_names) {
    df_all[[toupper(covariate)]] <- df_input[[covariate]]
  }

  df_all_na_omitted <- na.omit(df_all)

  formula_model1 <- as.formula(paste("PHENOTYPE ~ GENOTYPE +", paste(toupper(covariate_names), collapse = "+")))

  skip_variant <- FALSE
  tryCatch(
    model1 <- lm(formula = formula_model1, data=df_all_na_omitted),
    error = function(e) {
      print(e)
      skip_variant <<- TRUE
    }
  )
  if(skip_variant) { next }

  #Save Beta information
  output$BETA_MODEL1_0[i] <- model1$coefficients[1]
  output$BETA_TOTAL[i] <- model1$coefficients[2]

  # Save the variance covariance matrix to cluster SEs by family
  # Try and catch errors with generating variance covariance matrix

  if(!args$no_se) {
    tryCatch(
      vcv_matrix <- vcovCL(model1, cluster=df_all_na_omitted$FID),
      error = function(e) {
        print(e)
        skip_variant <<-TRUE
      }
    )
    if(skip_variant) { next }

    if(  is.na(output$BETA_MODEL1_0[i]) | is.na(output$BETA_TOTAL[i])) {
      output$VCV_MODEL1_0[i] <-NA
      output$VCV_MODEL1_0_TOTAL[i] <-NA
      output$VCV_MODEL1_TOTAL[i] <-NA
    } else {
      output$VCV_MODEL1_0[i] <- vcv_matrix[1,1]
      output$VCV_MODEL1_0_TOTAL[i] <- vcv_matrix[1,2]
      output$VCV_MODEL1_TOTAL[i] <- vcv_matrix[2,2]

    }

    #Derive the clustered SEs for the total effect and P-values
    #Try and catch errors with clustered standard errors

    tryCatch(
      test_matrix <- coeftest(model1, vcov.=vcv_matrix),
      error = function(e) {
        print(e)
        skip_variant <<-TRUE
      }
    )
    if(skip_variant) { next }

    if(  is.na(output$BETA_MODEL1_0[i]) | is.na(output$BETA_TOTAL[i])) {
      output$SE_BETA_MODEL1_0[i] <- NA
      output$SE_BETA_TOTAL[i] <- NA
      output$P_BETA_MODEL1_0[i] <- NA
      output$P_BETA_TOTAL[i] <- NA
    } else {
      output$SE_BETA_MODEL1_0[i] <- test_matrix[1,2]
      output$SE_BETA_TOTAL[i] <- test_matrix[2,2]
      output$P_BETA_MODEL1_0[i] <- test_matrix[1,4]
      output$P_BETA_TOTAL[i] <- test_matrix[2,4]
    }
  }

  df_fam_na_omitted <- df_all_na_omitted[df_all_na_omitted$IS_FS]
  df_fam_na_omitted$FAM_MEAN <- ave(as.numeric(unlist(df_fam_na_omitted$GENOTYPE)), df_fam_na_omitted$FID, FUN=mean)
  df_fam_na_omitted$CENTRED_GENOTYPE <- df_fam_na_omitted$GENOTYPE - df_fam_na_omitted$FAM_MEAN

  formula_model2 <- as.formula(paste("PHENOTYPE ~ FAM_MEAN + CENTRED_GENOTYPE +", paste(toupper(covariate_names), collapse = "+")))

  # Run unified regression
  skip_variant <- FALSE
  tryCatch(
    model2 <- lm(formula = formula_model2, data=df_fam_na_omitted),
    error = function(e) {
      print(e)
      skip_variant <<- TRUE
    }
  )
  if(skip_variant) { next }

  # Sample size in regression
  output$N_REG[i] <- length(resid(model2))

  # Save Beta information
  output$BETA_MODEL2_0[i] <- model2$coefficients[1]
  output$BETA_BF[i] <- model2$coefficients[2]
  output$BETA_WF[i] <- model2$coefficients[3]


  if(!args$no_se) {
    # save the variance covariance matrix
    tryCatch(
      vcv_matrix <- vcovCL(model2, cluster=df_fam_na_omitted$FID),
      error = function(e){
        print(e)
        skip_variant <<- TRUE
      }
    )
    if(skip_variant) { next }

    if(  is.na(output$BETA_MODEL2_0[i]) | is.na(output$BETA_BF[i]) | is.na(output$BETA_WF[i]) ) {
      output$VCV_MODEL2_0[i] <-NA
      output$VCV_MODEL2_0_BF[i] <-NA
      output$VCV_MODEL2_0_WF[i] <-NA
      output$VCV_MODEL2_BF[i] <-NA
      output$VCV_MODEL2_BF_WF[i] <-NA
      output$VCV_MODEL2_WF[i] <-NA
    } else {
      output$VCV_MODEL2_0[i] <- vcv_matrix[1,1]
      output$VCV_MODEL2_0_BF[i] <- vcv_matrix[1,2]
      output$VCV_MODEL2_0_WF[i] <- vcv_matrix[1,3]
      output$VCV_MODEL2_BF[i] <- vcv_matrix[2,2]
      output$VCV_MODEL2_BF_WF[i] <- vcv_matrix[2,3]
      output$VCV_MODEL2_WF[i] <- vcv_matrix[3,3]
    }

    # save the clustered SEs and corresponding P-values for WF/BF
    tryCatch(
      test_matrix <- coeftest(model2, vcov.=vcv_matrix),
      error = function(e) {
        print(e)
        skip_variant <<-TRUE
      }
    )
    if(skip_variant) { next }

    if(  is.na(output$BETA_MODEL2_0[i]) | is.na(output$BETA_BF[i]) | is.na(output$BETA_WF[i]) ) {
      output$SE_BETA_MODEL2_0[i] <- NA
      output$SE_BETA_BF[i] <- NA
      output$SE_BETA_WF[i] <- NA
      output$P_BETA_MODEL2_0[i] <- NA
      output$P_BETA_BF[i] <- NA
      output$P_BETA_WF[i] <- NA
    } else {
      output$SE_BETA_MODEL2_0[i] <- test_matrix[1,2]
      output$SE_BETA_BF[i] <- test_matrix[2,2]
      output$SE_BETA_WF[i] <- test_matrix[3,2]
      output$P_BETA_MODEL2_0[i] <- test_matrix[1,4]
      output$P_BETA_BF[i] <- test_matrix[2,4]
      output$P_BETA_WF[i] <- test_matrix[3,4]
    }
  }
}

if(args$no_se) {
  output <- subset(output, select = c(CHR, SNP, BP, A1, A2, N_REG, BETA_MODEL1_0, BETA_MODEL2_0, BETA_TOTAL, BETA_BF, BETA_WF))
}

fwrite(output, file = paste0(outprefix, "_within-family-table.tsv"), sep="\t")
