library(GEOquery)
library(dplyr)

outdir <- "data"
dir.create(outdir, showWarnings = FALSE)

# ── GSE50660 ──────────────────────────────────────────────────────────────────
# Liu et al. 2013 - 464 whole blood samples, 450k, smoking status available
cat("Downloading GSE50660...\n")
gse50660 <- getGEO("GSE50660", GSEMatrix = TRUE, destdir = outdir)
gse50660 <- gse50660[[1]]

# Extract phenotype data
pheno50660 <- pData(gse50660)
cat("GSE50660 samples:", nrow(pheno50660), "\n")
cat("Columns:", paste(colnames(pheno50660)[1:20], collapse=", "), "\n")

# Look for smoking status column
smoking_cols <- grep("smok|tobacco|cigarette", colnames(pheno50660), 
                     ignore.case=TRUE, value=TRUE)
cat("Smoking-related columns:", paste(smoking_cols, collapse=", "), "\n")

if (length(smoking_cols) > 0) {
    cat("Smoking status values:\n")
    print(table(pheno50660[[smoking_cols[1]]]))
}

saveRDS(gse50660, file.path(outdir, "GSE50660.rds"))
cat("Saved GSE50660.rds\n")

# ── GSE42861 ──────────────────────────────────────────────────────────────────
# Liu et al. 2013 rheumatoid arthritis cohort - 689 whole blood samples, 450k,
# smoking status available (never / ex / current / occasional). This is the
# second cohort used by the rest of the smoking pipeline (see 02_preprocess.R).
cat("Downloading GSE42861...\n")
gse42861 <- getGEO("GSE42861", GSEMatrix = TRUE, destdir = outdir)
gse42861 <- gse42861[[1]]

pheno42861 <- pData(gse42861)
cat("GSE42861 samples:", nrow(pheno42861), "\n")

smoking_cols <- grep("smok|tobacco|cigarette", colnames(pheno42861),
                     ignore.case=TRUE, value=TRUE)
cat("Smoking-related columns:", paste(smoking_cols, collapse=", "), "\n")

if (length(smoking_cols) > 0) {
    cat("Smoking status values:\n")
    print(table(pheno42861[[smoking_cols[1]]]))
}

saveRDS(gse42861, file.path(outdir, "GSE42861.rds"))
cat("Saved GSE42861.rds\n")

# ── Check methylation data availability ───────────────────────────────────────
cat("\nGSE50660 feature dimensions:", nrow(gse50660), "CpGs x", ncol(gse50660), "samples\n")
cat("GSE42861 feature dimensions:", nrow(gse42861), "CpGs x", ncol(gse42861), "samples\n")
