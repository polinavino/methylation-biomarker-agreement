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

# ── GSE110043 ─────────────────────────────────────────────────────────────────
# 94 samples, whole blood, 450k, smoking status available
cat("Downloading GSE110043...\n")
gse110043 <- getGEO("GSE110043", GSEMatrix = TRUE, destdir = outdir)
gse110043 <- gse110043[[1]]

pheno110043 <- pData(gse110043)
cat("GSE110043 samples:", nrow(pheno110043), "\n")

smoking_cols <- grep("smok|tobacco|cigarette", colnames(pheno110043),
                     ignore.case=TRUE, value=TRUE)
cat("Smoking-related columns:", paste(smoking_cols, collapse=", "), "\n")

if (length(smoking_cols) > 0) {
    cat("Smoking status values:\n")
    print(table(pheno110043[[smoking_cols[1]]]))
}

saveRDS(gse110043, file.path(outdir, "GSE110043.rds"))
cat("Saved GSE110043.rds\n")

# ── Check methylation data availability ───────────────────────────────────────
cat("\nGSE50660 feature dimensions:", nrow(gse50660), "CpGs x", ncol(gse50660), "samples\n")
cat("GSE50660 assay names:", assayNames(gse50660), "\n")

# Check what phenotype columns GSE110043 actually has
cat("\nGSE110043 phenotype columns:\n")
print(colnames(pData(gse110043)))
cat("\nFirst few rows:\n")
print(head(pData(gse110043)[, 1:15]))
