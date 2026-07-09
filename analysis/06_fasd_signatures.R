library(GEOquery)
library(dplyr)

# Assemble the four published FASD signature CpG sets, extract their beta values
# from both cohorts, and build tidy phenotype tables. Mirrors 02_preprocess.R.
#
# Signature provenance and nesting structure: see signatures/README.md.
#   portales  (657, buccal)  - superset of both Lussier sets
#   lussier161 (161, buccal) - validated subset of portales
#   lussier183 (183, buccal) - predictor-influence subset of portales
#   vanderlaan (204, blood)  - independent, ~disjoint from the buccal sets

sig_dir <- "signatures"
read_sig <- function(f) unique(grep("^cg", readLines(file.path(sig_dir, f)), value = TRUE))

sigs <- list(
    portales   = read_sig("portales_casamar_2016_658.txt"),
    lussier161 = read_sig("lussier_2018_validated_161.txt"),
    lussier183 = read_sig("lussier_2018_predictor_183.txt"),
    vanderlaan = read_sig("vanderlaan_2025_episig_204.txt")
)
cat("Signature sizes:\n"); print(sapply(sigs, length))

all_cpgs <- unique(unlist(sigs))
cat("Total unique CpGs across signatures:", length(all_cpgs), "\n")

# ── Phenotype tidy-up ─────────────────────────────────────────────────────────
# FASD clinical severity domains (facial / cns / growth) are scored per case and
# NA for controls. Kept as integers for monotonicity (C4) analysis.
tidy_pheno <- function(ph, cohort) {
    num <- function(x) suppressWarnings(as.integer(as.character(x)))
    data.frame(
        geo_accession = ph$geo_accession,
        diagnosis = ifelse(ph[["disease state:ch1"]] == "FASD", "FASD", "control"),
        facial = num(ph[["facial:ch1"]]),
        cns    = num(ph[["cns:ch1"]]),
        growth = num(ph[["growth:ch1"]]),
        age    = suppressWarnings(as.numeric(ph[["age:ch1"]])),
        sex    = ph[["gender:ch1"]],
        cohort = cohort,
        stringsAsFactors = FALSE
    )
}

process_cohort <- function(acc, cohort) {
    g <- readRDS(file.path("data", paste0(acc, ".rds")))
    present <- intersect(all_cpgs, rownames(exprs(g)))
    beta <- exprs(g)[present, ]
    ph <- tidy_pheno(pData(g), cohort)
    cat("\n", acc, "-", length(present), "/", length(all_cpgs),
        "signature CpGs present;", ncol(beta), "samples\n")
    cat("   diagnosis:\n"); print(table(ph$diagnosis))
    list(beta = beta, pheno = ph)
}

d_disc <- process_cohort("GSE112987", "GSE112987")
d_repl <- process_cohort("GSE113012", "GSE113012")

# Coverage assertion: every signature should be fully represented on the array.
for (nm in names(sigs)) {
    cov <- mean(sigs[[nm]] %in% rownames(d_disc$beta))
    if (cov < 1) cat("WARNING: signature", nm, "coverage", round(100 * cov), "%\n")
}

saveRDS(d_disc, "data/GSE112987_fasd.rds")
saveRDS(d_repl, "data/GSE113012_fasd.rds")
saveRDS(sigs,   "data/fasd_signature_cpgs.rds")

cat("\nSaved data/GSE112987_fasd.rds, data/GSE113012_fasd.rds, data/fasd_signature_cpgs.rds\n")
cat("Preprocessing complete.\n")
