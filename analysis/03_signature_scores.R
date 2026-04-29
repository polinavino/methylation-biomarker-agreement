library(dplyr)
library(tibble)

# ── Load data ─────────────────────────────────────────────────────────────────
d50660 <- readRDS("data/GSE50660_processed.rds")
d42861 <- readRDS("data/GSE42861_processed.rds")
sigs   <- readRDS("data/signature_cpgs.rds")

compute_scores <- function(beta, pheno, sigs) {

    # AHRR single CpG — lower beta = more likely smoker
    ahrr_score <- beta[sigs$ahrr, ]

    # EpiSmoke — mean of 4 CpGs (all hypomethylated in smokers)
    # lower score = more likely smoker
    epismoke_score <- colMeans(beta[sigs$epismoke, ])

    # Joehanes — weighted sum using effect size direction
    # all top CpGs are hypomethylated in smokers so mean works
    joehanes_score <- colMeans(beta[sigs$joehanes, ])

    # EpiTob — mean of signature CpGs
    epitob_score <- colMeans(beta[sigs$epitob, ])

    scores <- data.frame(
        sample        = colnames(beta),
        ahrr          = as.numeric(ahrr_score),
        epismoke      = epismoke_score,
        joehanes      = joehanes_score,
        epitob        = epitob_score,
        stringsAsFactors = FALSE
    ) %>%
        left_join(pheno, by = c("sample" = "geo_accession"))

    return(scores)
}

scores50660 <- compute_scores(d50660$beta, d50660$pheno, sigs)
scores42861 <- compute_scores(d42861$beta, d42861$pheno, sigs)

cat("GSE50660 scores (first 5 rows):\n")
print(head(scores50660[, c("sample","ahrr","epismoke","joehanes","epitob","smoking")], 5))

cat("\nGSE42861 scores (first 5 rows):\n")
print(head(scores42861[, c("sample","ahrr","epismoke","joehanes","epitob","smoking")], 5))

saveRDS(scores50660, "data/scores_GSE50660.rds")
saveRDS(scores42861, "data/scores_GSE42861.rds")

cat("\nCorrelations between signatures (GSE50660):\n")
print(round(cor(scores50660[, c("ahrr","epismoke","joehanes","epitob")]), 3))

cat("\nCorrelations between signatures (GSE42861):\n")
print(round(cor(scores42861[, c("ahrr","epismoke","joehanes","epitob")]), 3))

cat("\nDone.\n")
