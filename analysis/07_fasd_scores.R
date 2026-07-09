library(Biobase)
library(dplyr)
library(EpiDISH)

# Directional FASD signature scoring + blood cell-composition estimation.
# Mirrors 03_signature_scores.R, but FASD signatures are MIXED-direction
# (some CpGs hyper-, some hypomethylated in FASD), so a plain mean-beta score
# (valid for smoking) is invalid here. We use a directional z-score:
#
#     score(sample) = mean over signature CpGs of ( sign_i * z_i )
#       z_i    = within-cohort z-score of CpG i's beta
#       sign_i = +1 if hypermethylated in FASD, -1 if hypomethylated
#
# Higher score = more FASD-like. Per-CpG signs are estimated ONCE from the
# GSE112987 discovery cohort (FASD - control mean beta) and applied fixed to
# both cohorts, so the replication cohort (GSE113012) is scored with signs it
# did not inform -> replication is the honest out-of-sample test. Because signs
# are per-signature-consistent, between-signature agreement (C2/C3/cross-tissue)
# is invariant to the sign convention anyway.

d_disc <- readRDS("data/GSE112987_fasd.rds")
d_repl <- readRDS("data/GSE113012_fasd.rds")
sigs   <- readRDS("data/fasd_signature_cpgs.rds")

# ── Per-CpG direction, estimated on discovery ─────────────────────────────────
disc_dx <- d_disc$pheno$diagnosis[match(colnames(d_disc$beta),
                                        d_disc$pheno$geo_accession)]
mu_fasd <- rowMeans(d_disc$beta[, disc_dx == "FASD"],    na.rm = TRUE)
mu_ctrl <- rowMeans(d_disc$beta[, disc_dx == "control"], na.rm = TRUE)
cpg_sign <- sign(mu_fasd - mu_ctrl)            # +1 hyper / -1 hypo in FASD
cpg_sign[cpg_sign == 0] <- 1
cat("Per-CpG direction (discovery): ", sum(cpg_sign > 0), "hyper /",
    sum(cpg_sign < 0), "hypo\n")

# ── Directional scoring ───────────────────────────────────────────────────────
zscore_rows <- function(beta) {
    m <- rowMeans(beta, na.rm = TRUE)
    s <- apply(beta, 1, sd, na.rm = TRUE)
    s[s == 0 | is.na(s)] <- 1
    (beta - m) / s
}

compute_scores <- function(beta) {
    z <- zscore_rows(beta)
    zs <- z * cpg_sign[rownames(z)]            # orient toward FASD
    sapply(sigs, function(cpgs) {
        cpgs <- intersect(cpgs, rownames(zs))
        colMeans(zs[cpgs, , drop = FALSE], na.rm = TRUE)
    })
}

# ── Blood cell composition (EpiDISH, RPC, DHS blood reference) ─────────────────
data(centDHSbloodDMC.m)
estimate_cellcomp <- function(acc) {
    g <- readRDS(file.path("data", paste0(acc, ".rds")))
    beta_full <- exprs(g)
    ref_cpgs <- intersect(rownames(centDHSbloodDMC.m), rownames(beta_full))
    bm <- beta_full[ref_cpgs, ]
    # EpiDISH has no NA handling; impute the few NA reference cells with row means
    na_idx <- which(is.na(bm), arr.ind = TRUE)
    if (nrow(na_idx) > 0) bm[na_idx] <- rowMeans(bm, na.rm = TRUE)[na_idx[, 1]]
    out <- epidish(beta.m = bm, ref.m = centDHSbloodDMC.m[ref_cpgs, ], method = "RPC")
    as.data.frame(out$estF)
}

assemble <- function(d, acc) {
    scores <- as.data.frame(compute_scores(d$beta))
    scores$geo_accession <- rownames(scores)
    cell <- estimate_cellcomp(acc)
    cell$geo_accession <- rownames(cell)
    d$pheno %>%
        left_join(scores, by = "geo_accession") %>%
        left_join(cell,   by = "geo_accession")
}

scores_disc <- assemble(d_disc, "GSE112987")
scores_repl <- assemble(d_repl, "GSE113012")

saveRDS(scores_disc, "data/scores_GSE112987.rds")
saveRDS(scores_repl, "data/scores_GSE113012.rds")

# ── Sanity checks ─────────────────────────────────────────────────────────────
cat("\nMean directional score by diagnosis (should be FASD > control):\n")
for (acc in c("GSE112987", "GSE113012")) {
    s <- if (acc == "GSE112987") scores_disc else scores_repl
    cat(" ", acc, "\n")
    for (sig in names(sigs)) {
        mu <- tapply(s[[sig]], s$diagnosis, mean, na.rm = TRUE)
        cat(sprintf("    %-11s FASD=%+.3f  control=%+.3f  %s\n",
                    sig, mu["FASD"], mu["control"],
                    ifelse(mu["FASD"] > mu["control"], "ok", "INVERTED")))
    }
}

cat("\nInter-signature score correlation (discovery, Spearman):\n")
print(round(cor(scores_disc[, names(sigs)], method = "spearman",
                use = "pairwise.complete.obs"), 3))

cell_types <- colnames(centDHSbloodDMC.m)
cat("\nCell-fraction check (rows should sum ~1):\n")
cat("  cell types:", paste(cell_types, collapse = ", "), "\n")
cat("  discovery row sums range:",
    paste(round(range(rowSums(scores_disc[, cell_types])), 3), collapse = " to "), "\n")

cat("\nSaved data/scores_GSE112987.rds, data/scores_GSE113012.rds\n")
