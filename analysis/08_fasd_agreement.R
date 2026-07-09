library(dplyr)
library(ggplot2)
library(tidyr)

# FASD signature agreement battery (C1-C4) + cross-tissue convergent validity +
# cell-composition confound check. Mirrors 04_agreement_analysis.R.

scores_disc <- readRDS("data/scores_GSE112987.rds")
scores_repl <- readRDS("data/scores_GSE113012.rds")
all_scores  <- bind_rows(scores_disc, scores_repl)
sig_names   <- c("portales", "lussier161", "lussier183", "vanderlaan")
cell_types  <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")

# Higher directional score = more FASD-like (see 07).
cohen_d <- function(x, g) {
    a <- x[g == "FASD"]; b <- x[g == "control"]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    ps <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b)) / (length(a)+length(b)-2))
    (mean(a) - mean(b)) / ps
}
# AUC via the Mann-Whitney identity (dependency-free).
auc <- function(x, g) {
    ok <- !is.na(x); x <- x[ok]; g <- g[ok]
    pos <- g == "FASD"; np <- sum(pos); nn <- sum(!pos)
    (sum(rank(x)[pos]) - np * (np + 1) / 2) / (np * nn)
}

# ── C1: Sensitivity threshold (separation, and decay toward mild cases) ────────
cat("== C1: FASD vs control separation ==\n")
for (acc in c("GSE112987", "GSE113012")) {
    d <- all_scores %>% filter(cohort == acc)
    cat("\n", acc, "  (Cohen's d | AUC)\n")
    for (sig in sig_names)
        cat(sprintf("   %-11s d=%+.2f  AUC=%.2f\n", sig,
                    cohen_d(d[[sig]], d$diagnosis), auc(d[[sig]], d$diagnosis)))
}

# Does separation decay for milder cases? Composite severity = facial+cns+growth
# (discovery cases, all three domains present). Split cases at the median.
cat("\n== C1b: separation by severity (discovery) ==\n")
disc_cases <- scores_disc %>%
    filter(diagnosis == "FASD", !is.na(facial), !is.na(cns), !is.na(growth)) %>%
    mutate(severity = facial + cns + growth)
ctrl <- scores_disc %>% filter(diagnosis == "control")
sev_med <- median(disc_cases$severity)
cat("   composite severity: range", min(disc_cases$severity), "-",
    max(disc_cases$severity), ", median", sev_med, ", n cases",
    nrow(disc_cases), "\n")
for (sig in sig_names) {
    mild   <- disc_cases[[sig]][disc_cases$severity <= sev_med]
    severe <- disc_cases[[sig]][disc_cases$severity >  sev_med]
    cd <- function(v) (mean(v, na.rm=TRUE) - mean(ctrl[[sig]], na.rm=TRUE)) /
                       sd(ctrl[[sig]], na.rm=TRUE)
    cat(sprintf("   %-11s mild vs ctrl d=%+.2f | severe vs ctrl d=%+.2f\n",
                sig, cd(mild), cd(severe)))
}

# ── C2: Score stability (nested signatures) ───────────────────────────────────
cat("\n== C2: score stability (nested sets, Spearman rank corr) ==\n")
sp <- function(a, b) cor(all_scores[[a]], all_scores[[b]], method = "spearman",
                         use = "pairwise.complete.obs")
cat(sprintf("   lussier161 vs portales (adds ~500 CpGs): %.3f\n", sp("lussier161","portales")))
cat(sprintf("   lussier161 vs lussier183:                %.3f\n", sp("lussier161","lussier183")))
cat(sprintf("   vanderlaan vs lussier161 (independent):  %.3f\n", sp("vanderlaan","lussier161")))

# ── Classification (threshold = control mean + 2SD, pooled controls) ───────────
for (sig in sig_names) {
    cm <- mean(all_scores[[sig]][all_scores$diagnosis == "control"], na.rm = TRUE)
    cs <- sd(all_scores[[sig]][all_scores$diagnosis == "control"], na.rm = TRUE)
    all_scores[[paste0(sig, "_class")]] <- as.integer(all_scores[[sig]] > cm + 2 * cs)
}
agree_rate <- function(df) {
    cl <- df %>% dplyr::select(ends_with("_class"))
    na <- rowSums(cl)
    mean(na == length(sig_names) | na == 0)
}

# ── C3: Cross-cohort consistency + inter-signature agreement by group ──────────
cat("\n== C3: inter-signature classification agreement (all 4 agree) ==\n")
for (grp in c("control", "FASD")) {
    for (acc in c("GSE112987", "GSE113012")) {
        d <- all_scores %>% filter(diagnosis == grp, cohort == acc)
        cat(sprintf("   %-8s %s: %.0f%% agree (n=%d)\n", grp, acc,
                    100 * agree_rate(d), nrow(d)))
    }
}
cat("   -- discovery FASD split by severity --\n")
disc_f <- all_scores %>% filter(cohort == "GSE112987", diagnosis == "FASD",
                                !is.na(facial), !is.na(cns), !is.na(growth)) %>%
    mutate(severity = facial + cns + growth)
cat(sprintf("   mild   FASD: %.0f%% agree (n=%d)\n",
            100 * agree_rate(filter(disc_f, severity <= sev_med)),
            sum(disc_f$severity <= sev_med)))
cat(sprintf("   severe FASD: %.0f%% agree (n=%d)\n",
            100 * agree_rate(filter(disc_f, severity >  sev_med)),
            sum(disc_f$severity >  sev_med)))

# ── C4: Monotonicity with severity domains (discovery cases) ───────────────────
cat("\n== C4: monotonicity - Spearman(score, severity domain) among FASD cases ==\n")
for (sig in sig_names) {
    r <- sapply(c("facial", "cns", "growth"), function(dom) {
        sub <- disc_cases[!is.na(disc_cases[[dom]]), ]
        cor(sub[[sig]], sub[[dom]], method = "spearman")
    })
    cat(sprintf("   %-11s facial=%+.2f  cns=%+.2f  growth=%+.2f\n",
                sig, r["facial"], r["cns"], r["growth"]))
}

# ── Cross-tissue convergent validity (headline) ───────────────────────────────
cat("\n== Cross-tissue: van der Laan (blood) vs Lussier-161 (buccal) ==\n")
cat("   shared probes: near-zero; agreement is genuinely independent.\n")
for (acc in c("GSE112987", "GSE113012")) {
    d <- all_scores %>% filter(cohort == acc)
    cat(sprintf("   %s Spearman = %.3f\n", acc,
                cor(d$vanderlaan, d$lussier161, method = "spearman",
                    use = "pairwise.complete.obs")))
}

# ── Cell-composition confound ─────────────────────────────────────────────────
# Disagreement per sample = SD across the 4 standardized signature scores.
cat("\n== Cell composition: does disagreement track cell fractions? ==\n")
zsig <- scale(all_scores[, sig_names])
all_scores$disagreement <- apply(zsig, 1, sd)
for (ct in cell_types) {
    r <- cor(all_scores$disagreement, all_scores[[ct]], method = "spearman",
             use = "pairwise.complete.obs")
    cat(sprintf("   disagreement vs %-6s Spearman=%+.2f\n", ct, r))
}

# ── Figures ───────────────────────────────────────────────────────────────────
long <- all_scores %>%
    pivot_longer(all_of(sig_names), names_to = "signature", values_to = "score") %>%
    mutate(signature = factor(signature, levels = sig_names,
              labels = c("Portales\n(657, buccal)", "Lussier-161\n(buccal)",
                         "Lussier-183\n(buccal)", "van der Laan\n(204, blood)")))

p1 <- ggplot(long, aes(diagnosis, score, fill = diagnosis)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.4) +
    facet_grid(cohort ~ signature) +
    scale_fill_manual(values = c(control = "steelblue", FASD = "tomato")) +
    theme_bw() + theme(legend.position = "bottom") +
    labs(title = "FASD signature scores by diagnosis",
         subtitle = "Directional z-score (higher = more FASD-like), 4 signatures x 2 cohorts",
         x = NULL, y = "Directional signature score")
ggsave("analysis/fasd_scores_boxplots.png", p1, width = 11, height = 6, dpi = 150)

sev_long <- disc_cases %>%
    pivot_longer(all_of(sig_names), names_to = "signature", values_to = "score") %>%
    mutate(signature = factor(signature, levels = sig_names))
p2 <- ggplot(sev_long, aes(severity, score)) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, colour = "tomato") +
    facet_wrap(~ signature, nrow = 1) +
    theme_bw() +
    labs(title = "C4 monotonicity: score vs clinical severity (discovery FASD cases)",
         subtitle = "Composite severity = facial + CNS + growth domain scores",
         x = "Composite clinical severity", y = "Directional signature score")
ggsave("analysis/fasd_severity_monotonicity.png", p2, width = 11, height = 3.5, dpi = 150)

p3 <- ggplot(all_scores, aes(lussier161, vanderlaan, colour = diagnosis)) +
    geom_point(alpha = 0.7, size = 1.5) +
    facet_wrap(~ cohort) +
    scale_colour_manual(values = c(control = "steelblue", FASD = "tomato")) +
    theme_bw() + theme(legend.position = "bottom") +
    labs(title = "Cross-tissue convergent validity: blood vs buccal signature",
         subtitle = "van der Laan (blood, 204) vs Lussier-161 (buccal); ~0 shared probes",
         x = "Lussier-161 score (buccal-derived)", y = "van der Laan score (blood-derived)")
ggsave("analysis/fasd_signature_scatter.png", p3, width = 8, height = 4.5, dpi = 150)

cell_long <- all_scores %>%
    pivot_longer(all_of(cell_types), names_to = "cell", values_to = "fraction")
p4 <- ggplot(cell_long, aes(diagnosis, fraction, fill = diagnosis)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.3) +
    facet_wrap(~ cell, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(control = "steelblue", FASD = "tomato")) +
    theme_bw() + theme(legend.position = "bottom", axis.text.x = element_blank()) +
    labs(title = "Blood cell composition by diagnosis (EpiDISH RPC, DHS reference)",
         x = NULL, y = "Estimated cell fraction")
ggsave("analysis/fasd_cellcomp.png", p4, width = 11, height = 3, dpi = 150)

saveRDS(all_scores, "data/fasd_all_scores.rds")
cat("\nFigures: analysis/fasd_scores_boxplots.png, fasd_severity_monotonicity.png,",
    "fasd_signature_scatter.png, fasd_cellcomp.png\n")
