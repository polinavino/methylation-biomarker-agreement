library(dplyr)
library(ggplot2)
library(tidyr)

scores50660 <- readRDS("data/scores_GSE50660.rds")
scores42861 <- readRDS("data/scores_GSE42861.rds")
all_scores  <- bind_rows(scores50660, scores42861)

# ── C1: Sensitivity threshold ─────────────────────────────────────────────────
# Do signatures perform reliably at low exposure (former/occasional smokers)?
cat("── C1: Sensitivity threshold ──\n")

for (cohort in c("GSE50660", "GSE42861")) {
    cat("\n", cohort, "\n")
    d <- all_scores %>% filter(cohort == !!cohort, smoking %in% c("never","former","current"))
    for (sig in c("ahrr","epismoke","joehanes","epitob")) {
        means <- tapply(d[[sig]], d$smoking, mean, na.rm=TRUE)
        cat(sig, "- never:", round(means["never"],3),
            "former:", round(means["former"],3),
            "current:", round(means["current"],3), "\n")
    }
}

# ── C3: Cross-cohort consistency ──────────────────────────────────────────────
cat("\n── C3: Cross-cohort rank consistency ──\n")
# Classify each sample as smoker/non-smoker using each signature
# Threshold: mean of never smoker scores minus 2SD

classify <- function(scores, sig, threshold) {
    as.integer(scores[[sig]] < threshold)
}

for (sig in c("ahrr","epismoke","joehanes","epitob")) {
    never_mean <- mean(all_scores[all_scores$smoking=="never", sig], na.rm=TRUE)
    never_sd   <- sd(all_scores[all_scores$smoking=="never", sig], na.rm=TRUE)
    threshold  <- never_mean - 2 * never_sd
    all_scores[[paste0(sig,"_class")]] <- classify(all_scores, sig, threshold)
}

cat("\nClassification agreement (current smokers only):\n")
current <- all_scores %>% filter(smoking == "current")
agreement_current <- current %>%
    dplyr::select(ends_with("_class")) %>%
    mutate(n_agree = rowSums(.),
           all_agree = n_agree == 4 | n_agree == 0)
cat("Samples where all 4 signatures agree:", sum(agreement_current$all_agree),
    "of", nrow(agreement_current), "\n")

cat("\nClassification agreement (former smokers only):\n")
former <- all_scores %>% filter(smoking == "former")
agreement_former <- former %>%
    dplyr::select(ends_with("_class")) %>%
    mutate(n_agree = rowSums(.),
           all_agree = n_agree == 4 | n_agree == 0)
cat("Samples where all 4 signatures agree:", sum(agreement_former$all_agree),
    "of", nrow(agreement_former), "\n")

cat("\nClassification agreement (never smokers only):\n")
never <- all_scores %>% filter(smoking == "never")
agreement_never <- never %>%
    dplyr::select(ends_with("_class")) %>%
    mutate(n_agree = rowSums(.),
           all_agree = n_agree == 4 | n_agree == 0)
cat("Samples where all 4 signatures agree:", sum(agreement_never$all_agree),
    "of", nrow(agreement_never), "\n")

# ── C4: Monotonicity ──────────────────────────────────────────────────────────
cat("\n── C4: Monotonicity ──\n")
cat("Expected order: current < former < never (lower beta = more smoking)\n\n")
for (cohort in c("GSE50660","GSE42861")) {
    cat(cohort, ":\n")
    d <- all_scores %>%
        filter(cohort == !!cohort, smoking %in% c("never","former","current"))
    for (sig in c("ahrr","epismoke","joehanes","epitob")) {
        means <- tapply(d[[sig]], d$smoking, mean, na.rm=TRUE)
        monotone <- means["current"] < means["former"] & means["former"] < means["never"]
        cat(" ", sig, "- monotone:", monotone, "\n")
    }
}

# ── C2: Score stability ───────────────────────────────────────────────────────
cat("\n── C2: Score stability (AHRR alone vs multi-CpG) ──\n")
cat("Correlation between AHRR single-CpG and multi-CpG signatures:\n")
cat("AHRR vs EpiSmoke:", round(cor(all_scores$ahrr, all_scores$epismoke),3), "\n")
cat("AHRR vs Joehanes:", round(cor(all_scores$ahrr, all_scores$joehanes),3), "\n")
cat("AHRR vs EpiTob:  ", round(cor(all_scores$ahrr, all_scores$epitob),3), "\n")

# ── Plot ──────────────────────────────────────────────────────────────────────
plot_data <- all_scores %>%
    filter(smoking %in% c("never","former","current")) %>%
    pivot_longer(cols = c(ahrr, epismoke, joehanes, epitob),
                 names_to = "signature", values_to = "score") %>%
    mutate(smoking = factor(smoking, levels = c("never","former","current")),
           signature = factor(signature, 
                              levels = c("ahrr","epismoke","joehanes","epitob"),
                              labels = c("AHRR","EpiSmoke","Joehanes","EpiTob")))

p <- ggplot(plot_data, aes(x = smoking, y = score, fill = smoking)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_grid(cohort ~ signature) +
    scale_fill_manual(values = c("never"="steelblue","former"="orange","current"="tomato")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "Smoking methylation signature scores by smoking status",
         subtitle = "Four signatures across two independent cohorts",
         x = "Smoking status", y = "Mean beta value (lower = more methylation loss)",
         fill = "Smoking status")

ggsave("analysis/agreement_boxplots.png", p, width=12, height=7, dpi=150)

saveRDS(all_scores, "data/all_scores.rds")
cat("\nOutput: analysis/agreement_boxplots.png\n")
