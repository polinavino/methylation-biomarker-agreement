library(GEOquery)
library(dplyr)

ahrr_cpg <- c("cg05575921")

epismoke_cpgs <- c("cg05575921",
                   "cg21566642",
                   "cg04180046",
                   "cg03636183")

joehanes_cpgs <- c("cg05575921",
                   "cg03636183",
                   "cg21566642",
                   "cg06126421",
                   "cg23387569",
                   "cg01940273",
                   "cg19859270",
                   "cg12803068",
                   "cg04551776")

epitob_cpgs <- c("cg05575921",
                 "cg03636183",
                 "cg19859270",
                 "cg23387569",
                 "cg06126421",
                 "cg12803068",
                 "cg01940273",
                 "cg04551776",
                 "cg09935388",
                 "cg18316974",
                 "cg04180046",
                 "cg21566642")

all_cpgs <- unique(c(ahrr_cpg, epismoke_cpgs, joehanes_cpgs, epitob_cpgs))
cat("Total unique CpGs needed:", length(all_cpgs), "\n")

cat("Loading GSE50660...\n")
gse50660 <- readRDS("data/GSE50660.rds")
beta50660 <- exprs(gse50660)[all_cpgs, ]

cat("Loading GSE42861...\n")
gse42861 <- readRDS("data/GSE42861.rds")
beta42861 <- exprs(gse42861)[all_cpgs, ]

cat("Beta dimensions GSE50660:", dim(beta50660), "\n")
cat("Beta dimensions GSE42861:", dim(beta42861), "\n")

pheno50660 <- pData(gse50660) %>%
    dplyr::select(geo_accession,
                  smoking_raw = `smoking (0, 1 and 2, which represent never, former and current smokers):ch1`) %>%
    mutate(smoking_raw = as.integer(smoking_raw),
           smoking = case_when(
               smoking_raw == 0 ~ "never",
               smoking_raw == 1 ~ "former",
               smoking_raw == 2 ~ "current"),
           cohort = "GSE50660")

pheno42861 <- pData(gse42861) %>%
    dplyr::select(geo_accession,
                  smoking = `smoking status:ch1`) %>%
    mutate(smoking = case_when(
               smoking == "never"      ~ "never",
               smoking == "ex"         ~ "former",
               smoking == "current"    ~ "current",
               smoking == "occasional" ~ "occasional",
               TRUE ~ NA_character_),
           cohort = "GSE42861")

cat("\nGSE50660 smoking:\n")
print(table(pheno50660$smoking))
cat("\nGSE42861 smoking:\n")
print(table(pheno42861$smoking, useNA="ifany"))

saveRDS(list(beta=beta50660, pheno=pheno50660), "data/GSE50660_processed.rds")
saveRDS(list(beta=beta42861, pheno=pheno42861), "data/GSE42861_processed.rds")

saveRDS(list(
    ahrr     = ahrr_cpg,
    epismoke = epismoke_cpgs,
    joehanes = joehanes_cpgs,
    epitob   = epitob_cpgs
), "data/signature_cpgs.rds")

cat("\nPreprocessing complete.\n")
