library(GEOquery)
library(dplyr)

# FASD peripheral-blood 450k methylation cohorts (Henneman group, Amsterdam AMC).
# SuperSeries GSE113018; primary publication PMID 30873861.
#   GSE112987 - discovery,   n=103 (39 FASD / 64 control)
#   GSE113012 - replication, n= 35 ( 7 FASD / 28 control)
# Both expose a clean beta matrix via exprs(), same access pattern as the
# smoking cohorts in 01_download_data.R.

outdir <- "data"
dir.create(outdir, showWarnings = FALSE)

download_fasd <- function(acc) {
    cat("Downloading", acc, "...\n")
    g <- getGEO(acc, GSEMatrix = TRUE, destdir = outdir)[[1]]
    ph <- pData(g)
    e  <- exprs(g)

    cat("  ", acc, ":", nrow(g), "CpGs x", ncol(g), "samples\n")
    cat("   beta range:", paste(round(range(e, na.rm = TRUE), 3), collapse = " to "),
        "| NA:", sum(is.na(e)), "\n")

    # Disease state and FASD clinical severity domains (facial / cns / growth),
    # scored per case; NA for controls.
    for (col in c("disease state:ch1", "facial:ch1", "cns:ch1",
                  "growth:ch1", "gender:ch1")) {
        if (col %in% colnames(ph)) {
            cat("  ", col, ":\n")
            print(table(ph[[col]], useNA = "ifany"))
        }
    }

    saveRDS(g, file.path(outdir, paste0(acc, ".rds")))
    cat("   Saved", file.path(outdir, paste0(acc, ".rds")), "\n\n")
    invisible(g)
}

download_fasd("GSE112987")   # discovery
download_fasd("GSE113012")   # replication (may already be cached)

cat("Download complete.\n")
