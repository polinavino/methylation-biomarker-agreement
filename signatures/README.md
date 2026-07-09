# FASD DNA methylation signature CpG lists

One Illumina 450k probe ID (cg-ID) per line. These are the exposure/diagnosis
signatures scored by the FASD arm of the agreement analysis
(`analysis/06`–`08`). All four are 100% present on the 450k array.

| File | CpGs | Tissue derived | Source |
|------|------|----------------|--------|
| `portales_casamar_2016_658.txt` | 657 | buccal epithelium | Portales-Casamar et al. 2016, *Epigenetics & Chromatin* 9:25, Additional file 2 (Table S1, "DM Probes"). DOI 10.1186/s13072-016-0074-4 |
| `lussier_2018_validated_161.txt` | 161 | buccal epithelium | Lussier et al. 2018, *Clin Epigenetics* 10:5, Additional file 1, Suppl. Table 2 (CpGs validated in the KBHN cohort). DOI 10.1186/s13148-018-0439-6 |
| `lussier_2018_predictor_183.txt` | 183 | buccal epithelium | Lussier et al. 2018, Suppl. Table 3 (probes with non-zero influence in the gradient-boosting predictor). |
| `vanderlaan_2025_episig_204.txt` | 204 | **peripheral blood** | van der Laan et al. 2025, *Genet Med* 27(12):101586, Suppl. Table 2 (FAS episignature). DOI 10.1016/j.gim.2025.101586 |

## Structure that matters for the agreement analysis

- **Nested:** the Lussier validated-161 and predictor-183 sets are both 100%
  subsets of the Portales-Casamar set (they are its replicated / predictive
  fractions). Their mutual agreement is therefore partly guaranteed by
  construction — used here as the **C2 score-stability** test (does adding the
  ~500 extra Portales CpGs change sample rankings?).
- **Independent / cross-tissue:** the van der Laan blood episignature shares
  only 2 probes with Portales and 0 with Lussier-161. van der Laan (blood) vs
  Lussier (buccal), scored on the same blood samples, is a genuine
  **convergent-validity** test with no probe overlap inflating it.

## Direction (hyper/hypo in FASD)

These files contain probe IDs only. Per-CpG effect direction is estimated
empirically from the GSE112987 discovery cohort at scoring time (see
`analysis/07_fasd_scores.R`) rather than taken from the published buccal
Δβ columns, because buccal-derived directions may not transfer to blood.
Published-direction concordance is noted as a future refinement (requires an
xlsx reader to re-parse the supplementary tables).
