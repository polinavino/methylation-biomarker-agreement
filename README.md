# Epigenetic Biomarker Agreement Analysis

A formal framework for evaluating whether epigenetic biomarker signatures agree
across independent cohorts, demonstrated on smoking-associated DNA methylation data.

## Motivation

Multiple epigenetic biomarker signatures have been proposed for the same
environmental exposures — smoking, alcohol, air pollution, and prenatal alcohol
exposure (FASD) — often derived from different cohorts using different statistical
approaches. When applied to independent datasets, these signatures frequently
disagree: a sample classified as "exposed" by one signature may be classified as
"unexposed" by another. This disagreement is rarely characterized formally.

This project asks: under what conditions do biomarker signatures agree, and is
their disagreement predictable? We propose four formal criteria for evaluating
biomarker agreement, analogous to the measurement desiderata developed for kinase
inhibitor selectivity definitions in a companion study
(ChemRxiv: https://doi.org/10.26434/chemrxiv.15001618/v1). We demonstrate the
framework on smoking-associated DNA methylation data, where multiple validated
signatures exist and independent cohort data is publicly available.

The longer-term goal is to apply this framework to fetal alcohol spectrum disorder
(FASD) biomarkers, where cross-cohort validation is limited by data access
constraints. This project builds the methodological toolkit and demonstrates
technical competence in methylation array analysis prior to seeking access to
restricted FASD datasets.

## Framework

We define four criteria for biomarker agreement:

**C1 — Sensitivity threshold:** Does the signature perform reliably below a
minimum exposure level? A well-formed signature should have a defined reliability
threshold below which classifications are not meaningful. Analog of the D1
reliability desideratum from the companion selectivity paper.

**C2 — Score stability:** Does adding weakly-associated CpGs change the
classification? A well-formed signature should not change its ranking of samples
when CpGs with weak association to the exposure are added. Signatures that are
sensitive to the inclusion of near-threshold CpGs are unstable under distributional
perturbation. Analog of D3.

**C3 — Cross-cohort consistency:** Does the signature rank samples the same way
across independent cohorts? Signatures derived in one population may not generalize
due to cohort-specific technical factors (array batch, cell type composition) or
biological factors (age, sex, ancestry). We quantify cross-cohort rank correlation
and characterize the samples that receive inconsistent classifications.

**C4 — Monotonicity:** Does increasing exposure monotonically increase the score?
A well-formed signature should not assign a higher score to a never-smoker than to
a heavy smoker. Violations of monotonicity indicate that the signature is capturing
something other than the intended exposure. Analog of D4.

## Signatures evaluated

We evaluate four published smoking methylation signatures:

| Signature | CpGs | Source |
|-----------|------|--------|
| AHRR single-CpG | 1 (cg05575921) | Breitling et al. 2011 |
| EpiSmoke | 4 (AHRR, 2q37.1, 6p21.33, F2RL3) | Shenker et al. 2013 |
| EpiTob | ~20 | Caramaschi et al. 2022 |
| Joehanes top-10 | 10 | Joehanes et al. 2016 |

## Datasets

| GEO accession | N | Platform | Smoking groups | Reference |
|---------------|---|----------|----------------|-----------|
| GSE50660 | 464 | 450k | Never (179) / Former (263) / Current (22) | Liu et al. 2013 |
| GSE42861 | 689 | 450k | Never (193) / Ex (228) / Current (200) / Occasional (66) | Liu et al. 2013 (RA study) |

## Methods

### Data processing
- Download via GEOquery
- Normalization: functional normalization (minfi)
- Cell type correction: Houseman method (estimateCellCounts)
- Batch correction: ComBat

### Analysis pipeline
1. `analysis/01_download_data.R` — download and inspect GEO datasets
2. `analysis/02_preprocess.R` — normalize, QC, cell type correction
3. `analysis/03_signature_scores.R` — compute all four signature scores per sample
4. `analysis/04_agreement_analysis.R` — apply C1-C4 framework
5. `analysis/05_plots.R` — generate figures

## Connection to companion work

This project is methodologically parallel to:

- **Kinase selectivity definitions paper** (ChemRxiv: https://doi.org/10.26434/chemrxiv.15001618/v1)
  — applies the same formal desiderata framework to binding-based selectivity metrics
- **RNA-seq pipeline** (github.com/polinavino/rnaseq-pipeline)
  — demonstrates transcriptional specificity analysis using the same framework

Together these projects develop a general approach to the problem of definitional
instability in biological measurement — applicable wherever multiple definitions
of the same biological concept coexist without a principled basis for choosing
among them.

## Results

### Signature correlations

All four signatures are highly correlated across both cohorts (r = 0.81-0.97),
indicating broad agreement. This contrasts with the kinase selectivity setting
where definitions clustered into two distinct families (r = 0.27-0.48 between
families). Smoking methylation signatures are more concordant than kinase
selectivity definitions, likely because they are all derived from overlapping
CpG sets anchored by AHRR.

### C1 — Sensitivity threshold

All signatures show the expected never > former > current gradient. AHRR has
the largest dynamic range (gap = 0.244 between never and current smokers) while
EpiTob has the smallest (gap = 0.081). Former smoker scores cluster much closer
to never smokers than to current smokers across all signatures, suggesting a
reliability threshold issue for former smoker classification.

### C3 — Cross-cohort consistency

| Smoking group | N | Signatures agree | Agreement rate |
|---------------|---|-----------------|----------------|
| Current smokers | 222 | 183 | 82% |
| Former smokers | 491 | 350 | 71% |
| Never smokers | 372 | 358 | 96% |

Former smokers show the highest rate of inconsistent classification (29%) —
nearly 1 in 3 former smokers is classified differently by different signatures.
This is the primary failure mode of existing smoking methylation biomarkers and
directly parallels the instability observed for intermediate-selectivity compounds
in the kinase selectivity paper.

### C4 — Monotonicity

All four signatures pass monotonicity in both cohorts — current smokers always
score lower than former smokers, who always score lower than never smokers.
No violations observed.

### C2 — Score stability

AHRR single-CpG correlates r = 0.826-0.922 with multi-CpG signatures. Adding
more CpGs changes rankings enough to produce disagreements, particularly in the
former smoker range where scores are most ambiguous.

### Key finding

Signatures agree on clear cases (never smokers 96%, current smokers 82%) but
diverge substantially on former smokers (71% agreement). This is a C1 violation
— below a certain cumulative exposure history, classification becomes unreliable
and signature-dependent. The disagreement is consistent across cohorts (not
batch-specific) and is driven by biological heterogeneity in methylation recovery
after smoking cessation. Different signatures weight that recovery differently
because they use different CpG sets with different rates of methylation reversal.

## Status

- [x] Repository initialized
- [x] GSE50660 downloaded (n=464, never/former/current smokers)
- [x] GSE42861 downloaded (n=689, rheumatoid arthritis cohort with smoking labels)
- [x] Preprocessing and signature CpG extraction
- [x] Signature score computation (AHRR, EpiSmoke, Joehanes, EpiTob)
- [x] Agreement analysis (C1-C4 framework)
- [x] Figures
- [ ] Extension to air pollution or stress methylation datasets
- [ ] Formal writeup
