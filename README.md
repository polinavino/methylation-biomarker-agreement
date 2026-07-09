# Epigenetic Biomarker Agreement Analysis

A formal framework for evaluating whether epigenetic biomarker signatures agree
across independent cohorts, demonstrated on smoking-associated DNA methylation data
and extended to a reliability audit of fetal alcohol spectrum disorder (FASD)
signatures.

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

The framework is then applied to fetal alcohol spectrum disorder (FASD)
biomarkers (see "FASD extension" below), using the public blood cohorts and
signature CpG lists that are now available. Extension to the larger *restricted*
FASD datasets (e.g. the gated clinical episignature cohorts) remains future work;
this project establishes the methodology and the public-data audit first.

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

**C3 — Cross-cohort / cross-signature consistency:** Do independent signatures
classify the same samples the same way, and does that agreement hold across
independent cohorts? Signatures derived in one population may not generalize due
to cohort-specific technical factors (array batch, cell type composition) or
biological factors (age, sex, ancestry). Concretely, we classify each sample with
every signature and quantify the fraction that all signatures classify
identically, within each exposure/diagnosis group, reported across the independent
cohorts (and, for FASD, across tissues). This is an inter-signature agreement rate,
not a between-cohort correlation of matched samples (the cohorts contain different
individuals).

**C4 — Monotonicity:** Does increasing exposure monotonically increase the score?
A well-formed signature should not assign a higher score to a never-smoker than to
a heavy smoker. Violations of monotonicity indicate that the signature is capturing
something other than the intended exposure. Analog of D4.

## Signatures evaluated

We evaluate four published smoking methylation signatures:

| Signature | CpGs used here | Source |
|-----------|------|--------|
| AHRR single-CpG | 1 (cg05575921) | Breitling et al. 2011 |
| EpiSmoke | 4 (AHRR, 2q37.1, 6p21.33, F2RL3) | Shenker et al. 2013 |
| Joehanes | 9 (top-ranked CpG subset) | Joehanes et al. 2016 |
| EpiTob | 12 (CpG subset) | Caramaschi et al. 2022 |

The exact CpG lists are defined in [`analysis/02_preprocess.R`](analysis/02_preprocess.R);
their union is 12 CpGs (EpiTob is a superset of the others). Joehanes and EpiTob
are scored here on top-ranked CpG subsets, not their full published panels.

## Datasets

| GEO accession | N | Platform | Smoking groups | Reference |
|---------------|---|----------|----------------|-----------|
| GSE50660 | 464 | 450k | Never (179) / Former (263) / Current (22) | Liu et al. 2013 |
| GSE42861 | 689 | 450k | Never (193) / Ex (228) / Current (200) / Occasional (66) | Liu et al. 2013 (RA study) |

## Methods

### Data processing
- Download via GEOquery.
- Beta values are taken **as-published** from each GEO series matrix (already
  normalized by the original authors). This analysis does **not** re-normalize
  from IDATs, and does **not** apply ComBat batch correction — an agreement audit
  deliberately scores the signatures on the data as downstream users would obtain
  it. Cross-cohort robustness (C3) is assessed rather than corrected away.
- Cell type: the smoking arm applies no cell-type adjustment; the FASD arm
  *estimates* blood cell composition with EpiDISH (RPC, DHS reference) and uses it
  as a **confound check** (does signature disagreement track cell fractions?),
  not as a correction.
- This is a known limitation, discussed in `docs/writeup.md`.

### Analysis pipeline
Smoking demonstration (01–04):
1. `analysis/01_download_data.R` — download and inspect GEO datasets
2. `analysis/02_preprocess.R` — extract signature CpGs, tidy phenotype
3. `analysis/03_signature_scores.R` — compute all four signature scores per sample
4. `analysis/04_agreement_analysis.R` — apply C1–C4 framework and generate figures

FASD extension (05–08, see "FASD extension" below):
5. `analysis/05_fasd_download.R` — download the two FASD blood cohorts
6. `analysis/06_fasd_signatures.R` — assemble the four FASD signatures, extract CpGs
7. `analysis/07_fasd_scores.R` — directional scoring + EpiDISH cell composition
8. `analysis/08_fasd_agreement.R` — C1–C4 + cross-tissue agreement + confound checks

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

## Results (smoking demonstration)

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

## FASD extension (reliability audit)

The framework's intended target is FASD, where multiple signatures exist but
cross-cohort validation is limited by data access. Enough is now public to run
an honest reliability audit — deliberately framed as an audit of limits, not as
a prenatal screen (the evidence does not support the latter; see below).

### Datasets (peripheral blood, 450k; Henneman group, Amsterdam AMC)

| GEO accession | N | FASD / control | Role | Reference |
|---------------|---|----------------|------|-----------|
| GSE112987 | 103 | 39 / 64 | Discovery | PMID 30873861 |
| GSE113012 | 35 | 7 / 28 | Replication | PMID 30873861 |

Each FASD case carries the three clinical diagnostic domains (facial, CNS,
growth, scored 1–4), enabling severity-stratified analysis.

### Signatures evaluated

| Signature | CpGs | Tissue derived | Reference |
|-----------|------|----------------|-----------|
| Portales-Casamar | 657 | buccal | Portales-Casamar et al. 2016 |
| Lussier validated | 161 | buccal | Lussier et al. 2018 |
| Lussier predictor | 183 | buccal | Lussier et al. 2018 |
| van der Laan episignature | 204 | **blood** | van der Laan et al. 2025 |

Two comparison regimes fall out of the overlap structure: the Lussier sets are
100% nested inside Portales (their agreement is partly guaranteed — used as the
C2 test), while van der Laan shares ~0 probes with the buccal sets (a genuine,
independent cross-tissue convergent-validity test). CpG lists are vendored in
[`signatures/`](signatures/).

### Methods note

FASD signatures are **mixed-direction** (some CpGs hyper-, some hypomethylated
in FASD), so the smoking pipeline's plain mean-beta score is invalid. Scores are
directional z-scores (`sign_i · z_i`, higher = more FASD-like); per-CpG signs are
estimated on the discovery cohort and applied fixed to replication, so the
replication cohort is scored out-of-sample. Blood cell composition is estimated
with EpiDISH (RPC, DHS blood reference) as a confound check.

### Results

- **C1 (separation).** The blood-native van der Laan signature separates FASD
  from control strongly (AUC 0.93 discovery / 0.96 replication; Cohen's d 2.6 /
  3.2) — but see the non-independence caveat below: this is likely an upper bound.
  The buccal-derived signatures (independent of these cohorts) transfer poorly to
  blood (AUC 0.68–0.79, d 0.7–1.2). Separation is weaker for milder cases across
  all signatures.
- **C2 (stability).** Nested signatures agree almost perfectly (Spearman 0.96),
  but this is largely by construction; the independent cross-tissue pair agrees
  only moderately (0.38). Shared-probe overlap should not be mistaken for
  independent validation.
- **C3 (agreement).** The four signatures agree on controls (98% / 82% of samples
  classified identically) but diverge sharply on cases (36% / 0%). They agree on
  who is *not* FASD and disagree on who *is*.
- **C4 (monotonicity).** van der Laan tracks clinical severity (Spearman 0.42–0.54
  with facial/CNS/growth); the buccal signatures are essentially flat (0.02–0.26).
- **Cross-tissue convergent validity.** van der Laan (blood) vs Lussier (buccal),
  with near-zero shared probes, correlate only 0.43 (discovery) / 0.16
  (replication) despite both targeting FASD.
- **Confound check.** Inter-signature disagreement barely tracks cell fractions
  (|Spearman| ≤ 0.14), so the disagreement is not a cell-composition artifact —
  it reflects poor cross-tissue transfer. (FASD cases do differ in composition —
  higher B/CD4T, lower neutrophils — a separate caveat for the signal itself.)

### Caveat: the van der Laan signature is not independent of these cohorts

The corresponding author of van der Laan et al. 2025 (Peter Henneman, Amsterdam
UMC) is also the depositor of GSE112987/GSE113012 (PMID 30873861) — the same lab,
same tissue (blood 450K), overlapping timeframe. Its 93-case FAS cohort is the
natural superset of the earlier Amsterdam blood FASD samples. Per-sample reuse
could not be confirmed from the (access-restricted) main text, but **sample
overlap cannot be excluded**, so van der Laan's strong performance here is likely
**partly circular** — read it as an upper bound, not clean out-of-sample
validation. The two findings that do *not* depend on it stand regardless: the
buccal signatures are genuinely independent of these cohorts and still transfer
poorly to blood, and the between-signature agreement results (C2/C3, cross-tissue)
are invariant to which signature scores best.

### Key finding

Unlike the smoking case (where disagreement concentrated in an intermediate
within-tissue zone), the dominant FASD failure mode is **poor cross-tissue
transfer**: buccal-derived signatures carry only a weak, above-chance signal in
blood (AUC 0.68–0.79) and do not track severity there, while the blood-native
episignature separates cases strongly (AUC 0.93–0.96) and is monotone with
severity. This is a reliability-audit result — it says which signature to trust
in which tissue, not that a general prenatal screen is feasible.

### Responsible-use framing

This audit concerns **molecular classification in already-diagnosed children**
(postnatal blood), which is tractable. It does **not** support a prenatal or
at-birth exposure screen: the PACE cord-blood meta-analysis (Sharp et al. 2018)
found no replicating CpGs for maternal alcohol, and such a screen would carry
serious false-positive and stigma risks. See `docs/writeup.md` for the full
discussion, including the non-independence caveat above (the van der Laan
signature shares a lab, tissue, and likely samples with this Amsterdam dataset).

## Status

- [x] Repository initialized
- [x] GSE50660 downloaded (n=464, never/former/current smokers)
- [x] GSE42861 downloaded (n=689, rheumatoid arthritis cohort with smoking labels)
- [x] Preprocessing and signature CpG extraction
- [x] Signature score computation (AHRR, EpiSmoke, Joehanes, EpiTob)
- [x] Agreement analysis (C1-C4 framework)
- [x] Figures
- [x] FASD extension: blood cohorts (GSE112987, GSE113012) downloaded
- [x] FASD signatures assembled (Portales-Casamar, Lussier ×2, van der Laan)
- [x] FASD directional scoring + EpiDISH cell composition
- [x] FASD agreement analysis (C1-C4 + cross-tissue + confound checks) and figures
- [x] Formal writeup (smoking + FASD)
- [x] Checked van der Laan / Amsterdam-cohort independence — NOT independent
      (same lab and corresponding author, likely sample overlap; van der Laan
      performance treated as an upper bound, see caveat above)
- [ ] Extension to air pollution or stress methylation datasets
- [ ] (Optional) obtain van der Laan per-sample list to quantify exact overlap
- [ ] (Optional) published-direction concordance (re-parse supplementary Δβ tables)
