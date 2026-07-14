# A formal framework for evaluating epigenetic biomarker agreement: demonstration on smoking-associated DNA methylation

**Polina Vinogradova**

---

## Abstract

Multiple epigenetic biomarker signatures have been proposed for the same
environmental exposures, yet there is no formal framework for evaluating
whether they agree and under what conditions they diverge. We propose four
criteria — sensitivity threshold (C1), score stability (C2), cross-cohort
consistency (C3), and monotonicity (C4) — motivated by analogous desiderata
developed for kinase inhibitor selectivity definitions. We apply this framework
to four published smoking-associated DNA methylation signatures (AHRR, EpiSmoke,
Joehanes, EpiTob) across two independent whole-blood 450k cohorts (GSE50660,
n=464; GSE42861, n=689). Signatures are highly concordant for current smokers
(82% agreement) and never smokers (96%) but diverge substantially for former
smokers (71% agreement). All four signatures satisfy monotonicity (C4) but
show sensitivity threshold violations (C1) at low cumulative exposure. The
divergence is consistent across cohorts and driven by biological heterogeneity
in methylation recovery after smoking cessation. This framework is intended as
a generalizable tool for evaluating biomarker agreement in settings where
ground truth is unavailable. We then apply it to fetal alcohol spectrum disorder
signatures (Section 6), where the dominant failure mode is poor cross-tissue
transfer: buccal-derived signatures carry only a weak signal in blood and do not
track severity there, while a blood-native episignature does both.

---

## 1. Introduction

Environmental exposures leave epigenetic signatures in blood DNA methylation
that persist long after exposure ends. For tobacco smoking, dozens of
differentially methylated CpG sites have been identified in large cohort studies,
and multiple composite biomarker signatures have been developed for classifying
individuals by smoking status. These signatures — including single-CpG scores
based on AHRR methylation, multi-CpG indices such as EpiSmoke and EpiTob, and
large-scale meta-analysis-derived CpG sets — are routinely applied to independent
cohorts for exposure classification and epidemiological adjustment.

However, when applied to the same dataset, these signatures do not always agree.
A sample classified as a former smoker by one signature may be classified as a
never smoker by another. This disagreement is rarely characterized formally:
existing comparisons focus on overall predictive accuracy rather than on the
structural conditions under which signatures diverge.

This problem is not unique to smoking methylation. In kinase inhibitor
pharmacology, multiple selectivity definitions — the S-score, selectivity
entropy, Gini coefficient, and ratio-based measures — are applied
interchangeably to characterize drug-kinase binding profiles, yet they measure
fundamentally different properties and disagree systematically for compounds with
intermediate binding profiles (Vinogradova, 2025). The structural similarity
between these two problems — multiple definitions of the same biological concept,
each with different mathematical properties, producing inconsistent rankings —
motivates a common analytical framework.

Here we propose four formal criteria for evaluating biomarker agreement,
directly analogous to the desiderata developed for kinase selectivity definitions,
and demonstrate their application to smoking methylation data.

---

## 2. Framework

Let S = {s_1, ..., s_k} be a set of biomarker signatures, each mapping a
methylation profile M to a scalar score. We define four criteria:

**C1 — Sensitivity threshold.** A well-formed signature should have a defined
exposure level below which classification is not reliable. Specifically, the
score distributions of exposed and unexposed groups should be sufficiently
separated that a classification threshold can be set with defined sensitivity
and specificity. Signatures that fail to separate groups at low cumulative
exposure violate C1.

**C2 — Score stability.** Adding a CpG with weak association to the exposure
should not substantially change the ranking of samples. Formally, for any CpG
site c with effect size below a threshold delta, the rank correlation between
scores computed with and without c should exceed a stability bound rho.

**C3 — Cross-cohort consistency.** Signatures should rank samples consistently
across independent cohorts after adjustment for known technical confounders.
Formally, the Spearman rank correlation between signature scores in two cohorts,
computed on matched samples or on group means, should exceed a reproducibility
bound.

**C4 — Monotonicity.** Increasing cumulative exposure should monotonically
increase (or decrease) the signature score. Formally, if exposure group A has
higher cumulative exposure than group B, the mean score for A should be
directionally consistent across all signatures.

---

## 3. Methods

### 3.1 Datasets

We used two publicly available whole-blood 450k methylation datasets:

- **GSE50660** (Liu et al., 2013): 464 samples, never (n=179), former (n=263),
  current (n=22) smokers.
- **GSE42861** (Liu et al., 2013): 689 samples from a rheumatoid arthritis
  study with smoking status available: never (n=193), former/ex (n=228),
  current (n=200), occasional (n=66).

Beta values were extracted directly from GEO series matrix files using GEOquery.
The union of CpG sites across the four signatures was retained (n=12 CpGs; EpiTob
is a superset of the other three).

### 3.2 Signatures

Four smoking methylation signatures were evaluated:

| Signature | CpGs | Reference |
|-----------|------|-----------|
| AHRR | 1 (cg05575921) | Breitling et al. 2011 |
| EpiSmoke | 4 | Shenker et al. 2013 |
| Joehanes | 9 | Joehanes et al. 2016 |
| EpiTob | 12 | Caramaschi et al. 2022 |

Scores were computed as the mean beta value across signature CpGs (lower =
more methylation loss = more likely smoker).

### 3.3 Agreement analysis

For C1, we computed mean scores by smoking group and assessed separation.
For C2, we computed pairwise Spearman correlations between signatures.
For C3, we classified each sample as smoker/non-smoker using a threshold of
mean(never) - 2*SD(never) and computed the fraction of samples where all four
signatures agreed.
For C4, we tested whether group means followed the expected never > former >
current ordering.

---

## 4. Results

### 4.1 Signature correlations

All four signatures are highly correlated across both cohorts (r = 0.81-0.97),
indicating broad agreement. Joehanes and EpiTob are nearly identical (r = 0.97),
reflecting their shared CpG composition. AHRR shows the lowest correlation with
EpiTob (r = 0.83), suggesting that the single-CpG and multi-CpG approaches
capture partially distinct aspects of smoking-related methylation change.

### 4.2 C1 — Sensitivity threshold

All signatures discriminate never from current smokers with large effect sizes
(AHRR gap = 0.244 beta units; EpiTob gap = 0.081 beta units). However, former
smokers score intermediate to never and current in all signatures, with
substantial overlap with both groups. The AHRR single-CpG signature shows the
largest separation between groups (C1 best satisfied); EpiTob shows the smallest
separation (C1 least satisfied).

### 4.3 C3 — Cross-cohort consistency

| Smoking group | N | Agreement | Rate |
|---------------|---|-----------|------|
| Current smokers | 222 | 183/222 | 82% |
| Former smokers | 491 | 350/491 | 71% |
| Never smokers | 372 | 358/372 | 96% |

Signatures agree strongly on never smokers (96%) and moderately on current
smokers (82%), but diverge substantially for former smokers (71% agreement,
meaning 29% of former smokers receive inconsistent classifications). This
pattern is consistent across both cohorts, indicating it reflects biological
heterogeneity in methylation recovery rather than cohort-specific technical
factors.

### 4.4 C4 — Monotonicity

All four signatures satisfy monotonicity in both cohorts: current < former 
never in mean score. No violations observed.

### 4.5 C2 — Score stability

Pairwise correlations between signatures range from r = 0.83 to r = 0.97.
The lower bound (AHRR vs EpiTob, r = 0.83) indicates that adding weakly
associated CpGs to the AHRR single-site score changes sample rankings
meaningfully, particularly in the intermediate score range corresponding to
former smokers.

---

## 5. Discussion

The primary finding is that smoking methylation signatures agree well on
unambiguous cases — current heavy smokers and never smokers — but diverge
substantially for former smokers. This is a C1 violation: all four signatures
fail to reliably classify individuals whose cumulative exposure history places
them in a biologically ambiguous zone. The divergence is not random; it is
concentrated in the former smoker group and is consistent across independent
cohorts.

This pattern directly parallels the findings of the companion kinase selectivity
paper (Vinogradova, 2025), where selectivity definitions agreed on compounds
with extreme binding profiles (very selective or very promiscuous) but diverged
for compounds with intermediate profiles. In both cases, definitional instability
is concentrated in the biologically most interesting region — the intermediate
zone where classification is most consequential.

The practical implication for epidemiology is that studies adjusting for
smoking status using methylation-based classification should be cautious about
former smokers. Up to 29% of former smokers may be misclassified differently
depending on which signature is used, introducing systematic error that varies
by choice of biomarker.

**Limitation:** This analysis uses pre-processed beta values from GEO series
matrix files without re-normalization or cell-type correction. Cell type
composition is a known confounder of blood methylation and may contribute to
disagreement among signatures. Full re-analysis from IDAT files with cell-type
deconvolution is planned.

**Future directions:** The framework is applied to fetal alcohol spectrum
disorder (FASD) methylation signatures in Section 6. The smoking demonstration
establishes the technical methodology and provides a within-tissue benchmark
against which the FASD (cross-tissue) results are interpreted.

---

## 6. FASD extension: a reliability audit

### 6.1 Motivation and framing

FASD is the framework's intended target. Multiple blood/buccal DNAm signatures
now exist, but cross-cohort validation has been limited by data access. Enough is
now public to run a reliability audit. We frame it deliberately as an audit of
*limits*: it concerns molecular classification in **already-diagnosed children**
(postnatal blood), which is tractable, and it makes no claim toward a prenatal or
at-birth exposure screen. The latter is not supported by evidence — the PACE
cord-blood meta-analysis (Sharp et al., 2018; 1,147 exposed / 1,928 controls)
found no CpGs surviving multiple-testing correction for maternal alcohol — and
would carry serious false-positive and stigma/child-protection risks.

### 6.2 Data and signatures

Two peripheral-blood 450k cohorts from the Amsterdam AMC group (SuperSeries
GSE113018; PMID 30873861): **GSE112987** (discovery, n=103; 39 FASD / 64 control)
and **GSE113012** (replication, n=35; 7 FASD / 28 control). Each FASD case carries
the three clinical diagnostic domains (facial, CNS, growth; scored 1–4).

Four published signatures were scored: Portales-Casamar (657 CpGs, buccal),
Lussier validated (161, buccal), Lussier predictor (183, buccal), and van der
Laan (204, **blood** episignature). The Lussier sets are 100% nested within
Portales (agreement partly by construction → the C2 test); van der Laan shares
~0 probes with the buccal sets (a genuine cross-tissue convergent-validity test).

### 6.3 Methods

FASD signatures are mixed-direction, so the smoking mean-beta score is invalid.
We used a directional z-score, `score = mean_i(sign_i · z_i)`, higher = more
FASD-like, with per-CpG signs estimated on the discovery cohort and applied fixed
to replication (so replication is scored out-of-sample). Blood cell composition
was estimated with EpiDISH (RPC, DHS blood reference) as a confound check. AUC
was computed via the Mann-Whitney identity.

### 6.4 Results

- **C1 (separation).** The blood-native van der Laan signature separates FASD from
  control strongly (AUC 0.93 discovery / 0.96 replication; d 2.6 / 3.2) — but this
  is **not independent validation** (confirmed same-individual overlap with the
  2025 derivation cohort; see caveat iii), so read it as an upper bound. Buccal
  signatures (genuinely independent of these cohorts) transfer poorly to blood
  (AUC 0.68–0.79; d 0.7–1.2). Separation is weaker for milder cases across all
  signatures.
- **C2 (stability).** Nested signatures agree almost perfectly (Spearman 0.96) —
  largely by construction — whereas the independent cross-tissue pair agrees only
  moderately (0.38).
- **C3 (consistency).** Signatures classify controls concordantly (98% / 82%) but
  diverge on cases (36% / 0% of cases classified identically by all four).
- **C4 (monotonicity).** van der Laan tracks clinical severity (Spearman 0.42–0.54
  with facial/CNS/growth); buccal signatures are flat (0.02–0.26).
- **Cross-tissue convergent validity.** van der Laan vs Lussier (near-zero shared
  probes) correlate only 0.43 / 0.16, despite both targeting FASD.
- **Confound check.** Inter-signature disagreement barely tracks cell fractions
  (|Spearman| ≤ 0.14), so it is not a cell-composition artifact. (FASD cases do
  differ in composition — higher B/CD4T, lower neutrophils — a caveat for the
  signal itself, distinct from the disagreement.)

### 6.5 Interpretation

Where the smoking analysis found disagreement concentrated in an intermediate
*within-tissue* zone (former smokers), the dominant FASD failure mode is
**poor cross-tissue transfer**: buccal-derived signatures carry only a weak,
above-chance signal in blood (AUC 0.68–0.79) and do not track its severity there,
while the blood-native episignature separates cases strongly and is monotone with
severity. The practical read-out is not "FASD DNAm classification fails" but
"signature validity is tissue-bound" — a signature must be applied in, or
explicitly re-validated for, the tissue it is used on.

**Caveats.** (i) The replication cohort is small (7 FASD); effect sizes should be
read with the discovery cohort (n=39) carrying the weight. (ii) Per-CpG
directions were estimated empirically (blood) rather than from the published
buccal Δβ tables; published-direction concordance is a future refinement.
(iii) **Non-independence of the van der Laan signature (important; confirmed).** The
corresponding author of van der Laan et al. (2025), Peter Henneman (Amsterdam
UMC), is also the depositor of GSE112987/GSE113012 (PMID 30873861). **Overlap is
now confirmed by the corresponding author** (P. Henneman, personal communication,
July 2026): the earlier 450K cohort deposited as GSE112987/GSE113012 (published in
*Epigenomics* 2019) was included in the cohort used to derive the 2025
episignature, and a substantial number of those patients were re-profiled on the
EPIC array for the 2025 work; the only additional 2025 samples came from a
collaborator in Spain. The 2025 derivation cohort is thus essentially a **superset
of the individuals scored here.** The circularity is at the **individual** level:
the episignature was derived on (at least most of) the same people it is scored on.
One nuance — because many individuals were re-profiled on a *different* platform
(EPIC) for 2025 while we score the deposited **450K** data, this is
same-individual / cross-platform overlap rather than reuse of identical arrays, so
the technical (batch/array) noise is partly independent; but the signal an
episignature exploits is the individuals' biology, which is shared. Consequently
van der Laan's strong performance here (AUC 0.93–0.96) is **not independent
validation and should not be reported as such** — it is an **upper bound**
reflecting a partly circular, same-individual evaluation. An exact per-sample
overlap count, and whether *any* scored individuals were held out of the 2025
derivation (training and CpG selection), remain to be confirmed (follow-up in
progress); absent a clean held-out subset, no independent AUC can be computed on
these cohorts. This does *not* affect the two findings that do not depend on it:
(a) the buccal signatures (Portales-Casamar, Lussier) are independent of these
cohorts and still transfer poorly to blood; and (b) between-signature agreement
(C2/C3, cross-tissue) is invariant to which signature performs best. It does
temper the "blood-native signature clearly wins" framing and inflates part of the
cross-tissue disagreement magnitude.

---

## References

- Breitling LP et al. (2011). Tobacco-smoking-related differential DNA
  methylation: 27K discovery and replication. Am J Hum Genet.
- Caramaschi D et al. (2022). Blood DNA methylation signatures of lifestyle
  exposures: tobacco and alcohol consumption. Clin Epigenetics.
- Joehanes R et al. (2016). Epigenetic signatures of cigarette smoking.
  Circ Cardiovasc Genet.
- Liu Y et al. (2013). Epigenome-wide association data implicate DNA methylation
  as an intermediary of genetic risk in rheumatoid arthritis. Nat Biotechnol.
- Lussier AA et al. (2018). DNA methylation as a predictor of fetal alcohol
  spectrum disorder. Clin Epigenetics 10:5.
- Portales-Casamar E et al. (2016). DNA methylation signature of human fetal
  alcohol spectrum disorder. Epigenetics & Chromatin 9:25.
- Sharp GC et al. (2018; PACE Consortium). Maternal alcohol consumption and
  offspring DNA methylation: findings from six general population birth cohorts.
- van der Laan L et al. (2025). Discovery of a DNA methylation episignature as a
  molecular biomarker for fetal alcohol syndrome. Genet Med 27(12):101586.
- Shenker NS et al. (2013). Epigenome-wide association study in the European
  Prospective Investigation into Cancer and Nutrition (EPIC-Turin) cohort.
  PLoS One.
- Vinogradova P (2025). Towards a formal definition of kinase inhibitor
  selectivity. ChemRxiv. doi:10.26434/chemrxiv.15001618/v1
