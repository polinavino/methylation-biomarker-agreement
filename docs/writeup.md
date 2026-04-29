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
ground truth is unavailable, with planned application to fetal alcohol spectrum
disorder methylation signatures.

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
Only CpG sites present in all four signatures were retained (n=12 CpGs).

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

**Future directions:** We plan to apply this framework to fetal alcohol spectrum
disorder (FASD) methylation signatures, where multiple proposed biomarker
signatures exist but cross-cohort validation is limited by restricted data
access. The smoking methylation demonstration establishes the technical
methodology and provides a benchmark for interpreting FASD signature agreement.

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
- Shenker NS et al. (2013). Epigenome-wide association study in the European
  Prospective Investigation into Cancer and Nutrition (EPIC-Turin) cohort.
  PLoS One.
- Vinogradova P (2025). Towards a formal definition of kinase inhibitor
  selectivity. ChemRxiv. doi:10.26434/chemrxiv.15001618/v1
