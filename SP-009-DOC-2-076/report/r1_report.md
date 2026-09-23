# SP-009 R1 - Enzyme Families as Evolutionary Laboratories (DOC-2-076)
## Round 1 final report: evolutionary-flexibility add-on under distribution shift

**Verdict: FAIL - the flexibility add-on cannot graduate under the frozen five-family
estimand (3/5 estimable). Recommendation: CLOSE the method (R2 not triggered).**

- Protocol R1 lock sha256: `427a2ddd0441b03124e5e99ab72364d47c7e18c45776a17dad1649f3662c0d73` (locked 2026-09-22 03:41 IST, before any outcome inspection)
- Frozen alignment-source gate record sha256: `2ebffa75c579e55b69af6c005b690a10e1f84a23b7ef06d422a804253b4d00d3`
- Gate evaluation sha256: `837b2abb73fa6ac15ff14b44dd22711c872ee40838742bc1a0ff739b6e542937`
- Orchestrator route (2026-09-22 05:35 IST): proceed exactly as locked, no pre-outcome amendment; baselines only for non-estimable families; no imputation, no relaxed alignment, no AMIE substitute, no win/loss counting.

## 1. Estimand and frozen decision rule
Leave-one-family-out (LOFO) over five frozen families: OXDA, HXK4, BLAT, CBS, PPM1D
(AMIE backup-only). Frozen rule (verbatim): "all-features must beat conservation-only
AND profile-only with 95% CI excluding zero under BOTH bootstraps, in >=4/5 families,
for BOTH ROC-AUC and PR-AUC." Both uncertainty analyses: variant bootstrap and
position-block (hierarchical site) bootstrap. Heterogeneity clause: one severe reversal
(CI excluding zero negative) blocks any aggregate pass claim regardless of count.
Labels: ProteinGym frozen cutoffs, bin==1 iff score>=cutoff (0 exceptions, verified R0).
Preprocessing fit inside training folds only.

## 2. Alignment-source gate (frozen before modeling)
| family | Neff | coverage | gate |
|---|---|---|---|
| HXK4 | 412.2 | 0.940 | PASS |
| BLAT | 424.8 | 0.871 | PASS |
| CBS | 234.9 | 0.927 | PASS |
| OXDA | - | 0.736 | FAIL (<0.80) |
| PPM1D | 32.5 | 0.688 | FAIL |
| AMIE (backup) | 164.0 | 0.708 | FAIL |

3/5 estimable; the single permitted AMIE backup swap was consumed and still failed, so
graduation under the frozen five-family estimand was UNREACHABLE before modeling began.
OXDA/PPM1D were scored on baselines only with abstention reporting, per the frozen
failure rule and orchestrator route (a).

### Documented execution adaptations (gate text untouched)
1. EBI MAFFT REST 500-sequence cap: seeded stratified subsample to 499 homologs + target
   per family (recorded in r1_homolog_fetch_manifest.json with sha256s).
2. MAFFT jobs queued >90 min: Clustal Omega REST fallback used for all 6 families
   (r1_clustalo_jobs.json; alignment sha256s in manifest).
3. Global-permutation null diagnostic added alongside the frozen within-family
   shuffled control (section 5; analysis addition only, no gate change).

## 3. Per-family results (estimable families, LOFO)
ROC-AUC (held-out family):

| model | HXK4 | BLAT | CBS |
|---|---|---|---|
| conservation_only | 0.7654 | 0.6540 | 0.6949 |
| blosum_only (control) | 0.5809 | 0.6524 | 0.5815 |
| profile_only | 0.7497 | 0.6903 | 0.6921 |
| flexibility_only | 0.7728 | 0.6825 | 0.7030 |
| **all-features** | **0.7644** | **0.7149** | **0.7137** |
| class_prior (control) | 0.5000 | 0.5000 | 0.5000 |

Full battery (PR-AUC, MCC, balanced accuracy at training-only threshold, Brier,
reliability slope/intercept, ECE) is in results/r1_model_results.json - no success was
claimed on ROC-AUC alone; where ROC and PR disagree we report both (see HXK4 reversal).

### Frozen decision-rule evaluation (all-features vs each baseline, 95% CI under both bootstraps)
| comparison | HXK4 | BLAT | CBS |
|---|---|---|---|
| vs conservation ROC | -0.001 vb[-.009,+.007] pb[-.011,+.009] | +0.061 vb[+.052,+.070] pb[+.046,+.077] WIN | +0.019 vb[+.013,+.025] pb[+.011,+.027] WIN |
| vs conservation PR | **-0.018 vb[-.031,-.004] pb[-.034,-.001] REVERSAL** | +0.088 vb[+.074,+.102] pb[+.059,+.110] WIN | +0.021 vb[+.011,+.031] pb[+.010,+.034] WIN |
| vs profile ROC | +0.015 vb[+.008,+.021] pb[+.002,+.028] WIN | +0.025 vb[+.015,+.034] pb[-.001,+.050] no (pb) | +0.022 vb[+.016,+.027] pb[+.014,+.030] WIN |
| vs profile PR | +0.009 vb[-.001,+.020] pb[-.006,+.025] no | +0.045 vb[+.031,+.058] pb[+.012,+.078] WIN | +0.026 vb[+.017,+.035] pb[+.015,+.038] WIN |
| **family passes rule** | NO | NO | **YES** |

**1/3 estimable families pass; 1/5 under the frozen estimand; 4/5 required. FAIL.**
Heterogeneity clause TRIGGERED: HXK4 all-vs-conservation PR-AUC is a significant
negative reversal under both bootstraps, blocking any aggregate pass claim.

### Non-estimable families (baselines only, no win/loss counting)
| family | abstention coverage | conservation | blosum | profile | class_prior |
|---|---|---|---|---|---|
| OXDA | 0.737 | 0.6606 | 0.6098 | 0.6647 | 0.5 |
| PPM1D | 0.806 | 0.7227 | 0.5611 | 0.7157 | 0.5 |

### Continuous-score ridge transfer (secondary, descriptive)
Spearman on held-out family: HXK4 -0.387 (sign-inverted transport), BLAT +0.404,
CBS +0.052. Heterogeneous; the HXK4 sign inversion is a preserved negative.

## 4. Shuffled-label control: anomaly, investigation, resolution
The frozen within-family shuffled-label control (labels permuted within family,
seed 20260922, 20 reps) returned inflated nulls: HXK4 median 0.721 / p95 0.761;
BLAT median 0.522 / p95 0.696; CBS median 0.478 / p95 0.638. An un-explained
0.7-median null would have blocked shipping; it was investigated to root cause
(results/r1_shuffled_null_investigation.md, sha256
4d6a2c8a0e593cec6892ce4aec429b622facd205758386f836bd4045f8962af7).

**Mechanism (proven, not a bug):** training positive rates differ by family
(BLAT 0.5004, CBS 0.4589), so under within-family permutation E[y_perm] is the
per-family constant; pooled mean-centering leaves a deterministic offset mu that the
linear operator maps to a fixed direction w0 with cos(w0, true test-signal direction)
= +0.35; w0 alone scores AUC 0.7678 on real HXK4 labels. Between-family base-rate
differences leak through family-correlated feature axes (max off-diagonal feature
correlation 0.961) that also predict the held-out family. The permutation symmetry
theorem is not violated: within-family arrangements are not closed under
complementation (complement of a k_f-of-n_f arrangement has n_f-k_f ones, unreachable
by the sampler). Exhaustive enumeration of all 924 arrangements of a balanced 6-of-12
multiset through the same fitted operator gives mean cos = 0.0 exactly; global
permutation, synthetic data, and shuffled-Yte scoring all give clean nulls.

**Consequences:** the frozen control remains a valid RELATIVE comparison (every model
variant faces the identical inflated null - it is conservative, and all three
all-features models beat their family p95: 0.7644>0.761, 0.7149>0.696, 0.7137>0.638,
HXK4 only marginally). The inflation is itself a substantive finding: family-shift
base-rate leakage is large enough to reach AUC ~0.77 with zero label information,
which recalibrates how much of the raw LOFO performance is family-structure transport
rather than variant biology.

**Global-permutation diagnostic (added analysis):** HXK4 median 0.598 / p95 0.766;
BLAT median 0.543 / p95 0.651; CBS median 0.482 / p95 0.603 (20 reps, seed 20260922).
CBS nulls center at ~0.5 as expected; HXK4 retains elevated spread, consistent with
strongly family-correlated feature structure rather than any label leakage (features
never see test labels; preprocessing is fit inside training folds).

## 5. Interpretation
- The flexibility add-on's locked incremental claim fails: one clean win (CBS,
  +0.02 ROC / +0.03 PR over both baselines under both bootstraps), one partial win
  (BLAT: large gain over conservation, not separable from profile under the
  position-block bootstrap), one flat-to-negative (HXK4: no ROC gain, significant
  PR-AUC reversal vs conservation).
- flexibility_only alone is the best single model on HXK4 (0.7728) - descriptive,
  not a locked comparison; it does not rescue the add-on because all-features is the
  locked unit and it does not beat the conservation baseline there.
- Position-block bootstrap CIs are consistently wider than variant bootstrap CIs,
  as expected when variants at the same site are correlated; the BLAT profile
  comparison is decided by this difference, which is exactly why the amendment
  required both.
- The shuffled-null investigation shows raw LOFO AUCs up to ~0.77 are attainable
  from family-structure transport alone; absolute LOFO numbers on these families
  should not be read as variant-effect skill without the null decomposition.

## 6. R2 decision (orchestrator criterion, verbatim route (a))
"R2 is not automatic: only if consistent incremental signal under both bootstraps
does a new round preregister a narrower estimand on fresh families; otherwise close
the method." The signal is NOT consistent (1/3 pass, one significant reversal).
**Recommendation: CLOSE the evolutionary-flexibility add-on.** Per protocol
r2_preview, retain the best simple baseline as the useful negative: conservation-only
(HXK4 0.7654, CBS 0.6949) / profile-only (BLAT 0.6903), with the source-gate floors
(Neff, coverage) and the family-structure null decomposition as the durable method
contributions of this round.

## 7. Integrity disclosures
- Graduation was unreachable from the moment the source gate froze (3/5 estimable);
  all modeling was run and reported anyway per orchestrator route (a); no
  imputation, relaxed alignment, AMIE substitute, or win/loss counting occurred.
- The shuffled-null anomaly was detected, investigated to mechanism, and resolved
  BEFORE this report shipped; the full diagnostic trail is preserved.
- One compute ceiling preserved: torch-class frameworks exceed the 1GB sandbox
  (environment_ledger.json); 7-feature logistic models and dual-form ridge were used.
- Live sources only: ProteinGym v1 (sha256-verified shards), UniProt REST, EBI job
  services; Ensembl unreachable (never used); no simulated data, no fabricated
  citations.
