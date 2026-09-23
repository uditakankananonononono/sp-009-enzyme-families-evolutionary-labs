# SP-009 R1 — Shuffled-label control anomaly: investigation and resolution
Status: RESOLVED (mechanism identified, quantified, reproducible). Resolved 2026-09-22 ~06:55 IST, before R1 report ship.

## 1. Symptom
The frozen within-family shuffled-label control (protocol R1: labels permuted within each
training family, model refit, scored on the held-out family) returned null ROC-AUCs far
from 0.5: HXK4 median 0.715, all 30 reps >= 0.565; 100-draw experiments gave frac(AUC>0.5)
= 1.0 and mean cos(w, d_true) = +0.143 +/- 0.072 across three independent permutation
implementations (numpy Generator.permutation, argsort-of-random-keys, RandomState) and
multiple seeds, reproduced in a fresh process from saved arrays.

## 2. Hypotheses eliminated
- Shuffle mechanics broken: NO. corr(Y_perm, Y_train) ~ 0; complements verified exact
  (w(1-yp) = -w(yp) to machine precision; AUC(w)+AUC(-w) = 1.0 exactly).
- RNG/state corruption: NO. Three generators, many seeds, fresh process, saved arrays -
  identical values.
- Code-path bug in fit/score: NO. Synthetic single-family data gives null mean 0.4987.
- Scoring against the wrong label vector: NO. Scoring the same shuffled fits against a
  shuffled Yte gives mean AUC 0.5003 (clean).
- Floating-point/ill-conditioning artifact: NO. Bias is deterministic, large (20+ SE),
  and survives solver changes.
- Violation of the permutation symmetry theorem: NO. Exhaustive enumeration of all 924
  arrangements of a balanced 6-of-12 label multiset through the SAME fitted operator G
  gives mean cos = 0.0 exactly, pair sums 0.0. The theorem holds where its premise holds.

## 3. Mechanism (confirmed analytically and empirically)
The within-family permutation space is NOT closed under label complementation.
Training positive rates differ by family: BLAT p = 0.5004 (n = 5079), CBS p = 0.4589
(n = 7134), pooled 0.4761 (split point recovered at row 5079 by variance minimization).
Under within-family permutation, E[y_perm] is the per-family constant p_f, so the pooled
mean-centered label vector has a DETERMINISTIC nonzero component
  mu = (p_family - p_pooled) per row, ||mu|| = 12.60.
The OLS/dual operator G maps mu to a fixed weight direction w0 = G @ mu with
  cos(w0, d_true) = +0.3516,
where d_true is the real held-out-family signal direction. The deterministic component
ALONE scores ROC-AUC 0.7678 against the real HXK4 test labels - matching the observed
inflated null (HXK4 shuffled-null p95 = 0.7609).
Every within-family shuffled fit therefore equals w0 + permutation noise, centering the
null AUC distribution near 0.70-0.77 instead of 0.5. The permutation symmetry theorem is
not violated: the complement of a within-family arrangement has n_f - k_f ones in family
f, which the within-family sampler can never draw (k_f != n_f/2), so the exchangeability
premise fails and no zero-mean symmetry applies.
Root cause in scientific terms: between-family base-rate differences align with family
mean shifts on correlated feature axes (conservation/profile features, max off-diagonal
correlation 0.961), and those same axes predict the held-out family's labels. The
within-family null leaks family-shift structure into every null fit.

## 4. Corroborating controls
- GLOBAL permutation (pool rows across families, permute once): mean cos(w, d_true) =
  -0.0012, frac>0 = 0.55 over 60 draws - clean null. The within-family constraint is
  the sole cause.
- Synthetic single-family data: 0.4987 - clean (no between-family structure to leak).
- Shuffled-Yte scoring: 0.5003 - clean (w0 is orthogonal to a random label direction).

## 5. Implications for the frozen R1 gates (no gate text modified)
- The frozen control compares every model variant against the SAME within-family null
  procedure, so it remains a valid RELATIVE control, and the inflation makes it
  CONSERVATIVE (the bar is above 0.5, not below).
- The null values must not be read as "0.5-ish chance": they quantify family-shift
  leakage, which is itself a substantive R1 finding for the OOD estimand.
- Documented execution adaptation (analysis added, gates untouched): a global-permutation
  null diagnostic is reported alongside the frozen within-family null so readers see both
  the inflated family-structured null and the clean unstructured null.
- This investigation is a locked-gate integrity disclosure and ships with the R1 report.
