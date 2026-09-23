# Documented execution adaptation (gates untouched): global-permutation null diagnostic
# companion to the frozen within-family shuffled-label control. Same seed/reps.
# Mechanism reference: results/r1_shuffled_null_investigation.md
import json
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import importlib.util, sys, types, re

# Rebuild FEATS exactly as r1_models2.py does (same code path), without running modeling:
src = open('code/r1_models2.py').read()
cut = src.index("GROUPS=")  # stop right after FEATS construction
g = {'__name__':'__main__'}
exec(compile(src[:cut],'r1_models2.py','exec'), g)
FEATS = g['FEATS']; EST = g['EST']
GROUPS_ALL = [0,1,2,3,4,5,6]
rng = np.random.default_rng(20260922)
out = {}
for hold in EST:
    tr = [f for f in EST if f != hold]
    Xtr = np.vstack([FEATS[f]['X'][:,GROUPS_ALL] for f in tr])
    Ytr = np.concatenate([FEATS[f]['Y'] for f in tr])
    Yte = FEATS[hold]['Y']
    Xte = FEATS[hold]['X'][:,GROUPS_ALL]
    sh = []
    for rep in range(20):
        Ys = rng.permutation(Ytr)  # GLOBAL permutation across pooled training rows
        sc_ = StandardScaler().fit(Xtr)
        clf = LogisticRegression(max_iter=2000, class_weight='balanced').fit(sc_.transform(Xtr), Ys)
        ps = clf.predict_proba(sc_.transform(Xte))[:,1]
        sh.append(float(roc_auc_score(Yte, ps)))
    out[hold] = {'roc_auc_median': float(np.median(sh)), 'roc_auc_p95': float(np.quantile(sh,0.95)),
                 'roc_auc_min': float(np.min(sh)), 'roc_auc_max': float(np.max(sh)), 'reps': 20,
                 'permutation_scope': 'global_across_pooled_training_rows',
                 'seed': 20260922,
                 'rationale': 'diagnostic companion to frozen within-family null; see results/r1_shuffled_null_investigation.md'}
    print('GLOBALNULL', hold, out[hold], flush=True)
json.dump(out, open('results/r1_global_null_diagnostic.json','w'), indent=1)
print('WROTE results/r1_global_null_diagnostic.json', flush=True)
