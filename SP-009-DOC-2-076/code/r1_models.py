import json, re, math
import numpy as np
import pyarrow.parquet as pq
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef, balanced_accuracy_score, roc_curve, brier_score_loss
# BLOSUM62 (order ARNDCQEGHILKMFPSTWYV), from NCBI standard matrix
B62="""
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -3 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
AA='ARNDCQEGHILKMFPSTWYV'
BL={}
for line in B62.strip().splitlines():
    p=line.split(); BL[p[0]]={AA[i]:int(p[1+i]) for i in range(20)}
PC=0.5
def profile_feats(f, wt, mt):
    # f: 20-dim weighted column frequencies; add-0.5 pseudocount on normalized freqs (frozen)
    fs=(f+PC/20)/(1+PC)
    i_wt, i_mt = AA.index(wt), AA.index(mt)
    return fs[i_mt], fs[i_wt], math.log(fs[i_mt]/fs[i_wt])
SINGLE=re.compile(r'^([A-Z])(\d+)([A-Z])$')
PRIMARY={'OXDA':'OXDA_RHOTO_Vanella_2023_expression','HXK4':'HXK4_HUMAN_Gersing_2022_activity',
 'BLAT':'BLAT_ECOLX_Deng_2012','CBS':'CBS_HUMAN_Sun_2020','PPM1D':'PPM1D_HUMAN_Miller_2022'}
TARGETS=json=__import__('json').load(open('results/target_seqs.json'))
TSEQ={'OXDA':TARGETS['OXDA_RHOTO_Vanella_2023_expression'],'HXK4':TARGETS['HXK4_HUMAN_Gersing_2022_activity'],
 'BLAT':TARGETS['BLAT_ECOLX_Deng_2012'],'CBS':TARGETS['CBS_HUMAN_Sun_2020'],'PPM1D':TARGETS['PPM1D_HUMAN_Miller_2022']}
FAMS=list(PRIMARY)
# load variants
data={}
for i in range(5):
    t=pq.read_table(f'data/raw/dms_v1_shard{i}.parquet', columns=['DMS_id','mutant','DMS_score','DMS_score_bin'])
    d=t.to_pydict()
    for dms,mut,sc,b in zip(d['DMS_id'],d['DMS_score'],d['mutant'],d['DMS_score_bin']):
        for fam,assay in PRIMARY.items():
            if dms==assay:
                m=SINGLE.match(mut)
                if m: data.setdefault(fam,[]).append((m.group(1),int(m.group(2)),m.group(3),sc,b))
    del t,d
# feature matrix per family
FEATS={}
for fam,rows in data.items():
    al=np.load(f'results/align_{fam}_features.npy'); tgt=TSEQ[fam]
    X=[];Y=[];Z=[];POS=[]
    for wt,pos,mt,sc,b in rows:
        assert tgt[pos-1]==wt, f'{fam} coordinate assertion failed {wt}{pos}{mt}'
        f=al[pos-1,2:]; ent=al[pos-1,0]; rel=al[pos-1,1]
        f_m,f_w,llr=profile_feats(f,wt,mt)
        cons=1-ent/math.log2(20); bl=BL[wt][mt]
        X.append([cons,bl,llr,ent,rel,f_m,f_w]); Y.append(int(b)); Z.append(sc); POS.append(pos)
    FEATS[fam]={'X':np.array(X),'Y':np.array(Y),'Z':np.array(Z),'POS':np.array(POS)}
json.dump({k:len(v['Y']) for k,v in FEATS.items()},open('results/r1_variant_counts_loaded.json','w'))
# feature groups (indices): cons 0, blosum 1, profile 2, flexibility 3-6
GROUPS={'conservation_only':[0],'blosum_only':[1],'profile_only':[2],'flexibility_only':[3,4,5,6],'all':[0,1,2,3,4,5,6]}
def eval_metrics(yt,pt,pt_bin,thr):
    out={}
    out['roc_auc']=float(roc_auc_score(yt,pt)); out['pr_auc']=float(average_precision_score(yt,pt))
    pred=(pt>=thr).astype(int)
    out['mcc']=float(matthews_corrcoef(yt,pred)); out['bal_acc']=float(balanced_accuracy_score(yt,pred))
    out['brier']=float(brier_score_loss(yt,pt))
    # reliability slope/intercept
    A=np.vstack([pt,np.ones_like(pt)]).T
    sl,ic=np.linalg.lstsq(A,yt,rcond=None)[0]
    out['rel_slope']=float(sl); out['rel_intercept']=float(ic)
    # ECE 10-bin
    ece=0.0
    for k in range(10):
        m=(pt>=k/10)&(pt<(k+1)/10)
        if m.sum()>0: ece+=m.mean()*abs(yt[m].mean()-pt[m].mean())
    out['ece']=float(ece)
    return out
def delta_auc_boot(yt,p_all,p_base,pos,B=10000,seed=20260922,block=False):
    rng=np.random.default_rng(seed)
    upos=np.unique(pos); n=len(yt); d=[]
    for _ in range(B):
        if block:
            ps=rng.choice(upos,size=len(upos),replace=True)
            idx=np.concatenate([np.where(pos==p)[0] for p in ps])
        else:
            idx=rng.integers(0,n,n)
        if yt[idx].sum() in (0,len(idx)): continue
        d.append(roc_auc_score(yt[idx],p_all[idx])-roc_auc_score(yt[idx],p_base[idx]))
    d=np.array(d)
    return float(np.median(d)), float(np.quantile(d,0.025)), float(np.quantile(d,0.975))
results={}
rng_global=np.random.default_rng(20260922)
for hold in FAMS:
    tr=[f for f in FAMS if f!=hold]
    Xtr=np.vstack([FEATS[f]['X'] for f in tr]); Ytr=np.concatenate([FEATS[f]['Y'] for f in tr])
    Ztr=np.concatenate([FEATS[f]['Z'] for f in tr])
    Xte=FEATS[hold]['X']; Yte=FEATS[hold]['Y']; PosTe=FEATS[hold]['POS']
    # assertions: no held-out family data in training
    assert all(f!=hold for f in tr)
    res={}
    probs={}
    for g,cols in GROUPS.items():
        sc_=StandardScaler().fit(Xtr[:,cols])
        clf=LogisticRegression(max_iter=2000,class_weight='balanced').fit(sc_.transform(Xtr[:,cols]),Ytr)
        ptr=clf.predict_proba(sc_.transform(Xtr[:,cols]))[:,1]
        pte=clf.predict_proba(sc_.transform(Xte[:,cols]))[:,1]
        fpr,tpr,ths=roc_curve(Ytr,ptr); thr=float(ths[np.argmax(tpr-fpr)])
        res[g]=eval_metrics(Yte,pte,None,thr); probs[g]=pte
    # class prior control
    prior=Ytr.mean(); pp=np.full_like(Yte,prior,dtype=float)
    res['class_prior']=eval_metrics(Yte,pp,None,0.5); probs['class_prior']=pp
    # shuffled-label control (20 reps, all-features)
    sh=[]
    cols=GROUPS['all']
    for r in range(20):
        Ys=rng_global.permutation(Ytr)
        sc_=StandardScaler().fit(Xtr[:,cols])
        clf=LogisticRegression(max_iter=2000,class_weight='balanced').fit(sc_.transform(Xtr[:,cols]),Ys)
        pte=clf.predict_proba(sc_.transform(Xte[:,cols]))[:,1]
        sh.append(float(roc_auc_score(Yte,pte)))
    res['shuffled_control']={'roc_auc_median':float(np.median(sh)),'roc_auc_p95':float(np.quantile(sh,0.95)),'reps':20}
    # ridge spearman (all features)
    scz=StandardScaler().fit(Ztr.reshape(-1,1))
    # within-family z not needed for spearman; use raw scores per family
    scx=StandardScaler().fit(Xtr[:,cols])
    rg=Ridge().fit(scx.transform(Xtr[:,cols]),Ztr)
    from scipy.stats import spearmanr
    res['ridge_all']={'spearman':float(spearmanr(FEATS[hold]['Z'],rg.predict(scx.transform(Xte[:,cols]))).statistic)}
    # bootstraps: all vs conservation, all vs profile (ROC and PR deltas)
    boot={}
    for base in ['conservation_only','profile_only']:
        for met,pair in [('roc',None),('pr',None)]:
            key=f'all_vs_{base}_{met}'
            md_v,lo_v,hi_v=delta_auc_boot(Yte,probs['all'],probs[base],PosTe,block=False) if met=='roc' else (None,None,None)
            if met=='pr':
                rng=np.random.default_rng(20260922); upos=np.unique(PosTe); d=[]
                n=len(Yte)
                for _ in range(10000):
                    idx=rng.integers(0,n,n)
                    if Yte[idx].sum() in (0,len(idx)): continue
                    d.append(average_precision_score(Yte[idx],probs['all'][idx])-average_precision_score(Yte[idx],probs[base][idx]))
                d=np.array(d); md_v,lo_v,hi_v=float(np.median(d)),float(np.quantile(d,0.025)),float(np.quantile(d,0.975))
            md_b,lo_b,hi_b=delta_auc_boot(Yte,probs['all'],probs[base],PosTe,block=True) if met=='roc' else (None,None,None)
            if met=='pr':
                rng=np.random.default_rng(20260922); upos=np.unique(PosTe); d=[]
                for _ in range(10000):
                    ps=rng.choice(upos,size=len(upos),replace=True)
                    idx=np.concatenate([np.where(PosTe==p)[0] for p in ps])
                    if Yte[idx].sum() in (0,len(idx)): continue
                    d.append(average_precision_score(Yte[idx],probs['all'][idx])-average_precision_score(Yte[idx],probs[base][idx]))
                d=np.array(d); md_b,lo_b,hi_b=float(np.median(d)),float(np.quantile(d,0.025)),float(np.quantile(d,0.975))
            boot[key]={'variant_boot':[md_v,lo_v,hi_v],'position_block_boot':[md_b,lo_b,hi_b]}
    res['bootstraps']=boot
    results[hold]=res
    print(hold, json.dumps({g:round(res[g]['roc_auc'],4) for g in list(GROUPS)+['class_prior']}), 'shuf_p95', round(res['shuffled_control']['roc_auc_p95'],4))
json.dump(results,open('results/r1_model_results.json','w'),indent=1)
print('DONE')
