import json, re, math
import numpy as np
import pyarrow.parquet as pq
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef, balanced_accuracy_score, roc_curve, brier_score_loss
from scipy.stats import spearmanr
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
SINGLE=re.compile(r'^([A-Z])(\d+)([A-Z])$')
PRIMARY={'OXDA':'OXDA_RHOTO_Vanella_2023_expression','HXK4':'HXK4_HUMAN_Gersing_2022_activity',
 'BLAT':'BLAT_ECOLX_Deng_2012','CBS':'CBS_HUMAN_Sun_2020','PPM1D':'PPM1D_HUMAN_Miller_2022'}
TARGETS=json.load(open('results/target_seqs.json'))
TSEQ={'OXDA':TARGETS['OXDA_RHOTO_Vanella_2023_expression'],'HXK4':TARGETS['HXK4_HUMAN_Gersing_2022_activity'],
 'BLAT':TARGETS['BLAT_ECOLX_Deng_2012'],'CBS':TARGETS['CBS_HUMAN_Sun_2020'],'PPM1D':TARGETS['PPM1D_HUMAN_Miller_2022']}
EST=['HXK4','BLAT','CBS']; NONEST=['OXDA','PPM1D']
data={}
for i in range(5):
    t=pq.read_table(f'data/raw/dms_v1_shard{i}.parquet', columns=['DMS_id','mutant','DMS_score','DMS_score_bin'])
    d=t.to_pydict()
    for dms,mut,sc,b in zip(d['DMS_id'],d['mutant'],d['DMS_score'],d['DMS_score_bin']):
        for fam,assay in PRIMARY.items():
            if dms==assay:
                m=SINGLE.match(mut)
                if m: data.setdefault(fam,[]).append((m.group(1),int(m.group(2)),m.group(3),sc,b))
    del t,d
FEATS={}
for fam,rows in data.items():
    al=np.load(f'results/align_{fam}_features.npy'); tgt=TSEQ[fam]
    X=[];Y=[];Z=[];POS=[]
    for wt,pos,mt,sc,b in rows:
        assert tgt[pos-1]==wt, f'{fam} coordinate assertion failed'
        f=al[pos-1,2:]; ent=al[pos-1,0]; rel=al[pos-1,1]
        fs=(f+PC/20)/(1+PC); iwt,imt=AA.index(wt),AA.index(mt)
        llr=math.log(fs[imt]/fs[iwt])
        cons=1-ent/math.log2(20); bl=BL[wt][mt]
        X.append([cons,bl,llr,ent,rel,float(fs[imt]),float(fs[iwt])]); Y.append(int(b)); Z.append(sc); POS.append(pos)
    FEATS[fam]={'X':np.array(X),'Y':np.array(Y),'Z':np.array(Z),'POS':np.array(POS),
                'COV':np.load(f'results/align_{fam}_poscov.npy')}
GROUPS={'conservation_only':[0],'blosum_only':[1],'profile_only':[2],'flexibility_only':[3,4,5,6],'all':[0,1,2,3,4,5,6]}
def metrics(yt,pt,thr):
    pred=(pt>=thr).astype(int)
    A=np.vstack([pt,np.ones_like(pt)]).T
    sl,ic=np.linalg.lstsq(A,yt,rcond=None)[0]
    ece=0.0
    for k in range(10):
        m=(pt>=k/10)&(pt<(k+1)/10)
        if m.sum()>0: ece+=m.mean()*abs(yt[m].mean()-pt[m].mean())
    return {'roc_auc':float(roc_auc_score(yt,pt)),'pr_auc':float(average_precision_score(yt,pt)),
     'mcc':float(matthews_corrcoef(yt,pred)),'bal_acc':float(balanced_accuracy_score(yt,pred)),
     'brier':float(brier_score_loss(yt,pt)),'rel_slope':float(sl),'rel_intercept':float(ic),'ece':float(ece)}
def boot_ci_metric(yt,pt,pos,B,seed,block):
    rng=np.random.default_rng(seed); upos=np.unique(pos); n=len(yt); rocs=[];prs=[]
    for _ in range(B):
        if block:
            ps=rng.choice(upos,size=len(upos),replace=True)
            idx=np.concatenate([np.where(pos==p)[0] for p in ps])
        else: idx=rng.integers(0,n,n)
        if yt[idx].sum() in (0,len(idx)): continue
        rocs.append(roc_auc_score(yt[idx],pt[idx])); prs.append(average_precision_score(yt[idx],pt[idx]))
    return {'roc':[float(np.quantile(rocs,0.025)),float(np.quantile(rocs,0.975))],
            'pr':[float(np.quantile(prs,0.025)),float(np.quantile(prs,0.975))]}
def boot_delta(yt,pa,pb,pos,B,seed,block,metric):
    rng=np.random.default_rng(seed); upos=np.unique(pos); n=len(yt); d=[]
    fn=roc_auc_score if metric=='roc' else average_precision_score
    for _ in range(B):
        if block:
            ps=rng.choice(upos,size=len(upos),replace=True)
            idx=np.concatenate([np.where(pos==p)[0] for p in ps])
        else: idx=rng.integers(0,n,n)
        if yt[idx].sum() in (0,len(idx)): continue
        d.append(fn(yt[idx],pa[idx])-fn(yt[idx],pb[idx]))
    d=np.array(d)
    return [float(np.median(d)),float(np.quantile(d,0.025)),float(np.quantile(d,0.975))]
def fit_pred(train_fams,test_fam,cols,seed=0):
    Xtr=np.vstack([FEATS[f]['X'][:,cols] for f in train_fams]); Ytr=np.concatenate([FEATS[f]['Y'] for f in train_fams])
    sc_=StandardScaler().fit(Xtr)
    clf=LogisticRegression(max_iter=2000,class_weight='balanced',random_state=seed).fit(sc_.transform(Xtr),Ytr)
    ptr=clf.predict_proba(sc_.transform(Xtr))[:,1]
    fpr,tpr,ths=roc_curve(Ytr,ptr); thr=float(ths[np.argmax(tpr-fpr)])
    pte=clf.predict_proba(sc_.transform(FEATS[test_fam]['X'][:,cols]))[:,1]
    return pte,thr,Ytr
results={'estimable':{},'non_estimable':{}}
rng=np.random.default_rng(20260922)
for hold in EST:
    tr=[f for f in EST if f!=hold]
    Yte=FEATS[hold]['Y']; Pte=FEATS[hold]['POS']
    res={}; probs={}
    for g,cols in GROUPS.items():
        pte,thr,Ytr=fit_pred(tr,hold,cols)
        res[g]=metrics(Yte,pte,thr); probs[g]=pte
    prior=float(np.concatenate([FEATS[f]['Y'] for f in tr]).mean())
    pp=np.full(len(Yte),prior); res['class_prior']=metrics(Yte,pp,0.5); probs['class_prior']=pp
    sh=[]
    for r in range(20):
        pte,_t,_=fit_pred(tr,hold,GROUPS['all'],seed=100+r)
        # shuffle within training families (frozen: labels permuted within family)
        Xtr=np.vstack([FEATS[f]['X'][:,GROUPS['all']] for f in tr])
        Ys=np.concatenate([rng.permutation(FEATS[f]['Y']) for f in tr])
        sc_=StandardScaler().fit(Xtr)
        clf=LogisticRegression(max_iter=2000,class_weight='balanced').fit(sc_.transform(Xtr),Ys)
        ps=clf.predict_proba(sc_.transform(FEATS[hold]['X'][:,GROUPS['all']]))[:,1]
        sh.append(float(roc_auc_score(Yte,ps)))
    res['shuffled_control']={'roc_auc_median':float(np.median(sh)),'roc_auc_p95':float(np.quantile(sh,0.95)),'reps':20}
    Xtr=np.vstack([FEATS[f]['X'][:,GROUPS['all']] for f in tr])
    Ztr=np.concatenate([FEATS[f]['Z'] for f in tr])
    scx=StandardScaler().fit(Xtr)
    rg=Ridge().fit(scx.transform(Xtr),Ztr)
    res['ridge_all']={'spearman':float(spearmanr(FEATS[hold]['Z'],rg.predict(scx.transform(FEATS[hold]['X'][:,GROUPS['all']]))).statistic)}
    boot={}
    for base in ['conservation_only','profile_only']:
        for met in ['roc','pr']:
            boot[f'all_vs_{base}_{met}']={
              'variant_boot':boot_delta(Yte,probs['all'],probs[base],Pte,10000,20260922,False,met),
              'position_block_boot':boot_delta(Yte,probs['all'],probs[base],Pte,10000,20260922,True,met)}
    res['bootstraps']=boot
    results['estimable'][hold]=res
    print('EST',hold,{g:round(res[g]['roc_auc'],4) for g in list(GROUPS)+['class_prior']},'shuf95',round(res['shuffled_control']['roc_auc_p95'],4),flush=True)
for hold in NONEST:
    tr=EST
    Yte=FEATS[hold]['Y']; Pte=FEATS[hold]['POS']; cov=FEATS[hold]['COV']
    covered=np.array([cov[p-1]>=0.5 for p in Pte])
    res={'abstention_coverage':float(covered.mean()),'note':'baselines only; flexibility non-estimable per frozen source gate; no win/loss counting'}
    for g in ['conservation_only','blosum_only','profile_only']:
        cols=GROUPS[g]
        pte,thr,Ytr=fit_pred(tr,hold,cols)
        if g=='blosum_only':
            m=metrics(Yte,pte,thr)
            m['bootstrap_ci_variant']=boot_ci_metric(Yte,pte,Pte,10000,20260922,False)
            m['bootstrap_ci_position_block']=boot_ci_metric(Yte,pte,Pte,10000,20260922,True)
            res[g]=m
        else:
            if covered.sum()<100: res[g]={'abstained':'<100 covered variants','covered_n':int(covered.sum())}; continue
            m=metrics(Yte[covered],pte[covered],thr)
            m['covered_n']=int(covered.sum())
            m['bootstrap_ci_variant']=boot_ci_metric(Yte[covered],pte[covered],Pte[covered],10000,20260922,False)
            m['bootstrap_ci_position_block']=boot_ci_metric(Yte[covered],pte[covered],Pte[covered],10000,20260922,True)
            res[g]=m
    prior=float(np.concatenate([FEATS[f]['Y'] for f in tr]).mean())
    res['class_prior']=metrics(Yte,np.full(len(Yte),prior),0.5)
    results['non_estimable'][hold]=res
    print('NONEST',hold,'cov',round(float(covered.mean()),3),{g:(round(res[g]['roc_auc'],4) if isinstance(res.get(g),dict) and 'roc_auc' in res[g] else res.get(g)) for g in ['conservation_only','blosum_only','profile_only','class_prior']},flush=True)
json.dump(results,open('results/r1_model_results.json','w'),indent=1)
print('DONE')
