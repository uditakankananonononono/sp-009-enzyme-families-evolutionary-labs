import json, hashlib, math, os
import numpy as np
AA='ARNDCQEGHILKMFPSTWYV'
BG=dict(zip('ARNDCQEGHILKMFPSTWYV',[.078,.051,.045,.054,.019,.038,.053,.072,.022,.053,.091,.059,.023,.039,.042,.068,.059,.014,.032,.066]))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()
def read_fa(p):
    recs=[]; name=None; buf=[]
    for l in open(p):
        l=l.strip()
        if l.startswith('>'):
            if name is not None: recs.append((name,''.join(buf)))
            name=l[1:]; buf=[]
        elif l: buf.append(l)
    if name is not None: recs.append((name,''.join(buf)))
    return recs
targets=json.load(open('results/target_seqs.json'))
tseq={'OXDA':targets['OXDA_RHOTO_Vanella_2023_expression'],'HXK4':targets['HXK4_HUMAN_Gersing_2022_activity'],
      'BLAT':targets['BLAT_ECOLX_Deng_2012'],'CBS':targets['CBS_HUMAN_Sun_2020'],
      'PPM1D':targets['PPM1D_HUMAN_Miller_2022'],'AMIE':targets['AMIE_PSEAE_Wrenbeck_2017']}
out={}
os.makedirs('data/alignments',exist_ok=True)
for fam,tgt in tseq.items():
    p=f'data/alignments/{fam}.aln.fasta'
    if not os.path.exists(p): out[fam]={'status':'alignment missing'}; continue
    recs=read_fa(p)
    tname=[n for n,_ in recs if n.startswith('TARGET')][0]
    trow=dict(recs)[tname]
    keep_cols=[i for i,c in enumerate(trow) if c!='-']
    assert ''.join(trow[i] for i in keep_cols)==tgt, f'{fam}: target-spanning residues do not reconstruct target'
    homs=[(n,s) for n,s in recs if not n.startswith('TARGET')]
    # drop >50% gaps over target-spanning columns
    surv=[]
    for n,s in homs:
        gaps=sum(1 for i in keep_cols if s[i]=='-')/len(keep_cols)
        if gaps<=0.5: surv.append((n,s))
    L=len(tgt)
    # Henikoff weights on target-spanning columns
    n=len(surv)
    w=np.ones(n)
    cols=[[s[i] for _,s in surv] for i in keep_cols]
    for col in cols:
        from collections import Counter
        cnt=Counter(c for c in col if c!='-')
        for j,c in enumerate(col):
            if c!='-' and cnt[c]>0: w[j]+=0 # placeholder
    # proper Henikoff: w_j = sum over positions 1/(n_pos * count(aa_j at pos))
    w=np.zeros(n)
    for col in cols:
        from collections import Counter
        ng=sum(1 for c in col if c!='-')
        cnt=Counter(c for c in col if c!='-')
        for j,c in enumerate(col):
            if c!='-' and ng>0: w[j]+=1.0/(ng*cnt[c])
    w=w/w.sum()
    neff=float(1.0/np.sum(w**2))
    # per-position weighted occupancy + column stats
    pos_cov=np.zeros(L); feats=np.zeros((L,22))  # entropy, rel_ent, 20-spectrum
    for pi,i in enumerate(keep_cols):
        col=cols[pi]
        wg=np.array([w[j] for j in range(n) if col[j]!='-'])
        if wg.sum()>0: pos_cov[pi]=wg.sum()
        f=np.zeros(20)
        for j in range(n):
            c=col[j]
            if c in AA: f[AA.index(c)]+=w[j]
        if f.sum()>0: f=f/f.sum()
        ent=float(-(f[f>0]*np.log2(f[f>0])).sum())
        bg=np.array([BG[a] for a in AA])
        rel=float((f[f>0]*np.log2(f[f>0]/bg[f>0])).sum())
        feats[pi,0]=ent; feats[pi,1]=rel; feats[pi,2:]=f
    cov80=float((pos_cov>=0.5).mean())
    np.save(f'results/align_{fam}_features.npy', feats)
    np.save(f'results/align_{fam}_poscov.npy', pos_cov)
    np.save(f'results/align_{fam}_weights.npy', w)
    out[fam]={'n_input_homologs':len(homs),'n_surviving':n,'Neff':round(neff,1),
        'target_spanning_cols':L,'position_coverage_ge50':round(cov80,4),
        'gate_pass':bool(cov80>=0.80 and neff>=200),
        'alignment_sha256':sha(p)}
json.dump(out,open('results/r1_alignment_gate.json','w'),indent=1)
print(json.dumps(out,indent=1))
