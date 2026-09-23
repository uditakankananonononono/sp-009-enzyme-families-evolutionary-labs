import json, subprocess, hashlib, random, re, csv, io, os
def sha_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()
targets=json.load(open('results/target_seqs.json'))
tlen={'OXDA':364,'HXK4':465,'BLAT':286,'CBS':551,'PPM1D':605,'AMIE':346}
tseq={'OXDA':targets['OXDA_RHOTO_Vanella_2023_expression'],'HXK4':targets['HXK4_HUMAN_Gersing_2022_activity'],
      'BLAT':targets['BLAT_ECOLX_Deng_2012'],'CBS':targets['CBS_HUMAN_Sun_2020'],
      'PPM1D':targets['PPM1D_HUMAN_Miller_2022'],'AMIE':targets['AMIE_PSEAE_Wrenbeck_2017']}
queries={'OXDA':'xref:pfam-PF01266','HXK4':'(xref:pfam-PF00349) AND (xref:pfam-PF03727)',
 'BLAT':'xref:pfam-PF13354','CBS':'(xref:pfam-PF00571) AND (xref:pfam-PF00291)',
 'PPM1D':'xref:pfam-PF00481','AMIE':'xref:pfam-PF00795'}
BAD=set('XBZJUO')
rng=random.Random(20260922)
manifest={}
os.makedirs('data/homologs',exist_ok=True)
for fam,q in queries.items():
    L=tlen[fam]; lo,hi=int(L*0.5),int(L*2)
    full=f'({q}) AND NOT fragment:true AND length:[{lo} TO {hi}]'
    url=("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query="
         + subprocess.run(['python3','-c',f"import urllib.parse;print(urllib.parse.quote('''{full}'''))"],capture_output=True,text=True).stdout.strip()
         +"&fields=accession,length,organism_name,lineage,sequence")
    raw=f'data/homologs/{fam}_raw.tsv.gz'
    subprocess.run(['curl','-s','--max-time','600','-L','-o',raw,url],check=True)
    import gzip
    try:
        txt=gzip.open(raw,'rt').read()
    except Exception:
        txt=open(raw).read()  # server may send plain
    rows=list(csv.DictReader(io.StringIO(txt),delimiter='\t'))
    n0=len(rows)
    # quality filters
    tgt=tseq[fam]
    def best_window_ident(s):
        short,lon=(s,tgt) if len(s)<=len(tgt) else (tgt,s)
        best=0
        for off in range(0,len(lon)-len(short)+1):
            w=lon[off:off+len(short)]
            ident=sum(1 for x,y in zip(short,w) if x==y)/len(short)
            if ident>best: best=ident
            if best>=0.97: break
        return best
    def km3(x): return set(x[i:i+3] for i in range(len(x)-2))
    Kt=km3(tgt)
    kept=[]
    for r in rows:
        s=r['Sequence']
        if any(c in BAD for c in s): continue
        if s==tgt: continue
        K=km3(s); j=len(K&Kt)/len(K|Kt)
        if j>=0.85 and best_window_ident(s)>=0.97: continue
        lin=r.get('Taxonomic lineage') or ''
        m=re.search(r', (\w+) \(domain\)',lin)
        sk=m.group(1) if m else 'Unknown'
        kept.append({'acc':r['Entry'],'seq':s,'sk':sk})
        pass
    n1=len(kept)
    # taxonomic cap 60% of 1500 = 900 per superkingdom; sample seeded
    by={}
    for k in kept: by.setdefault(k['sk'],[]).append(k)
    cap=int(1500*0.6)
    pool=[]
    for sk,lst in by.items():
        rng.shuffle(lst)
        pool+=lst[:cap]
    rng.shuffle(pool)
    pool=pool[:1500]
    # greedy kmer3-jaccard 0.8 redundancy removal
    def km(s): return set(s[i:i+3] for i in range(len(s)-2))
    surv=[]; kmers=[]
    for k in pool:
        K=km(k['seq'])
        red=False
        for K2 in kmers:
            inter=len(K&K2)
            if inter and inter/len(K|K2)>=0.8: red=True; break
        if not red:
            surv.append(k); kmers.append(K)
    fa=f'data/homologs/{fam}_submit.fasta'
    with open(fa,'w') as f:
        f.write(f'>TARGET_{fam}\n{tgt}\n')
        for k in surv: f.write(f">{k['acc']}|{k['sk']}\n{k['seq']}\n")
    manifest[fam]={'query':full,'raw_rows':n0,'post_quality':n1,'post_taxcap':len(pool),'post_cluster':len(surv),
        'sk_counts':{sk:sum(1 for k in surv if k['sk']==sk) for sk in by},
        'raw_sha256':sha_file(raw),'submit_sha256':sha_file(fa)}
    print(fam, manifest[fam])
json.dump(manifest,open('results/r1_homolog_fetch_manifest.json','w'),indent=1)
