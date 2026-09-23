import subprocess, json, os, time
EMAIL='agent@instinct.co'
jobs={}
for fam in ['OXDA','HXK4','BLAT','CBS','PPM1D','AMIE']:
    fa=f'data/homologs/{fam}_submit.fasta'
    r=subprocess.run(['curl','-s','--max-time','120','-X','POST',
        'https://www.ebi.ac.uk/Tools/services/rest/mafft/run',
        '-F',f'email={EMAIL}','-F','stype=protein','-F','format=fasta',
        '-F',f'sequence=@{fa};type=text/plain'],capture_output=True,text=True)
    jid=r.stdout.strip()
    jobs[fam]=jid
    print(fam, jid, r.stderr[:100])
json.dump(jobs,open('results/r1_mafft_jobs.json','w'),indent=1)
