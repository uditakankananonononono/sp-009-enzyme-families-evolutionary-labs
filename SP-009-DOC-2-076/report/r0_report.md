# SP-009 / DOC-2-076 Enzyme Families as Evolutionary Laboratories - R0 Report

## Verdict: R0 GO (6/6 gates pass). Estimand identifiable.

## What was frozen (protocol sha256 49292b006e75c1066d9df96fcc5ab9a45cfc54b714dbe3be75f4c24fbeb45918, locked 03:39 IST before any outcome inspection)
- 5 primary families, leave-one-family-out CV: OXDA (EC 1.4.3.3, FAD oxidase, PF01266), HXK4 (EC 2.7.1.1, glucokinase, PF00349+PF03727), BLAT (EC 3.5.2.6, serine beta-lactamase, PF13354), CBS (EC 4.2.1.22, PLP lyase, PF00571+PF00291), PPM1D (EC 3.1.3.16, metal phosphatase, PF00481). AMIE (PF00795) frozen as backup-only, outside CV.
- Primary assay per family pre-declared; secondary assays of the same protein never cross holdout boundaries.
- Labels: ProteinGym binarized classes with per-assay frozen cutoffs. Semantics proven from data: bin==1 iff score>=cutoff (exhaustive check, zero exceptions). Raw scores never pooled across assays.
- Leakage: no shared Pfam domains (structural), pairwise identity screen all <0.30 (measured max 0.112), target sequence excluded from its own alignment, z-scores computed within training families only.
- Flexibility features (frozen definitions): site entropy, relative entropy, 20-dim substitution spectrum, PSSM log-likelihood ratio. Conservation baselines: column conservation, BLOSUM62, PSSM log-odds.
- Primary R1 question is the ablation: all-features must beat conservation-only AND profile-only with variant-bootstrap 95% CI excluding zero in >=4/5 families. Calibration gate + per-family abstention rules frozen. Useful-negative outcome explicitly in scope.

## Feasibility evidence (all live-sourced)
- ProteinGym v1 official HF repo (OATML-Markslab/ProteinGym_v1): 5 parquet shards, sha256 of shards 0-2 byte-identical to upstream LFS oids; 2.4M substitution rows total.
- Variant floors: primary assays carry 4,996-8,570 single substitutions each, both classes present (worst balance OXDA 1329/5440, still >250/class floor).
- Coordinates: 9/11 assay sequences byte-identical to UniProt canonical; OXDA exact substring 2..365. No ambiguous mapping anywhere.
- Homologs: 19k-250k UniProt entries per Pfam family (live counts).
- Licenses: ProteinGym code MIT + open data with citations; UniProt CC BY 4.0; BRENDA online CC BY 4.0 (audited; commercial DSI/Cali-Fund caveat flagged; not a dependency); Rhea direct API Cloudflare-blocked (limitation preserved; UniProt xrefs route used); InterPro API reachable.
- Environment: pyarrow 17.0.0 wheel added (sha256 f7ae2de6...), full ledger in environment_ledger.json. No torch; parquet streaming within 1GB RAM.

## Non-estimability preserved
Frozen clauses cover mapping failure, alignment column coverage <80%, and identity-screen failure. None triggered at R0.
