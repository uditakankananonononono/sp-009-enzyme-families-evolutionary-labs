# SP-009 / DOC-2-076 - Enzyme Families as Evolutionary Laboratories

## Summary (from `SP-009-DOC-2-076/README.md`)


## Status

Closed useful negative after R1. The evolutionary-flexibility add-on did not meet the frozen five-family graduation rule.

## Useful result

Only three of five primary families passed the alignment-source gate. Of those, CBS alone passed the both-baseline, both-bootstrap incremental-value rule. BLAT showed a partial gain that did not separate from the profile baseline under the position-block bootstrap. HXK4 showed a significant PR-AUC reversal against conservation. Family-specific base rates coupled to family-correlated features also inflated a within-family permutation null, demonstrating why raw leave-one-family-out AUC can overstate variant-level biological transfer.

## What is new

This project tested whether natural substitution-pattern features add transferable information beyond conservation and profile baselines across mechanistically distinct enzyme families. It combined pre-outcome alignment estimability gates, leave-one-family-out validation, site-block uncertainty, and a mechanistic null-anomaly investigation.

## Why it matters

The negative prevents an unreliable evolutionary-flexibility score from being promoted into enzyme-engineering mutation prioritization. The durable application is a benchmark/audit workflow: verify alignment depth and coverage, use position-block uncertainty, compare against simple family-profile baselines, and test family-shift base-rate leakage before claiming cross-family skill.

## Claim boundary

No enzyme-engineering product or catalytic-activity claim graduates. Best simple conservation/profile baselines remain the operational reference.

## Top-lab next question

Would a prospectively selected cohort restricted to alignment-estimable enzyme families reproduce the small CBS signal on fresh families, or does family-shift structure explain all apparent incremental value? This is not an automatic R2 because the current method failed its frozen rule.

## Contents

- `protocol/`: frozen R0 and R1 protocols and lock records
- `code/`: retrieval, alignment, modeling and null-diagnostic code
- `results/`: source gates, model results, uncertainty and diagnostic artifacts
- `report/`: R0 and R1 technical reports
- `environment_ledger.json`: compute and dependency record

## Drive artifact

https://drive.google.com/file/d/1B8CeMKTI2yxfipX4p-mC8Jw2Wn5MOWvT/view?usp=drivesdk&authuser=uditakankana%40gmail.com

## Contents

- `SP-009-DOC-2-076/` - migrated unchanged from `science-program/projects/SP-009-DOC-2-076` (46 files)

## Provenance

Split out of the `science-program` repository (source commit `028a7141ed5f951a7b6e6517d4e72768d414a560`) on 2026-09-23. Every file is byte-identical to the source; `MIGRATION_MANIFEST.tsv` lists sha256, original path and new path for each of the 46 files.

Part of Udita Phookan's computational science program: every experiment locks its question, validation design, success gate and failure policy before outcome analysis, and negative results are preserved. Program-wide ledgers and standards live in the `science-program-ledger` repository.
