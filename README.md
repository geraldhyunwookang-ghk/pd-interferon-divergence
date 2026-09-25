# pd-interferon-divergence

**Testing cell-type-divergent Type-I interferon transcription in the Parkinsonian substantia nigra.**

A pre-registered reanalysis of public single-nucleus RNA-seq data from post-mortem human substantia nigra. It asks whether, in Parkinson's disease, **microglia show increased interferon-responsive transcription** while **dopaminergic neurons show reduced expression of the machinery that receives interferon signals**.

- **Microglial arm (H1):** tested in three independent cohorts — Kamath 2022, Smajić 2022, Martirosyan 2024.
- **Neuronal arm (H2):** testable only in Kamath 2022; reported as **unreplicated**.
- **Unit of analysis:** the donor (pseudobulk). No statistical test is performed on individual cells.

## Pre-registration

`PRE_ANALYSIS_PLAN_v2.0.md` was committed **before any expression data were loaded**. The commit timestamp is the pre-registration date. Any later change is logged in `DEVIATIONS.md`.

## Repository

| Path | Contents |
|---|---|
| `PRE_ANALYSIS_PLAN_v2.0.md` | Frozen analysis plan |
| `DEVIATIONS.md` | Post-freeze changes |
| `harmonisation.md` | How cell types are matched across cohorts |
| `audit/` | Donor table, source checksums, and the metadata-only audit script |
| `genesets/` | Every gene set used, with MSigDB version |
| `scripts/` | Analysis scripts, run in numbered order |

## Data

All data are public: Broad Single Cell Portal SCP1768; GEO GSE157783 and GSE243639. See `audit/sources.csv` for exact files and SHA-256 checksums. Raw data are not redistributed here.

## Author

Gerald Kang, University of Cambridge.
