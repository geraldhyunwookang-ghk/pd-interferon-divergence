# Cell-class harmonisation

Fixed before any disease effect was estimated (plan §6.2). Each cohort's **published** annotation is used; nothing is re-clustered.

| Cell class | Kamath (SCP1768) | Smajić (GSE157783) | Martirosyan (GSE243639) |
|---|---|---|---|
| **Microglia — primary H1 class** | all `mg_UMAP` states **except** `Macro_CD200R1` | `cell_ontology == "Microglia"` | high-level `IDENT == "Micro"` |
| Microglia incl. macrophages — sensitivity | all `mg_UMAP` states | same as primary (no macrophage class) | same as primary (no macrophage class) |
| DA neurons — H2 | all `da_UMAP` subtypes | not used (74 nuclei in total) | not used (median 4 per PD donor) |
| DA neurons, SOX6 / CALB1 family | first token of `Cell_Type` | — | — |
| S-6 comparators | `nonda_UMAP` (primary), `astro_UMAP`, `olig_UMAP`, `opc_UMAP` | — | — |

## Decisions and reasons

1. **Macrophages.** Kamath separates a macrophage cluster; Smajić and Martirosyan do not, so their microglia may include border-associated macrophages. The primary H1 class excludes Kamath's macrophages; a sensitivity analysis includes them so that all three cohorts are compared like-for-like.
2. **Martirosyan uses the high-level `Micro` label (12,995 nuclei), not the microglia sub-cluster sheet (10,154 nuclei).** The high-level label matches the paper's Supplementary Table 2 exactly for all 29 donors; the sub-cluster sheet is missing ~2,840 nuclei without explanation.
3. **Smajić donor IDs.** The cell file uses `PD1`–`PD5`; GEO uses `IPD1`–`IPD5`. They are mapped by number. This is checked after loading expression data: sex-linked genes must identify `PD1` and `C1` as the two female donors (plan §4.3).
4. **Sorted fractions (Kamath).** DA and non-DA neurons come mainly from the NURR-positive fraction; glia mainly from the negative fraction. Effect sizes are compared only within a fraction (plan §4.2, §6.6).
