# Pre-Analysis Plan

## Testing cell-type-divergent Type-I interferon transcription in the Parkinsonian substantia nigra

| | |
|---|---|
| **Author** | Gerald Kang, University of Cambridge |
| **Version** | 2.0 |
| **Date frozen** | The timestamp of the Git commit that adds this file to the repository |
| **Status** | **FROZEN** — any later change is recorded in `DEVIATIONS.md` |
| **Repository** | `github.com/geraldhyunwookang-ghk/pd-interferon-divergence` |

> **Freeze statement.** This plan was committed before any individual-level expression matrix was loaded and before any disease-group effect was estimated. Pre-freeze inspection is limited to donor metadata, nucleus-annotation counts, gene-set membership, and published aggregate marker summaries (disclosed in §4.4). No donor-level expression value has been inspected. Every change after the freeze is recorded in `DEVIATIONS.md` with its date, reason, and whether it was made before or after the relevant result was seen.

---

### The study in five lines

1. **Question.** In Parkinson's disease (PD), do microglia show *increased* interferon-responsive transcription while dopaminergic (DA) neurons in the same donor cohort show *reduced* expression of the machinery that receives interferon signals?
2. **Data.** Public single-nucleus RNA-seq from post-mortem human substantia nigra. Three cohorts for the microglial question; one for the neuronal question.
3. **Method.** One summed ("pseudobulk") profile per donor per cell type; a pathway score per donor; a linear model adjusting for age, post-mortem interval and sex. Donors, never cells, are the unit of analysis.
4. **Primary tests.** H1 (microglia, three cohorts pooled) and H2 (DA neurons, one cohort), Holm-corrected.
5. **Honest scope.** H2 cannot be externally replicated with current public data and will be labelled **UNREPLICATED** in every output.

---

## 0. What changed in v2.0, and why

v2.0 is the result of a line-by-line audit of v1.9 against primary sources. The architecture is unchanged. The corrections fall into four groups. Full version history is in Appendix B.

### 0.1 Factual and citation errors

| # | Error in v1.9 | Correction | Verified against |
|---|---|---|---|
| F1 | Human single-cell IFN-I study cited as "Yan Y, *et al.*, 22:172" | First author is **Quan P**; article **22:212** | *Cell Commun Signal* article page and PubMed 38566100 |
| F2 | Hinkle 2022 cited with the title of a *commentary* about it | Correct title: "STING mediates neurodegeneration and neuroinflammation in nigrostriatal α-synucleinopathy", PNAS 119:e2118819119 | PNAS; Johns Hopkins record |
| F3 | Hinkle 2022 presented as uncontested | Its central finding — STING-deficient mice protected from DA loss — **was not reproduced** by Klæstrup *et al.* 2026 | *npj Parkinsons Dis* 12:17 |
| F4 | Sliter 2018 retraction reason paraphrased loosely | Retracted because Fig. 3a–c data were discrepant with raw data and Fig. 1d was not reproducible | *Nature* 644:1116 (2025) |
| F5 | Villanueva 2026 cited as "33:57" | Article number not verifiable; cited by DOI | Springer article page |
| F6 | "Fixed 21-gene canonical core" | The list contains **20** genes | Count |
| F7 | Gene-set files referenced as `genesets/manifest.yml` | Actual files are `genesets/membership.csv` and `genesets/provenance.txt` | Repository |
| F8 | "Only `OAS1` of the core survives; the other 19 are shared with IFN-γ" | **Computed (C3): only 12 are shared; 5 sit in the IFN-γ set only and 2 in neither.** This exposed a second problem — Axis E contained canonical Type-I genes — now fixed (§5.6) | `00_genesets.R` output |

### 0.2 Omissions in the evidence base

| # | Omission | Consequence |
|---|---|---|
| E1 | The neuronal-protective arm (Ejlerskov 2015, Tresse 2021, Magalhaes 2021, Villanueva 2026) comes **entirely from one laboratory** | Declared in §2.4 and limitation 14 |
| E2 | The IFNAR-dependent glial-damaging arm (Main 2016, 2017) also comes predominantly from one laboratory, as does the review used to frame the field (Chen 2025) | Declared in §2.4 |
| E3 | Cytosolic mtDNA release is **not established in sporadic PD**; cell-free mtDNA is elevated in PRKN/PINK1 carriers but not in sporadic PD (Borsche 2020) | Plain-English summary corrected |

### 0.3 Statistical corrections

| # | Problem in v1.9 | Correction |
|---|---|---|
| S1 | **H3 was said to "cost no additional alpha" but did not.** Kamath-H1 was tested at nominal 0.05, outside the Holm family, creating an extra route to a false claim | Both H3 components now require p < 0.025. H3 is a separately controlled secondary claim (§6.5) |
| S2 | **The plan's own covariate rule required PMI in the primary model, but the model omitted it.** PMI is complete in every cohort | PMI enters the primary model. Fallback threshold changed from 8 to 7 residual df, because the old threshold would have removed PMI from Smajić — the cohort where PMI is most imbalanced (§6.3) |
| S3 | The population used to standardise genes was unspecified, so LBD donors could influence the PD-vs-control estimate through z-scoring | Standardisation is performed within each model's own analysis set (§6.3) |
| S4 | **The primary H1 pooled test used Hartung–Knapp random effects with k = 3**, which uses a t-distribution with 2 df (critical value 4.30) and has very low power | Primary pooled test is now common-effect inverse-variance; Hartung–Knapp random effects is a mandatory sensitivity analysis. Protection against a single-cohort-driven result comes from the concordance criterion (§6.9) |

### 0.4 Design gaps

| # | Gap | Correction |
|---|---|---|
| D1 | **Wang was excluded for lacking deposited nucleus-level annotations, but Smajić and Martirosyan were never checked against the same standard** | Blocking pre-freeze check C1 |
| D2 | Kamath's microglia file includes a macrophage cluster (`Macro_CD200R1`); v1.9 did not say whether H1 includes it | Primary H1 excludes macrophages; inclusion is a sensitivity analysis (§6.2) |
| D3 | S-6 (cell-type specificity) required astrocyte, oligodendrocyte, OPC and non-DA annotations the pipeline did not load; it also compared effect sizes across sorted fractions and at unequal depth | Annotations loaded (C4); S-6 redesigned with non-DA neurons as a same-fraction primary comparator and depth matching (§6.6) |
| D4 | E2 compared GPNMB states against "homeostatic" microglia, which Kamath's annotation never defines | E2 redefined; feasibility computed (§3.4) |
| D5 | Overlaps between axes (e.g. `STAT2`, `IRF9`, `USP18` possibly in Hallmark IFN-α) not checked | Pre-freeze check C3 |
| D6 | Symbol aliases incomplete (`GBA1`/`GBA`, `PRKN`/`PARK2`, `RIGI`/`DDX58`, `TENT5A`/`FAM46A`) | Alias table added (§5.1) |
| D7 | Stale references to "v1.7" in several sections; duplicated headings in the checklist | Removed |

---

## 1. Plain-English summary

Interferon is the body's antiviral alarm. Cells release it when they detect something that looks like infection — including, in laboratory models, DNA leaking from damaged mitochondria. Whether that particular trigger operates in the brains of people with ordinary (sporadic) Parkinson's disease has not been established.

Two groups of experiments point in opposite directions. In one, blocking interferon signalling in mice protects dopamine neurons from toxins, suggesting interferon made by immune cells is harmful. In the other, removing interferon or its receptor specifically from neurons makes those neurons degenerate, suggesting neurons need a baseline level of interferon signalling to keep their mitochondria healthy.

These need not contradict each other: they manipulate different cells. This study asks whether both patterns are visible, in different cell types, in brain tissue donated by people with and without Parkinson's disease. Specifically: do microglia show **more** interferon-responsive gene activity, while dopamine neurons show **less** of the machinery needed to receive interferon signals?

This is observational. It can detect whether the pattern is present. It cannot show that interferon causes or prevents neuronal death.

---

## 2. Background and evidential status

### 2.1 Evidence that glial IFN-I signalling is damaging

- **IFNAR-dependent toxicity.** In the MPTP model, *Ifnar1*<sup>−/−</sup> mice show attenuated neuroinflammation and reduced dopaminergic loss, reproduced with a blocking anti-IFNAR1 antibody (Main 2016). In co-culture, IFN-I produced by **glia** drives rotenone-induced neuronal death, and the effect depends on IFNAR1 in the glia rather than the neurons (Main 2017). This is the most direct prior evidence localising a damaging arm to the glial compartment.
- **Chronic neurodegeneration.** In murine prion disease, a STING-dependent IFN-β response arises in microglia and IFNAR1 deficiency slows progression (Nazmi 2019).
- **STING-specific evidence is contested.**
  - Sliter 2018 reported STING-dependent inflammation and dopaminergic loss in PINK1/Parkin models. **It was retracted in 2025** because the data in Fig. 3a–c were discrepant with the raw data and Fig. 1d could not be reproduced. It carries no evidentiary weight here.
  - Hinkle 2022 reported that STING-deficient mice were protected from α-synuclein-fibril-induced dopaminergic loss. **Klæstrup 2026 did not reproduce this**: using the same strains and injection site, STING-deficient mice showed faster but similar-degree degeneration at six months. The authors discuss differences in fibril dose and preparation and in time points.
- **Human single-cell evidence.** In GSE157783, IFN-I activity scores were highest in microglia, with NFATc2 proposed as an upstream regulator (Quan 2024). That analysis scored individual cells rather than donors.

### 2.2 Evidence that neuronal IFN-I signalling is protective

- **Mitochondrial maintenance.** Neuronal IFN-β acts through STAT5 and PGAM5 to phosphorylate DRP1, supporting mitochondrial fission; *Ifnb*<sup>−/−</sup> neurons accumulate damaged mitochondria, and IFN-β rescues dopaminergic death in other PD models (Tresse 2021). *Ifnb*<sup>−/−</sup> mice develop spontaneous neurodegeneration with Lewy-body-like pathology (Ejlerskov 2015).
- **Cell-type-specific deletion.** Neuronal *Ifnar1* deletion alone induces nigral DA neuron loss and progressive motor and cognitive deficits, with defective mitophagy reversible by IFNAR1 overexpression (Villanueva 2026).
- **Human expression, read precisely.** Villanueva's human analysis displayed scaled mean expression per donor from the same Kamath dataset used here (n = 3–8 per diagnostic group). Their reported pattern:
  - *IFNAR1* lower in **Lewy body dementia (LBD)** than in both controls and PD, in DA neurons and microglia.
  - *IFNB1* lower in DA neurons in both LBD and PD.
  - *IFNAR1* and related genes also lower in astrocytes, oligodendrocytes and OPCs.

> **Three consequences, recorded before analysis.**
>
> **(i)** The reported receptor reduction is **LBD-predominant**. H2's primary contrast is PD versus control — the comparison where Villanueva's descriptive data showed the least change. H2 is therefore a genuinely open test. An LBD contrast is pre-specified as secondary (§3.2) so that a null PD result cannot later be replaced by it.
>
> **(ii)** The reduction is reported across **many cell types**. A DA-specific interpretation is available only if S-6 (§6.6) shows it.
>
> **(iii)** Microglial reception machinery (Axis A) may fall alongside neuronal Axis A. That would be concordant, not divergent — but H1 concerns interferon-*responsive transcription* (Axis B), not reception. The divergence reading depends on the A/B distinction being real, so microglial Axis A is reported alongside H1.

### 2.3 The compartment-specific model

| | Damaging arm | Protective arm |
|---|---|---|
| Compartment manipulated | Global, or glia-restricted | Neuron-restricted |
| Manipulation | Pathway **blocked** → protection | Pathway **removed** → damage |
| Models | MPTP, rotenone, prion | *Ifnb*<sup>−/−</sup>, neuronal *Ifnar1* deletion, 6-OHDA |

Read together, these studies are consistent with a model in which **induced IFN-I-responsive activity in glia accompanies damaging inflammatory states, while tonic IFN-I signalling in neurons supports mitochondrial maintenance.** That model has been inferred across separate experimental systems. It has not been examined in both compartments of unmanipulated human tissue.

### 2.4 Evidential status of each arm

| Claim | Principal evidence | Independence | Status |
|---|---|---|---|
| Glial IFN-I contributes to DA damage (IFNAR-dependent) | Main 2016, Main 2017 | Predominantly one laboratory | Supported in toxin models; not independently replicated in PD models to our knowledge |
| …specifically via STING | Sliter 2018; Hinkle 2022 | Two groups | **Contested**: one retracted, one not reproduced |
| Microglial IFN-I activation in human PD | Quan 2024 | One study | Cell-level inference; donor-level test not yet reported |
| Neuronal IFN-I/IFNAR required for DA maintenance | Ejlerskov 2015; Tresse 2021; Magalhaes 2021; Villanueva 2026 | **All one laboratory** | Internally consistent; no independent replication identified |

This table matters for how results are interpreted. Because **each arm rests largely on a single laboratory's work**, a donor-level test in human tissue by an unrelated analyst is an independent line of evidence for both. That is part of this study's value — and it also means a null result would be informative rather than merely disappointing.

### 2.5 The measurement distinction

- **Axis A — reception and transduction.** Transcript abundance of receptors and signal transducers required to *receive* an interferon signal.
- **Axis B — interferon-responsive transcription.** Transcript abundance of interferon-stimulated genes, consistent with a response to interferon or related inflammatory signalling.

**Wording rule, binding throughout.** Axis A measures the expression of signalling components, not receptor function. Axis B measures a transcriptional state, not ligand concentration or secretion. The hypothesis is therefore stated as:

> PD is associated with divergent transcriptional states of the Type-I interferon pathway across cell types: **increased interferon-responsive transcription in microglia** alongside **reduced expression of interferon reception/transduction machinery in dopaminergic neurons**.

### 2.6 Novelty, stated against the closest prior work

**H1 alone is not novel.** Quan 2024 already analysed GSE157783 — this study's second H1 cohort — for IFN-I activity across cell types. **The gene panel is not novel either.** Villanueva 2026 displayed five of the seven Axis A genes, in the same cells, from the same Kamath donors.

The contribution is:

1. **Testing co-occurrence.** Whether microglial interferon-responsive transcription is accompanied by an opposing neuronal reception phenotype in the same donor cohort. Neither prior study could address this: GSE157783 contains too few DA neurons, and Villanueva made no formal paired contrast.
2. **Inference at the right unit.** Donor-level modelling with covariate adjustment, confidence intervals, explicit correction for the ten-fold difference in recoverable DA nuclei between groups, and pre-registration. Quan scored individual cells; Villanueva displayed descriptive means.
3. **Replication by hypothesis.** H1 tested in every public cohort that can support it.
4. **An independent line of evidence** for two arms that each rest largely on one laboratory (§2.4).

The contribution is inferential and methodological, not the choice of genes, and will be described in exactly those terms.

### 2.7 Relation to the author's prior work

The author previously identified *IFNE* as recurrently deleted in glioblastoma and showed that restoring its expression suppressed proliferation *in vitro*. That work motivated an interest in context-dependent interferon biology. It does not constrain this study: *IFNE* is in no tested gene set, and both directions are tested on both axes.

---

## 3. Research questions and hypotheses

### 3.1 Primary hypotheses

| ID | Hypothesis | Estimand | Predicted sign | Cohorts | Supported if |
|---|---|---|---|---|---|
| **H1** | Microglial interferon-responsive transcription is higher in PD | β<sub>diagnosis</sub>, Axis B score, microglia | > 0 | Kamath, Smajić, Martirosyan | §6.9 criterion met and Holm-adjusted pooled p rejects |
| **H2** | DA-neuronal interferon reception/transduction machinery is lower in PD | β<sub>diagnosis</sub>, Axis A score, DA neurons | < 0 | Kamath only | Holm-adjusted p rejects with β < 0 |

**H3 — the divergence claim (secondary, Kamath only, unreplicated).** Satisfied if and only if, in Kamath, the H1 contrast has p < 0.025 with β > 0 **and** H2 has p < 0.025 with β < 0.

- Because it requires both components to reject, H3's own Type I error is at most 0.025 by the intersection–union principle.
- It is a **separate secondary claim**, not part of the primary familywise error statement.
- It establishes two opposite-signed group-level contrasts in cell classes from the same donor cohort. It does **not** establish that the two states covary within individual donors, their temporal order, or causation.

> **Why not a single divergence score?** A per-donor index D = z(Axis B, microglia) − z(Axis A, DA) was considered and rejected as the test of H3. D rises if microglial activity increases while neuronal machinery stays flat, so a significant D cannot show that *both* arms moved. D is plotted descriptively (Figure 5) with no p-value.

### 3.2 Pre-specified secondary contrasts (Kamath only; nominal p-values; labelled *secondary*)

- **H2-SOX6.** H2 restricted to SOX6-family DA neurons — the vulnerable population. Also used in sensitivity analysis S-2.
- **H2-LBD.** H2 as LBD versus control. Declared now because the prior human data are LBD-predominant (§2.2). **Never promoted to primary.** A null H2 with a significant H2-LBD is reported in exactly these words: *the pre-registered PD contrast was null; the secondary LBD contrast was not.*

*Why PD stays primary despite the prior.* The study concerns idiopathic PD, and only Kamath has usable LBD numbers. Choosing a primary contrast by its probability of significance is precisely what pre-registration exists to prevent.

### 3.3 RQ2 — what accompanies reduced neuronal capacity? (exploratory; gated)

**Gate.** RQ2 is evaluated only if H2 (i) passes its Holm-adjusted test, (ii) passes the detectability-parity check, (iii) meets the S-1 robustness rule, and (iv) is not reversed by S-2. If any condition fails, RQ2 is reported as *not evaluable* and no pattern analysis is shown.

| Pattern | Description | Marker |
|---|---|---|
| **P1** Receptor-proximal | Axis A effect carried by receptor-proximal genes | A1 effect ≈ A effect |
| **P2** Negative-regulator elevation | Axis D raised concurrently | Axis D β > 0 **and** Axis D without `USP18` β > 0 (§5.6b) |
| **P3** Mitochondrial/autophagic reduction | Mitophagy–fission and lysosomal sets reduced concurrently | Mito set β < 0 |

Patterns may co-occur and may be causally ordered; this design cannot distinguish order. Reported as estimates with confidence intervals in one table. **No model-selection statistic; no p-value declares a winner.** Elevated `SOCS1` does not demonstrate feedback shutdown, and reduced `PGAM5` does not demonstrate failed fission.

### 3.4 Exploratory analyses

- **E1 — vulnerability gradient.** Axis A disease effect in SOX6-family versus CALB1-family DA neurons. Only 7 control / 3 PD / 3 LBD donors carry ≥20 nuclei of both families: **underpowered, descriptive only**, shown as estimates with confidence intervals.
- **E2 — microglial states.** Axis B disease effect estimated separately within each microglial state. *(Redefined in v2.0: Kamath's annotation does not designate a "homeostatic" state, so no such comparator is used.)* Feasibility, computed pre-freeze at ≥20 nuclei per donor with ≥4 donors per group:

| State | Ctrl | PD | LBD | Feasible (PD vs Ctrl) |
|---|---|---|---|---|
| MG_GPNMB_LPL | 7 | 6 | 3 | Yes |
| MG_GPNMB_SULT1C2 | 8 | 6 | 4 | Yes |
| MG_GPNMB_SUSD1 | 7 | 6 | 4 | Yes |
| MG_CECR2_FGL1 | 8 | 6 | 4 | Yes |
| MG_FOSL2 | 8 | 6 | 4 | Yes |
| MG_OPRM1 | 8 | 6 | 4 | Yes |
| MG_SPON1 | 7 | 6 | 4 | Yes |
| MG_TSPO_VIM | 8 | 6 | 4 | Yes |
| MG_CCL3 | 4 | 1 | 1 | No |
| MG_MGAM | 4 | 3 | 2 | No |
| MG_MKI67 | 1 | 0 | 0 | No |
| Macro_CD200R1 *(macrophage)* | 8 | 6 | 4 | Excluded from E2 |

---

## 4. Datasets and pre-freeze audit

> Every count below was generated from deposited metadata or supplementary tables, not from published prose (§4.6 explains why). The script-generated `audit/donor_table.csv` is authoritative for analysis.

### 4.1 Cohort summary

| Cohort | Accession | Tissue | Preparation | Donors (Ctrl / PD / other) | H1 | H2 |
|---|---|---|---|---|---|---|
| **Kamath 2022** | SCP1768 / GSE178265 | SNpc | NURR-sorted + unsorted | 8 / 6 / 4 LBD | ✓ | ✓ |
| **Smajić 2022** | GSE157783 | Ventral midbrain | Unenriched | 6 / 5 / — | ✓ | ✗ |
| **Martirosyan 2024** | GSE243639 | SNpc | Unenriched | 14 / 15 / — | ✓ | ✗ |
| Wang 2024 | GSE184950 | SN | Unenriched | 9 / 6 / 17 PDD | ✗ | ✗ |

### 4.2 Kamath — primary cohort for both arms

Post-mortem SNpc; NR4A2/NURR-based enrichment for DA nuclei. Metadata and annotations from Broad SCP1768.

**Donors passing ≥20 nuclei per cell class:**

| Threshold | DA neurons (Ctrl / PD / LBD) | Microglia (Ctrl / PD / LBD) |
|---|---|---|
| **≥20 (primary)** | **8 / 5 / 4** | **8 / 6 / 4** |
| ≥50 | 8 / 4 / 3 | 8 / 6 / 4 |
| ≥100 | 7 / 3 / 3 | 8 / 6 / 4 |

Median DA nuclei per donor: 1,464 (Ctrl), **139 (PD)**, 884 (LBD). Median microglia per donor (macrophage cluster excluded, as in the primary H1 class): 2,093 / 1,351 / 1,579. Donor counts at every threshold are unchanged by the exclusion.

**H2 primary donor set — verified complete for age, sex and PMI:**

| Group | n | Ages | PMI (h) | Sex |
|---|---|---|---|---|
| Control | 8 | 49, 49, 79, 79, 82, 90, 91, 92 | 7.0–23.3 | 6F / 2M |
| PD | 5 | 76, 78, 82, 83, 90 | 2.0–22.3 | 2F / 3M |

One PD donor (15 DA nuclei in total) is excluded by the threshold; another contributes 26.

**Recorded before analysis:**
- **The neuronal arm is threshold-fragile; the microglial arm is not.** At ≥100, H2 retains 3 PD donors. S-3 is expected to show H2 attenuating; attenuation is neither evidence for nor against H2.
- **Two controls are aged 49**, against a PD group aged 76–90. Age is adjusted for, but this is a real imbalance in the smallest analysis.
- **Sorting confound.** DA nuclei derive ~97% from the NURR-positive fraction and microglia ~76% from the negative fraction. Each disease-vs-control comparison is within a fraction and valid, but **effect sizes cannot be compared between compartments**. H3 concerns direction only.
- **LBD handling.** The four LBD donors are excluded from all primary models. They enter only H2-LBD and the PD+LBD sensitivity analysis.
- No formal power calculation has been done. Three thousand microglia from six donors is an effective n of six.

### 4.3 Smajić — H1 replication

Unenriched ventral midbrain; 11 donors; ~41,000 nuclei.

**Donor split resolved from GEO metadata: 5 idiopathic PD (`IPD1`–`IPD5`), 6 controls (`C1`–`C6`).** The paper's abstract states the reverse; its Methods and the deposited metadata agree with each other. Quan 2024 independently analysed the same dataset as five PD and six controls.

| | n | Age, mean (range) | PMI h, mean (range) | Sex |
|---|---|---|---|---|
| IPD | 5 | 77.4 (66–84) | 22.2 (13–25) | 1F / 4M |
| Control | 6 | 83.0 (66–93) | 15.3 (5–29) | 1F / 5M |

**Recorded before analysis:**
1. **PMI is ~7 h longer in cases.** This is why PMI is in the primary model (§6.3).
2. **Controls are older than cases — the reverse of Kamath.** A concordant H1 effect across both cohorts would therefore be unlikely to reflect age.
3. **Sex is 2F / 9M overall.** Sex is dropped here under the fallback rule.
4. **Cannot test H2.** Fewer than 200 DA neurons across the whole dataset (Martirosyan 2024).
5. **Not an independent discovery.** Quan 2024 analysed this dataset for IFN-I; the result here is a donor-level re-test.

**C1 and C2 — PASSED.** `GSE157783_IPDCO_hg_midbrain_cell.tsv` gives a cell-type label (`cell_ontology`) and a donor (`patient`) for every one of 41,435 nuclei across all 11 donors.

| | Microglia per donor | Median | ≥20 | ≥50 | ≥100 |
|---|---|---|---|---|---|
| Control (C1–C6) | 192, 206, 109, 321, 140, 237 | 199 | 6/6 | 6/6 | 6/6 |
| PD (PD1–PD5) | 207, 653, 556, 129, 1,153 | 556 | 5/5 | 5/5 | 5/5 |

**No threshold attrition at any level.** The file labels only **74** nuclei as DA neurons across all 11 donors, confirming that Smajić cannot test H2.

**Three decisions, fixed now:**
- **Donor ID mapping is an assumption, and will be checked.** The cell file uses `PD1`–`PD5`; GEO uses `IPD1`–`IPD5`. Covariates come from GEO, so `PDk` is mapped to `IPDk` by number. **Pre-specified QC check (not a disease effect):** after loading the expression matrix, sex-linked expression (*XIST*; Y-chromosome genes such as *RPS4Y1*, *DDX3Y*) must identify `PD1` and `C1` as the two female donors, as GEO records. If it does not, the mapping is reported as unverified and Smajić's covariate-adjusted estimate is replaced by an unadjusted one, clearly labelled.
- **No separate macrophage class exists** in this annotation. As with Martirosyan, Smajić's `Microglia` is compared against Kamath's including-macrophage class in the sensitivity analysis.
- **PD donors contribute about 2.8 times more microglial nuclei than controls** (median 556 vs 199). This could reflect microgliosis or dissection, and it means PD pseudobulks are deeper. Detectability parity (§6.1) is therefore **also run for Axis B in microglia**, in every H1 cohort.

### 4.4 Martirosyan — H1 replication

SNpc; 15 PD / 14 controls, confirmed from the GEO series matrix; per-donor cell counts from Supplementary Table 2.

| | n | Age, mean (range) | PMI h, mean (range) | RIN, mean (range) | Sex |
|---|---|---|---|---|---|
| PD | 15 | 82.3 (57–99) | 23.3 (3–61) | 6.99 (6.3–7.8) | 4F / 11M |
| Control | 14 | 75.5 (30–93) | 17.2 (3.25–52) | 7.24 (6.6–8.0) | 5F / 9M |

**Microglia per donor:** median 469 (Ctrl) and 364 (PD); all 29 donors pass ≥20 and ≥50; one control falls below ≥100. **No threshold attrition.**

**DA neurons per donor:** median **4** in PD (range 0–276); 5 of 15 PD donors reach ≥20. The paper's aggregate of "more than 2,000 DA nuclei" is accurate but concentrated in a few control donors. **H1 only.** Within the same tissue and protocol, microglial yield per PD donor exceeds DA yield roughly 91-fold.

**Recorded before analysis:**
- **One control is aged 30** in a cohort otherwise 61–93. A sensitivity analysis excluding this donor is reported alongside the full estimate; neither is designated primary afterwards.
- **PMI extends to 61 h.**
- **RIN is available only here.** It enters a Martirosyan-specific sensitivity analysis, never the pooled model, because differently specified models cannot be pooled coherently.

**C1 — PASSED.** `GSE243639_UMAP_coordinates.xlsx` deposits a cell-type label for every one of 83,484 nuclei across all 29 samples. The sample prefix of each nucleus ID (e.g. `s.0096`) matches the Sample ID column of Supplementary Table 2, giving nucleus → donor → diagnosis. Per-donor microglia counts derived from this file are **identical to Supplementary Table 2 for all 29 donors**.

**Two harmonisation decisions, fixed now:**
- **H1 uses the high-level `Micro` label** (12,995 nuclei). The file's separate microglia sub-clustering sheet contains only 10,154 nuclei; the ~2,840 missing nuclei were evidently removed during sub-clustering, for reasons not documented. Using the high-level label keeps H1 consistent with Supplementary Table 2.
- **The high-level labels have no separate macrophage class** (classes: Oligo, Astro, Micro, OPC, Neurons, VC, T cells). Martirosyan's `Micro` may therefore include border-associated macrophages. It is matched to Kamath's *including-macrophage* microglia in a pre-specified sensitivity analysis, and the discrepancy is recorded in `harmonisation.md` (limitation 15).

### 4.5 Wang — excluded from both arms

Table S1 lists 32 sequenced donors: **9 controls, 6 PD, 17 PD dementia** — not the "23 idiopathic PD" of the article text. The published cluster assignments are not deposited; GEO holds per-donor raw matrices only. Re-deriving the clustering would break the inherited-annotation rule (§6.2). Wang may enter a future version only if author-derived nucleus-to-cluster assignments are obtained before its disease effects are inspected.

*Conditional rule, fixed now in case that happens:* published marker summaries (Fig. S10B) show c9 expresses *TH* in 5.0% of cells, versus 29.6% (c6_2) and 54.4% (c7_3). The harmonised DA class would be c6_2 + c7_3; c9 would enter only as a sensitivity analysis. These aggregate percentages are the only expression information inspected before freezing.

### 4.6 A reportable observation: prose versus deposited metadata

In **every** cohort examined, published prose misdescribed the deposited metadata in a way that mattered for eligibility or disease definition:

| Cohort | Published description | Deposited metadata |
|---|---|---|
| Smajić | Abstract: 6 PD / 5 control | Series matrix: 5 IPD / 6 control (Methods agree) |
| Martirosyan | A 2025 reanalysis: 11 male PD / 9 control | Series matrix: 15 PD / 14 control |
| Wang | Article text: 23 idiopathic PD / 9 control | Table S1: 6 PD / 17 PDD / 9 control |
| Wang | HCA portal: 21 IPD + 1 LBD + 9 control (n = 31) | Table S1 sequenced set: n = 32 |

This is not a criticism of the source studies. It is a methodological finding: **cohort composition for reanalysis should come from deposited sample metadata, never from prose — including the source paper's own abstract.** It is reported regardless of hypothesis outcome.

### 4.7 Replication matrix

| | Kamath | Smajić | Martirosyan | Wang |
|---|---|---|---|---|
| **H1** microglial Axis B | ✓ | ✓ | ✓ | ✗ annotations not deposited |
| **H2** DA-neuronal Axis A | ✓ | ✗ DA yield | ✗ DA yield per donor | ✗ annotations not deposited |
| **H3** conjunction | ✓ unreplicated | ✗ | ✗ | ✗ |

**H1 draws on 28 controls and 26 PD donors across three cohorts, with no threshold attrition in any cohort at any threshold (all verified).**

**Heterogeneity, declared.** The H1 cohorts differ in preparation (sorted-negative vs unenriched), dissection (SNpc vs ventral midbrain), tissue bank, chemistry and depth. Directional concordance is therefore treated as more informative than any pooled point estimate.

---

## 5. Gene sets

All sets are generated by `00_genesets.R` and written to `genesets/membership.csv` by `scripts/00_genesets.R`, with versions recorded in `genesets/provenance.txt`. Unmatched genes are logged, never silently dropped.

### 5.1 Symbol aliases

Published matrices use different gene-annotation vintages. Every set is matched on current and legacy symbols:

| Current (HGNC) | Legacy |
|---|---|
| `CGAS` | `MB21D1` |
| `STING1` | `TMEM173` |
| `RIGI` | `DDX58` |
| `PRKN` | `PARK2` |
| `GBA1` | `GBA` |
| `TENT5A` | `FAM46A` |

### 5.2 Axis A — IFN-I reception and transduction (H2 estimand)

| Set | Genes | Role |
|---|---|---|
| **A (primary)** | `IFNAR1, IFNAR2, JAK1, TYK2, STAT1, STAT2, IRF9` | H2 estimand |
| A1 — receptor-proximal | `IFNAR1, IFNAR2, JAK1, TYK2` | Sensitivity S-5 |
| A2 — ISGF3 transducers | `STAT1, STAT2, IRF9` | Descriptive |
| A3 — STAT5 branch | `STAT5A, STAT5B` | Descriptive only; never summed into A |
| A-γ — Type-II receptor | `IFNGR1, IFNGR2, JAK2` | Specificity comparison (§5.6) |

- **Why A3 exists.** The neuronal-protective mechanism (§2.2) runs through STAT5, not the canonical ISGF3 complex. Without A3, the plan would omit the transducer of the pathway H2 is meant to detect. With two genes it cannot support set-level testing and carries no hypothesis weight. A concordant A3-down / mitophagy-down pattern in DA neurons would be the most mechanistically specific result this design can produce — stated now so it cannot later be presented as an unexpected discovery.
- **Why A1 is the key sensitivity.** `STAT1`, `STAT2` and `IRF9` are transducers but are also themselves interferon-inducible, blurring Axes A and B. **Rule:** the phrase *reduced reception/transduction machinery* is used only if A1 is also negative. If A is negative but A1 is near zero or positive, the finding is described as driven by inducible transducers.
- **Ligands** (`IFNB1`, `IFNE`, IFN-α genes) are not in any tested set. They are reported in a per-gene detection table only; nuclear RNA and post-mortem tissue make them near-undetectable. `IFNE` is epithelially restricted and included solely for continuity with the author's prior work.
- **Acknowledged constraint.** Seven genes (four in A1) is a small set. The canonical IFN-I reception module is small; this is not padded with correlated interferon-stimulated genes.

### 5.3 Axis B — interferon-responsive transcription (H1 estimand)

- **Primary:** `HALLMARK_INTERFERON_ALPHA_RESPONSE` (97 genes, msigdbr 26.1.1), unmodified.
- **Sensitivity 1:** `REACTOME_INTERFERON_ALPHA_BETA_SIGNALING` *(confirmed present in msigdbr 26.1.1, C3)*.
- **Sensitivity 2:** a fixed **20-gene** canonical core — `ISG15, IFI6, MX1, MX2, OAS1, OAS2, OAS3, OASL, IFIT1, IFIT2, IFIT3, IFI44, IFI44L, BST2, XAF1, RSAD2, IFITM3, HERC5, EIF2AK2, SAMD9L`.

Only the Hallmark set enters primary correction and synthesis. Sensitivity sets cannot replace a null primary result. Axis B is never described as ligand production, secretion, or interferon "output".

### 5.4 Axes C and D (descriptive)

- **C1 — cytosolic DNA sensing:** `CGAS, STING1, TBK1, IRF3, IKBKE, TREX1`. The set of interest for the mtDNA hypothesis.
- **C2 — cytosolic RNA sensing:** `RIGI, IFIH1, MAVS`. A contrast that should *not* track the mtDNA hypothesis.
- **D — negative regulators:** `SOCS1, SOCS3, USP18, PTPN2, PIAS2`. `PIAS2` is included because PIAS2-mediated blockade of IFN-β signalling has been proposed in PD dementia (Magalhaes 2021).

### 5.5 Comparator programmes (RQ2)

- **Mitophagy / fission:** `DNM1L, PGAM5, PINK1, PRKN, OPTN, INF2, MFN1, MFN2, OPA1, FIS1, BNIP3, BNIP3L, SQSTM1, CALCOCO2`. `TBK1` is excluded because it belongs to C1; the overlap is reported. *`DNM1L` encodes DRP1, the terminal effector of the neuronal IFN-β → STAT5 → PGAM5 → DRP1 axis, and also appears among the top genes in a UCL Genomics pilot methylation study of levodopa-induced dyskinesia. That pilot was underpowered and blood-derived; the convergence is noted as hypothesis-generating only.*
- **Oxidative phosphorylation:** `HALLMARK_OXIDATIVE_PHOSPHORYLATION`.
- **Lysosomal / autophagic:** `KEGG_LYSOSOME` *(confirmed present, C3)* plus `GBA1, LAMP1, TFEB, CTSD, CTSB`.
- **Proteostasis:** `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`.

### 5.6 Where Type-I specificity can and cannot be tested

**Computed membership (msigdbr 26.1.1):**

| Set | n |
|---|---|
| Hallmark IFN-α response | 97 |
| Hallmark IFN-γ response | 200 |
| Shared | **73** |
| Alpha-only | 24 |
| IFN-γ minus IFN-α | 127 |
| **Axis E = IFN-γ minus IFN-α minus canonical core** | **122** |

**73 of the 97 IFN-α response genes are also IFN-γ response genes.** Type-I and Type-II interferon converge on STAT1, so the canonical interferon-stimulated programme is largely shared.

**Where the 20-gene canonical core actually sits (C3, computed):**

| Location | n | Genes |
|---|---|---|
| Shared by both Hallmark sets | 12 | — |
| IFN-α set only | 1 | `OAS1` |
| **IFN-γ set only** | **5** | **`MX2, OAS2, OAS3, IFIT1, XAF1`** |
| In neither Hallmark set | 2 | `IFI6, HERC5` |

*v2.0's earlier statement that the 19 non-`OAS1` genes were all shared was an unverified inference and was wrong; this table replaces it.*

**Two consequences, fixed now:**

- **Axis E as first defined contained canonical Type-I genes.** Subtracting the IFN-α set from the IFN-γ set leaves `MX2, OAS2, OAS3, IFIT1` and `XAF1` inside Axis E, although they are textbook Type-I-inducible genes. A contrast set meant to represent IFN-γ-associated transcription should not contain them. **Axis E is therefore redefined as IFN-γ minus IFN-α minus the 20-gene canonical core: 122 genes.** The five removed genes are listed in the methods.
- **The 20-gene core is not a subset of the primary Axis B set.** Seven of its genes (`MX2, OAS2, OAS3, IFIT1, XAF1, IFI6, HERC5`) are absent from Hallmark IFN-α. As a sensitivity analysis it is therefore a partly *independent* canonical interferon signature, which makes it a more informative robustness check, not a less valid one.

**Three consequences, fixed now:**

1. **Axis B measures interferon-responsive transcription, not Type-I-specific transcription.** No version of H1 supports a Type-I-selective claim.
2. **The 24 alpha-only genes are not a purified Type-I set.** They are a peripheral group (`CCRL2, CD47, CNP, CSF1, ELF1, GBP2, GMPR, HLA-C, IFITM1, LAMP3, LPAR6, MOV10, MVB12A, NCOA7, NUB1, OAS1, PARP9, PROCR, SAMD9, SELL, TENT5A, TMEM140, TRIM5, UBA7`) with weak canonical credentials. They are reported descriptively only.
3. **Specificity is testable only at the receptor level**, where ligand-binding chains are exclusive:

| Type-I exclusive | Shared | Type-II exclusive |
|---|---|---|
| `IFNAR1, IFNAR2, TYK2` | `JAK1, STAT1` | `IFNGR1, IFNGR2, JAK2` |

Axis A-γ is curated rather than taken from MSigDB, because Hallmark sets describe responses, not receptors.

**Interpretation rule, fixed now:**

| Pattern | Reported as |
|---|---|
| A1 reduced, A-γ null | **Type-I-receptor-selective** reduction — the strongest specificity claim available |
| A1 reduced, A-γ reduced less | Reduction with **Type-I predominance** |
| A1 and A-γ reduced comparably | **Interferon-receptor-wide** reduction; not Type-I-selective |
| Axis E null while Axis B rises | Interferon-responsive transcription without a broader IFN-γ-associated programme |
| Axis E rises with Axis B | Interferon-family or broad inflammatory response. Expected in part given the 73 shared genes; **uninformative** about Type-I selectivity rather than evidence against it |

Boundaries are judged on direction and confidence-interval overlap and reported transparently.

### 5.6b Overlaps between axes (C3, computed)

| Curated axis | Genes also in Axis B (Hallmark IFN-α) |
|---|---|
| A | `STAT2, IRF9` |
| A2 | `STAT2, IRF9` |
| **A1** | **none** |
| A3 | none |
| C2 | `IFIH1` |
| D | `USP18` |

These overlaps are expected — `STAT2`, `IRF9`, `USP18` and `IFIH1` are themselves interferon-inducible — but they matter wherever two axes are compared **within the same cell class**, because shared genes would make the axes correlate by construction.

**Rules, fixed now:**
1. **H1 uses the unmodified Hallmark set.** H1 and H2 are measured in different cell types, so the overlap cannot create a spurious H1–H2 relationship.
2. **Whenever Axis B is compared with Axis A or D in the same cell class** — microglial Axis A reported alongside H1 (§2.2 iii), and the RQ2 pattern table in DA neurons (§6.8) — Axis B is recomputed **excluding `STAT2, IRF9, USP18, IFIH1`** (*Axis B-excl*).
3. **A1 is confirmed as the clean capacity measure**: it shares no genes with Axis B. This strengthens its role in S-5 and in the receptor-level specificity test.
4. **RQ2 pattern P2 (negative-regulator elevation)** is described as such only if Axis D **excluding `USP18`** also rises; otherwise a rise in Axis D is reported as consistent with ordinary interferon induction.

### 5.7 Negative control — empirical null

For each cell class, 1,000 random gene sets are drawn, matched to the target set on size, mean pseudobulk expression and detection rate, and run through the identical pipeline.

1. Their p-values should be approximately uniform (Kolmogorov–Smirnov against U(0,1)).
2. The observed Axis A and B statistics are reported as a rank within this null, alongside the parametric p-value.

**Stop rule.** Analysis halts for debugging only if the null is **grossly non-uniform (KS p < 0.001)**. A single random set with p < 0.05 is the expected false-positive rate and triggers nothing.

Housekeeping genes (`ACTB, GAPDH, TBP, RPL13A, PPIA`) are a QC readout only, not a statistical null: their expression shifts with RNA integrity in post-mortem tissue.

---

## 6. Analysis plan

### 6.1 Gate 0 — can the genes be measured at all?

Run first; reported regardless of outcome.

A gene is **analysable** in a cell class if it exceeds 1 logCPM in ≥50% of donors **in each diagnostic group separately**. A set is analysable if **≥5** members pass; **≥3** for A1 and A-γ. A2 and A3 are descriptive and reported gene-wise regardless.

| If this fails… | …then |
|---|---|
| Axis A in DA neurons | H2 is **not evaluable** — never reported as null. A1 becomes the reported estimand if it passes its own gate; if not, the neuronal arm is abandoned and reported as an assay-sensitivity finding. H3 becomes unevaluable |
| Axis B in microglia | H1 not evaluable in that cohort |
| A-γ | Receptor-level specificity comparison not evaluable |
| Axis E | Reported gene-wise only |

**Detectability parity (mandatory).** Regress the per-donor number of detectable Axis A genes in DA neurons on diagnosis, with the covariates of §6.3. If detection itself differs by diagnosis, an apparent Axis A reduction may be a sparsity artefact. **Reported alongside H2 every time H2 is presented**, not in supplementary material. The same check is run for **Axis B in microglia in every H1 cohort**, because PD donors contribute more microglial nuclei than controls in Smajić (§4.3); there, a detection difference would bias *towards* H1.

The detectability map is Figure 1.

### 6.2 Pseudobulk construction

- Raw UMI counts are summed per **donor × cell class**. **The donor is the unit of analysis. No statistical test is ever performed at the level of individual cells.**
- **Cell classes:**
  - *DA neurons:* all; SOX6-family; CALB1-family.
  - *Microglia (primary):* all microglial states **excluding** `Macro_CD200R1`. Including it is a sensitivity analysis.
  - *Microglial states:* for E2.
  - *For S-6:* non-DA neurons, astrocytes, oligodendrocytes, OPCs.
- Each dataset's **published** annotation is used; no re-clustering. For Smajić and Martirosyan, whether their microglial labels include border-associated macrophages is recorded in `harmonisation.md` before any effect is estimated.
- **Inclusion:** ≥20 nuclei per donor × class (primary); 50 and 100 as sensitivity thresholds.
- Gene filtering with `edgeR::filterByExpr` within each class; TMM normalisation.
- **No donor is excluded for any reason other than the nucleus threshold or a missing covariate.** There is no outlier removal. This is fixed so no donor can be dropped because of its effect on a result.

### 6.3 Primary estimand — the donor-level pathway score

For each cohort × cell class × gene set, **within each model's own analysis set** (for example, PD and control donors only in the primary model):

1. Compute donor × gene logCPM from the TMM-normalised pseudobulk.
2. Standardise each gene across the donors in that analysis set (z-score). *(v2.0: standardising across a wider set — for example including LBD donors — would let donors outside the model influence its estimate.)*
3. Score = mean z across analysable set members: one number per donor.
4. Fit the donor-level linear model below and report **β<sub>diagnosis</sub>, its standard error and 95% CI**.

**Primary model:**

```
score ~ diagnosis + age + PMI + sex
```

**Why PMI is in the primary model.** The pre-specified rule was that a covariate enters the primary model if it is complete in every contributing cohort. The audit shows PMI is complete in all three cohorts, so the rule requires it. PMI is also imbalanced in Smajić (22.2 vs 15.3 h) and affects RNA integrity.

**Fallback rule, fixed now.** Each model must retain **at least 7 residual degrees of freedom**. If it does not, covariates are removed in this order: **sex first, then PMI**. Age is never removed. The model actually used is printed on every forest plot.

*Why the threshold is 7, not the 8 used in v1.x:* with PMI now in the model, a threshold of 8 would have removed PMI from Smajić — the cohort where it is most imbalanced. The change was made before any expression data were seen.

**Resulting models, computed in advance:**

| Analysis | Donors | Model | Residual df |
|---|---|---|---|
| H1, Kamath | 8 + 6 | full | 9 |
| H1, Smajić | 6 + 5 | sex dropped | 7 |
| H1, Martirosyan | 14 + 15 | full | 24 |
| H2, Kamath | 8 + 5 | full | 8 |
| H2-LBD | 8 + 4 | full | 7 |
| S-2 (adds SOX6 fraction) | 8 + 5 | full + SOX6 fraction | 7 |
| S-3 at ≥50 (H2) | 8 + 4 | full | 7 |
| S-3 at ≥100 (H2) | 7 + 3 | sex and PMI dropped | 7 |


**Other covariates.**
- **RIN:** Martirosyan-specific sensitivity analysis only.
- **log10(nuclei count):** excluded from all primary models. In DA neurons it is a consequence of the disease and nearly collinear with diagnosis; adjusting for it would adjust away part of the effect being measured. Reported as a sensitivity analysis.

**Interpretation of β.** Because genes are standardised across the analysis set, the denominator includes between-group variation. β is therefore a between-donor standardised difference but **not Cohen's d**, and its scale depends partly on each cohort's case–control balance. Hedges' g is reported as a scaling sensitivity (§6.9).

**Diagnosis contrasts.** Primary: idiopathic PD vs control, in every cohort. Sensitivity (Kamath only): PD + LBD vs control — reported alongside, never substituted.

### 6.4 Supporting analyses

Run in parallel as corroborating evidence, not as the primary estimand:

- Gene-level `limma-voom` with sample quality weights, same design.
- `limma::camera` (competitive) and `limma::fry` (self-contained) gene-set tests per cohort.
- Per-gene forest plots for every Axis A gene in DA neurons, so a reader can see whether an effect is broad or driven by one gene.

Agreement with the primary score strengthens confidence. Disagreement is reported, never resolved by choosing the more favourable result.

### 6.5 Multiple testing

| Claim | Test | Threshold |
|---|---|---|
| **H1** | Pooled Axis B test across completed cohorts (§6.9) | Holm, familywise α = 0.05 |
| **H2** | Kamath Axis A test | Holm, familywise α = 0.05 |
| **H3** | Kamath H1 contrast **and** H2, predicted opposite signs | **Both p < 0.025** |
| Everything else | — | Nominal, uncorrected, labelled *exploratory* or *secondary* |

Holm correction applies to exactly the two primary p-values.

**H3's error control, stated correctly** *(v2.0)*. v1.9 claimed H3 "costs no additional alpha" while testing its microglial component at 0.05 outside the Holm family. That created a separate route to a false claim. H3 now requires each component to reach p < 0.025, so by the intersection–union principle H3's own Type I error is at most 0.025. It is reported as a separately controlled secondary claim.

If only one H1 cohort has been analysed when preliminary results are shown, its p-value is labelled *interim* and does not replace the planned pooled test.

### 6.6 Mandatory sensitivity analyses

All are reported whether or not H2 rejects.

**S-1 — Sparsity matching.** Control donors have roughly ten times more DA nuclei than PD donors (median 1,464 vs 139). Sparser pseudobulk detects fewer low-abundance transcripts, which could mimic a real reduction.
- *Procedure.* Let T = the median DA nucleus count among PD donors passing ≥20. **Computed from annotations: PD donors contribute 26, 77, 201, 436 and 1,960 DA nuclei, so T = 201.** Every included donor, in either group, with more than T nuclei is randomly downsampled to T; donors with 20 to T nuclei are kept unchanged. Pseudobulk is rebuilt and H2 re-run, with z-scores recomputed each time, across **100 fixed, recorded random seeds**.
- *H2 passes S-1 if* (i) the median resampled β is negative, (ii) ≥80% of iterations are negative, and (iii) the absolute median resampled β is ≥50% of the original. Iteration-level significance is not required.
- *Failure* is reported as **unresolved from sparsity**, not as support.

**S-2 — Subtype composition.** SOX6-family neurons degenerate preferentially, so the PD "all DA" pseudobulk contains a different mixture of subtypes than the control one.
- *Outputs:* (i) per-donor SOX6/CALB1 proportions; (ii) H2 refitted with SOX6 fraction as a covariate; (iii) H2-SOX6, which is composition-controlled by construction.
- *H2 passes S-2 if* neither adjusted estimate reverses sign and at least one retains ≥50% of the primary absolute β. The SOX6-fraction model is small and may be collinear; its condition number and variance-inflation factors are reported, and H2-SOX6 carries more weight if it is unstable.

**S-3 — Nucleus threshold.** H1 and H2 at ≥20, ≥50 and ≥100. A stability description, not a test.

**S-4 — Ambient RNA.** Cell-type-exclusive markers checked in the wrong compartment (`TH`, `SLC6A3` in microglia; `P2RY12`, `CSF1R` in DA neurons). Contamination biases towards the null, but its size should be visible.

**S-5 — Receptor-proximal only.** Axis A recomputed as A1 (§5.2).

**S-6 — Cell-type specificity of the Axis A effect.** Villanueva reports interferon-pathway reductions in several cell types. If Axis A falls everywhere, "reduced neuronal interferon capacity" is the wrong description, and a technical explanation (RNA quality, global transcriptional decline) becomes as plausible as a biological one.

*Feasibility (C4, computed).* Every comparator class is available for all 18 Kamath donors:

| Class | Nuclei | Donors ≥20 (Ctrl / PD / LBD) | Median per donor (Ctrl / PD / LBD) | Sorted fraction |
|---|---|---|---|---|
| DA neurons | 22,048 | 8 / 5 / 4 | 1,464 / **139** / 884 | 97% NURR-positive |
| **Non-DA neurons** | 91,479 | 8 / 6 / 4 | 3,011 / 6,864 / 2,324 | **95% NURR-positive** |
| Microglia | 31,811 | 8 / 6 / 4 | 2,093 / 1,351 / 1,579 | 76% negative |
| Astrocytes | 33,506 | 8 / 6 / 4 | 1,750 / 1,931 / 1,570 | 89% negative |
| Oligodendrocytes | 178,815 | 8 / 6 / 4 | 7,630 / 10,446 / 8,030 | 90% negative |
| OPCs | 13,691 | 8 / 6 / 4 | 756 / 751 / 684 | 92% negative |

*Two problems this table exposes, and how S-6 handles them (fixed now):*

1. **Depth.** In PD donors, every comparator class has 5–75 times more nuclei than DA neurons. DA neurons are therefore the class most exposed to sparsity artefacts, so an Axis A reduction confined to DA neurons could partly reflect sparsity alone. **All comparator classes are downsampled per donor to the same T nuclei used in S-1**, so every class is compared at matched depth.
2. **Sorted fraction.** Glial classes come mainly from the NURR-negative fraction; DA neurons from the positive fraction. As §4.2 states, effect *sizes* cannot be compared across fractions. **Non-DA neurons, however, come from the same NURR-positive fraction as DA neurons.** They are therefore the **primary comparator**: same fraction, same broad cell type, depth-matched. Glial classes are compared on **direction only**.

*Procedure.* The identical Axis A model in each class, at matched depth, plotted side by side with confidence intervals.

| Pattern | Reported as |
|---|---|
| Down in DA neurons; flat in non-DA neurons; glial directions mixed or flat | **DA-selective**; supports H2 as framed |
| Down in DA neurons; smaller reduction in non-DA neurons | **DA-predominant**; gradient reported (magnitude comparison valid only against non-DA neurons) |
| Down comparably in DA **and** non-DA neurons | **Neuronal, not DA-specific**; described as such |
| Down in the same direction across most classes, glial included | **Global**; interferon-capacity interpretation withdrawn |

This is the single check most likely to change how a positive H2 is described.

### 6.7 How H2 is reported — the robustness ladder

| Level | Requirement | Permitted conclusion |
|---|---|---|
| **Statistical** | Holm-adjusted Kamath test rejects with β < 0 | Association detected in Kamath |
| **Technical** | Detectability parity acceptable; S-1 passes; S-3 does not clearly reverse | Not readily explained by differential sparsity |
| **Biological** | S-2 passes; A1 negative; S-6 shows DA selectivity or predominance | Consistent with reduced DA reception/transduction machinery |
| **External** | Same direction in an independent DA-capable cohort or DA-specific assay | Externally supported |

The strongest this design can reach is the biological level. **Even a result meeting all three internal levels is labelled UNREPLICATED** until the external level is met.

### 6.8 RQ2 pattern table (exploratory; gated)

Effect estimates and confidence intervals for A, A1, A3, B, D, and the mitophagy–fission and lysosomal sets in DA neurons, tabulated against P1–P3 (§3.3). Supporting: Spearman correlation of donor-level Axis A with mitophagy–fission score in the same cells, with confidence interval, uncorrected.

### 6.9 Combining the H1 cohorts

Cohorts are **never merged** into one dataset. The identical pipeline runs in each, producing β and SE per cohort, which are then combined.

**Primary pooled test — common-effect inverse-variance** *(changed in v2.0)*.

With three cohorts, between-cohort heterogeneity cannot be estimated reliably. The Hartung–Knapp random-effects method correctly reflects that uncertainty but has very low power (Bender 2018): with k = 3 it uses a t-distribution with 2 degrees of freedom, whose critical value is 4.30. Used as the primary test, it would probably fail to detect a real and consistent effect.

The common-effect model answers a narrower question — *is there an average effect across these cohorts?* That is a conditional inference about these cohorts, not a claim about all possible cohorts. **Generalisation rests on the concordance criterion below, not on the pooled p-value.**

**Mandatory sensitivity analyses:**
- REML random effects with Hartung–Knapp adjustment, always reported alongside.
- Hedges' g calculated from the donor-level pathway score with conventional pooled within-group standardisation.
- τ² and I² are reported but not interpreted: with k ≤ 3 they are close to uninformative.

**Per-cohort estimates are always shown**, so a pooled result driven by one cohort is visible.

**H1 is considered supported only if all three conditions hold:**

> (i) the point estimate is positive in **at least two** independent cohorts;
> (ii) **no** cohort has a 95% confidence interval lying entirely below zero; and
> (iii) the common-effect pooled p-value passes its Holm-adjusted threshold.

If (iii) passes but the random-effects estimate does not, H1 is reported as *supported across these cohorts, with limited evidence about generalisability*. Disagreement between the β/SE and Hedges' g syntheses is reported, never resolved by choosing the favourable scale.

H2 and H3 are single-cohort and are reported as **unreplicated**, in those words, in the abstract and slides, regardless of their Kamath results.

---

## 7. Prior forecasts

Personal forecasts recorded for later calibration. They play no role in the frequentist analysis.

| Outcome | P(supported) |
|---|---|
| H1 microglial interferon-responsive transcription ↑ | ~0.60 |
| H2 neuronal reception machinery ↓ | ~0.40 |
| H3 conjunction | ~0.25 |
| H2 survives S-1, given H2 | ~0.60 |
| Gate 0 fails for Axis B in DA neurons | ~0.25 |
| RQ2 evaluable at all | ~0.40 |

These are unchanged from v1.9. The v2.0 literature corrections concern mechanistic causality in animal models; they do not directly bear on the transcriptional predictions. The change to a common-effect primary test increases H1's power relative to v1.9, which would argue for a modest upward revision; the forecast is deliberately left unchanged rather than adjusted in the favourable direction.

A forecast-versus-outcome comparison appears in the final write-up.

---

## 8. Declared limitations

Written before results are known.

1. **Small donor numbers.** H2 rests on 8 controls and 5 PD donors; fewer for some analyses. No formal power calculation.
2. **Survivorship bias.** Every DA neuron in a PD brain is one that survived. Survivors may retain protective capacity, mount compensatory responses, or differ in subtype. The direction of this bias cannot be identified. H2 concerns surviving end-stage DA nuclei, not the cells already lost.
3. **Sorting confound.** Neurons and microglia come from different sorted fractions; directions are comparable, magnitudes are not.
4. **Differential sparsity.** ~10-fold difference in recoverable DA nuclei; addressed by S-1 but not eliminated.
5. **Subtype composition shift.** Addressed by S-2 but not eliminated.
6. **Nuclear RNA.** Under-represents cytoplasmic transcripts; low-abundance genes may be undetectable regardless of biology.
7. **Small gene sets.** Axis A has 7 genes, A1 4, A3 2. Set-level estimates are sensitive to single genes; per-gene plots are shown for this reason.
8. **Post-mortem, cross-sectional, end-stage tissue.** No temporal ordering; agonal state and PMI effects are only partly controlled.
9. **Inherited annotations.** Cell labels come from source publications. This improves reproducibility but excludes Wang despite its size.
10. **Transcript abundance is not signalling activity.** Protein, localisation and phosphorylation are not measured.
11. **Observational.** No perturbation; no causal claim will be made.
12. **Overlap with prior work.** Quan 2024 analysed Tier 2 for IFN-I; Villanueva 2026 displayed five of seven Axis A genes in Tier 1.
13. **Contradictory prior evidence.** Agarwal 2020 reported no association between PD genetic risk and microglia; this is addressed directly.
14. **Concentrated evidence base.** Each arm of the motivating model rests largely on one laboratory, and the STING-specific evidence for the damaging arm is contested (one retraction, one failed reproduction). This study tests a transcriptional prediction of the model; it cannot rescue or refute the mechanistic claims.
15. **Macrophage harmonisation.** Other cohorts' microglial labels may include border-associated macrophages that Kamath separates.
16. **No external H2 replication.** H2 and H3 are single-cohort results regardless of significance or internal robustness.

---

## 9. What would change the conclusion

| Observation | Consequence |
|---|---|
| Empirical null grossly non-uniform | Pipeline fault; halt and debug |
| Axis A fails Gate 0 in DA neurons | H2 **not evaluable** — never reported as null |
| Axis B fails Gate 0 in DA neurons | Neuronal responsive-transcription arm abandoned, not reframed |
| Detectability parity check positive | H2 reported as confounded by detection, regardless of its p-value |
| H2 fails S-1 | Reported as unresolved from sparsity |
| H2 fails S-2 | Reported as attributable to subtype loss rather than per-cell regulation |
| S-6 shows comparable reduction across most cell types | Interferon-capacity interpretation withdrawn; reported as a global finding |
| A negative but A1 not | Described as inducible-transducer-driven |
| Axis E moves comparably to Axis B | Reported as broad inflammatory signalling, not Type-I-selective |
| Same-direction effects in both compartments | H3 not supported |
| Any H1 cohort with 95% CI entirely below zero | H1 not supported, regardless of the pooled p-value |
| H2 null | RQ2 not evaluable |
| Always | H2 and H3 labelled unreplicated in the abstract |

---

## 10. Deliverables

### Before the October meeting (two-week pilot)
1. Frozen plan, `audit/donor_table.csv`, `harmonisation.md`, and executable Kamath scripts.
2. **Figure 1:** donor audit, Gate 0 detectability and Axis A detectability parity.
3. **Figure 2:** Kamath H1 and H2 estimates, with individual Axis A genes.
4. **Figure 3:** H2 robustness panel — S-1, S-2, A1 and S-6.
5. At least one completed independent H1 replication (Martirosyan preferred; Smajić if blocked).
6. One-page summary and five-slide presentation, marked *preliminary* wherever H1 synthesis is incomplete.

### Full project
1. Reproducible repository: numbered scripts `00_genesets` … `08_synthesis`, `DEVIATIONS.md`, `harmonisation.md`, gene-set files and a session-info lockfile.
2. H1 forest plot across all eligible cohorts, with common-effect, random-effects and Hedges' g syntheses.
3. Complete empirical null, all mandatory sensitivity analyses, and any RQ2 analyses the gate permits.
4. Manuscript-ready figures, summary, presentation, and a proposed external H2 validation study.

**Proposed follow-up, in order of priority:** (i) independent transcriptomic replication of H2 in a DA-capable cohort; (ii) protein and spatial validation in substantia nigra (RNAscope; immunostaining for IFNAR1 and phosphorylated STAT1); (iii) functional work in neuronal and microglial models. A longer-horizon idea — whether IFN-pathway CpGs (*IFNAR1*, *STAT1/2*, *DNM1L*, *PGAM5*) show coordinated methylation differences in PPMI — is a point of contact with UCL Genomics' work, explicitly flagged for modality mismatch: PPMI methylation is blood-derived, while this phenotype is nigral and cell-type-specific.

---

## 11. Timeline

| Days | Work | Output |
|---|---|---|
| 1–2 | Complete Appendix A; freeze and commit | Frozen plan; audit table |
| 3–5 | Kamath pseudobulk; Gate 0; detectability parity | Figure 1 |
| 6–9 | Kamath H1 and H2; gene-level results; S-1, S-2, S-3, S-5, S-6 | Core effects; H2 robustness panel |
| 10–12 | Martirosyan H1 (Smajić if blocked after one focused session) | ≥1 independent H1 estimate |
| 13–14 | Figures, summary, slides, repository clean-up | Meeting package |

**After the meeting:** third H1 cohort and pooled synthesis; full empirical null; remaining sensitivity analyses; RQ2 if permitted; approach Wang authors for annotations; manuscript only after result wording is fixed by the robustness ladder.

**Compute note.** S-1 must be complete before H2 is presented. The full 1,000-set empirical null may follow the meeting, but must first be tested on one cell class, and no confirmatory conclusion is final until it passes.

---

## 12. Amendments

Any change after freezing is recorded in `DEVIATIONS.md` with the date, the change, the reason, and whether it was made before or after the relevant result was seen. Post-hoc changes are permitted but must be labelled and reported in the write-up.

---

## 13. References

*Verified against primary sources in the v2.0 audit unless marked †. Every DOI should be checked by the author before submission.*

**Datasets**
1. Kamath T, *et al.* Single-cell genomic profiling of human dopamine neurons identifies a population that selectively degenerates in Parkinson's disease. *Nat Neurosci* 2022;25:588–595.
2. Smajić S, *et al.* Single-cell sequencing of human midbrain reveals glial activation and a Parkinson-specific neuronal state. *Brain* 2022;145:964–978.
3. Martirosyan A, *et al.* Unravelling cell type-specific responses to Parkinson's disease at single cell resolution. *Mol Neurodegener* 2024;19:7.
4. Wang Q, *et al.* Molecular profiling of human substantia nigra identifies diverse neuron types associated with vulnerability in Parkinson's disease. *Sci Adv* 2024;10:eadi8287.
5. Agarwal D, *et al.* A single-cell atlas of the human substantia nigra reveals cell-specific pathways associated with neurological disorders. *Nat Commun* 2020;11:4183.

**Glial / damaging arm**
6. Main BS, *et al.* Type-1 interferons contribute to the neuroinflammatory response and disease progression of the MPTP mouse model of Parkinson's disease. *Glia* 2016;64:1590–1604.
7. Main BS, *et al.* Type-1 interferons mediate the neuroinflammatory response and neurotoxicity induced by rotenone. *J Neurochem* 2017. †
8. Nazmi A, *et al.* Chronic neurodegeneration induces type I interferon synthesis via STING, shaping microglial phenotype and accelerating disease progression. *Glia* 2019;67:1254–1276.
9. Hinkle JT, Patel J, Panicker N, *et al.* STING mediates neurodegeneration and neuroinflammation in nigrostriatal α-synucleinopathy. *Proc Natl Acad Sci USA* 2022;119:e2118819119.
10. Klæstrup IH, Reinert LS, Ferreira SA, *et al.* Lack of functional STING modulates immunity but does not protect dopaminergic neurons in the alpha-synuclein pre-formed fibrils Parkinson's disease mouse model. *npj Parkinsons Dis* 2026;12:17. doi:10.1038/s41531-025-01228-0
11. Sliter DA, *et al.* Parkin and PINK1 mitigate STING-induced inflammation. *Nature* 2018;561:258–262. **RETRACTED** — Retraction Note: *Nature* 2025;644:1116. doi:10.1038/s41586-025-09346-8. *Cited for transparency only; carries no evidentiary weight.*
12. Quan P, Li X, Si Y, *et al.* Single cell analysis reveals the roles and regulatory mechanisms of type-I interferons in Parkinson's disease. *Cell Commun Signal* 2024;22:212. doi:10.1186/s12964-024-01590-1
13. Chen S, Crack PJ, Taylor JM. The contribution of type-I IFN-mediated neuroinflammation to Parkinson's disease progression. *Brain Behav Immun Health* 2025;101017. †

**Neuronal / protective arm**
14. Ejlerskov P, *et al.* Lack of neuronal IFN-β-IFNAR causes Lewy body- and Parkinson's disease-like dementia. *Cell* 2015;163:324–339.
15. Tresse E, *et al.* IFN-β rescues neurodegeneration by regulating mitochondrial fission via STAT5, PGAM5, and Drp1. *EMBO J* 2021;40:e106868.
16. Magalhaes J, *et al.* PIAS2-mediated blockade of IFN-β signaling: a basis for sporadic Parkinson disease dementia. *Mol Psychiatry* 2021. †
17. Villanueva EB, *et al.* Distinct and combined interferon-α/β-receptor-1 loss in neurons and astrocytes disrupt brain energy metabolism and drive Parkinsonian dementia. *J Biomed Sci* 2026;33. doi:10.1186/s12929-026-01257-8

**Other biology**
18. Borsche M, *et al.* Mitochondrial damage-associated inflammation highlights biomarkers in PRKN/PINK1 parkinsonism. *Brain* 2020;143:3041–3051.

**Methods**
19. Ritchie ME, *et al.* limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 2015;43:e47.
20. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 2012;40:e133.
21. Squair JW, *et al.* Confronting false discoveries in single-cell differential expression. *Nat Commun* 2021;12:5692.
22. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw* 2010;36:1–48.
23. IntHout J, Ioannidis JPA, Borm GF. The Hartung–Knapp–Sidik–Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian–Laird method. *BMC Med Res Methodol* 2014;14:25.
24. Bender R, Friede T, Koch A, *et al.* Methods for evidence synthesis in the case of very few studies. *Res Synth Methods* 2018;9:382–392. doi:10.1002/jrsm.1297

---

## Appendix A — Pre-freeze checklist

**All pre-freeze checks are complete.** None requires inspecting expression values.

| # | Check | Why it blocks | How |
|---|---|---|---|

**Completed:**
- ✅ Kamath donor counts at three thresholds (DA and microglia)
- ✅ Kamath H2 primary donor set verified for age, sex and PMI
- ✅ Kamath microglial-state feasibility (E2)
- ✅ Smajić donor split and covariates from GEO metadata
- ✅ Martirosyan donor split, covariates, DA yield and microglial yield
- ✅ Martirosyan nucleus-level annotations deposited; counts match Supplementary Table 2 exactly (C1)
- ✅ Smajić nucleus-level annotations deposited; all 11 donors ≥109 microglia (C1, C2)
- ✅ Gene-set membership, axis overlaps and set names verified; Axis E redefined (C3)
- ✅ Kamath comparator classes for S-6: all 18 donors pass; sorted fractions and depth imbalance recorded; S-6 redesigned (C4)
- ✅ SHA-256 checksums for every metadata and annotation file in `audit/sources.csv`; donor table regenerated by `audit/00_audit.py` and matching every count in §4 (C5)

*The pre-freeze audit script is written in Python because it only reads metadata; the analysis pipeline itself is in R (§6.4).*
- ✅ Wang cohort composition and annotation status
- ✅ Hallmark IFN-α / IFN-γ overlap (msigdbr 26.1.1)

---

## Appendix B — Version history

| Version | Principal changes (all pre-freeze; no expression data seen) |
|---|---|
| 1.0–1.1 | Initial hypotheses; RQ2 triage; Axis E added |
| 1.2 | Background reorganised by compartment manipulated |
| 1.3 | Third cohort promoted; H2 split |
| 1.4 | Statistical reconstruction: donor-level score regression replaced a non-poolable `camera`-based design; empirical null; S-1, S-2 |
| 1.5 | Full-text reading of Villanueva 2026: LBD-predominance, multi-cell-type effect, H2-LBD, Axis A3, S-6 |
| 1.6 | Thirteen internal corrections (Axis E contradiction, Gate 0 failure rules, LBD handling, S-1 specification) |
| 1.7 | Third-cohort audit: no DA replication cohort; H2 and H3 unreplicated; Wang excluded; Sliter 2018 retraction identified |
| 1.8 | Audit placeholders filled from metadata; Smajić split resolved; prose-vs-metadata observation |
| 1.9 | Gene-set membership computed; H1 wording constrained; Axis A-γ added |
| **2.0** | **Line-by-line audit against primary sources: five citation errors corrected; one failed replication and two single-laboratory evidence bases disclosed; H3 error control, PMI promotion, standardisation set and synthesis method corrected; four design gaps closed (§0)** |
