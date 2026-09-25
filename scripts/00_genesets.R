#!/usr/bin/env Rscript
# =============================================================================
# 00_genesets.R  (plan v2.0)
# Generates every gene set in the pre-analysis plan and runs pre-freeze
# check C3. Writes genesets/membership.csv, genesets/provenance.txt,
# genesets/core20_status.csv, genesets/axis_overlaps.csv.
# RUN: Rscript 00_genesets.R   (or Source in RStudio)
# =============================================================================
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(msigdbr)
dir.create("genesets", showWarnings = FALSE)

# ---- Hallmark (known to work without msigdbdf) ------------------------------
h <- msigdbr(species = "Homo sapiens", collection = "H")
gs <- function(df, nm) unique(df$gene_symbol[df$gs_name == nm])
a <- gs(h, "HALLMARK_INTERFERON_ALPHA_RESPONSE")
g <- gs(h, "HALLMARK_INTERFERON_GAMMA_RESPONSE")
oxphos <- gs(h, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
upr    <- gs(h, "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
stopifnot(length(a) > 0, length(g) > 0)

# ---- C2 sets: C3(c) confirm names exist -------------------------------------
# May require: install.packages("msigdbdf", repos="https://igordot.r-universe.dev")
c2 <- tryCatch(msigdbr(species = "Homo sapiens", collection = "C2"),
               error = function(e) { message("C2 unavailable: ", conditionMessage(e)); NULL })
find_set <- function(pattern) {
  if (is.null(c2)) return(NA_character_)
  hits <- unique(c2$gs_name[grepl(pattern, c2$gs_name)])
  if (length(hits) == 0) NA_character_ else hits[1]
}
reactome_name <- find_set("^REACTOME_INTERFERON_ALPHA_BETA_SIGNALING$")
kegg_name     <- find_set("^KEGG.*LYSOSOME$")
reactome <- if (is.na(reactome_name)) character(0) else gs(c2, reactome_name)
kegg_lys <- if (is.na(kegg_name))     character(0) else gs(c2, kegg_name)

# ---- Curated axes (NOT from MSigDB) -----------------------------------------
A   <- c("IFNAR1","IFNAR2","JAK1","TYK2","STAT1","STAT2","IRF9")
A1  <- c("IFNAR1","IFNAR2","JAK1","TYK2")
A2  <- c("STAT1","STAT2","IRF9")
A3  <- c("STAT5A","STAT5B")
Ag  <- c("IFNGR1","IFNGR2","JAK2")
C1  <- c("CGAS","STING1","TBK1","IRF3","IKBKE","TREX1")
C2s <- c("RIGI","IFIH1","MAVS")
D   <- c("SOCS1","SOCS3","USP18","PTPN2","PIAS2")
core20 <- c("ISG15","IFI6","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1",
            "IFIT2","IFIT3","IFI44","IFI44L","BST2","XAF1","RSAD2","IFITM3",
            "HERC5","EIF2AK2","SAMD9L")
ligands <- c("IFNB1","IFNE")
mito <- c("DNM1L","PGAM5","PINK1","PRKN","OPTN","INF2","MFN1","MFN2","OPA1",
          "FIS1","BNIP3","BNIP3L","SQSTM1","CALCOCO2")
lyso_extra <- c("GBA1","LAMP1","TFEB","CTSD","CTSB")
housekeep  <- c("ACTB","GAPDH","TBP","RPL13A","PPIA")
stopifnot(length(core20) == 20)

alpha_only <- setdiff(a, g); gamma_only <- setdiff(g, a); shared <- intersect(a, g)
# Plan v2.0 §5.6: Axis E excludes canonical Type-I genes that sit in the IFN-gamma set
axisE      <- setdiff(gamma_only, core20)
# Plan v2.0 §5.6b: Axis B without genes shared with curated axes (within-cell comparisons)
axisB_excl <- setdiff(a, c(A, D, C2s))

# ---- C3(a): classify each core-20 gene --------------------------------------
core_status <- data.frame(
  gene = core20,
  in_hallmark_alpha = core20 %in% a,
  in_hallmark_gamma = core20 %in% g,
  status = ifelse(core20 %in% alpha_only, "alpha-only",
           ifelse(core20 %in% shared, "shared alpha+gamma",
           ifelse(core20 %in% g, "gamma-only", "absent from Hallmark alpha"))))
write.csv(core_status, "genesets/core20_status.csv", row.names = FALSE)

# ---- C3(b): overlaps between curated axes and Axis B ------------------------
axes <- list(A = A, A1 = A1, A2 = A2, A3 = A3, A_gamma = Ag, C1 = C1,
             C2 = C2s, D = D, mito = mito)
ov <- do.call(rbind, lapply(names(axes), function(n) data.frame(
  axis = n,
  overlap_with_AxisB_hallmark_alpha = paste(intersect(axes[[n]], a), collapse = ";"),
  overlap_with_AxisE = paste(intersect(axes[[n]], gamma_only), collapse = ";"))))
write.csv(ov, "genesets/axis_overlaps.csv", row.names = FALSE)

# ---- Membership file ---------------------------------------------------------
b <- function(nm, v) if (length(v)) data.frame(set = nm, gene = sort(unique(v))) else NULL
out <- rbind(
  b("AxisA_reception", A), b("AxisA1_receptor_proximal", A1),
  b("AxisA2_ISGF3", A2), b("AxisA3_STAT5", A3), b("AxisAgamma_typeII_receptor", Ag),
  b("AxisB_hallmark_alpha", a), b("AxisB_reactome_sens", reactome),
  b("AxisB_core20_sens", core20), b("alpha_only_descriptive", alpha_only),
  b("AxisC1_DNA_sensing", C1), b("AxisC2_RNA_sensing", C2s),
  b("AxisD_negative_regulators", D), b("AxisD_minus_USP18", setdiff(D, "USP18")),
  b("AxisE_gamma_minus_alpha_minus_core20", axisE), b("gamma_minus_alpha_unfiltered", gamma_only),
  b("AxisB_excl_curated_overlaps", axisB_excl),
  b("shared_alpha_gamma", shared), b("ligands_descriptive", ligands),
  b("mitophagy_fission", mito), b("lysosomal", c(kegg_lys, lyso_extra)),
  b("oxphos", oxphos), b("UPR", upr), b("housekeeping_QC", housekeep))
write.csv(out, "genesets/membership.csv", row.names = FALSE)

# ---- Alias table (legacy symbols accepted at matching time) -----------------
write.csv(data.frame(
  current = c("CGAS","STING1","RIGI","PRKN","GBA1","TENT5A"),
  legacy  = c("MB21D1","TMEM173","DDX58","PARK2","GBA","FAM46A")),
  "genesets/aliases.csv", row.names = FALSE)

# ---- Provenance --------------------------------------------------------------
writeLines(c(
  paste("generated:", Sys.time()), paste("R:", R.version.string),
  paste("msigdbr:", as.character(packageVersion("msigdbr"))),
  paste("alpha:", length(a), "| gamma:", length(g), "| shared:", length(shared),
        "| alpha_only:", length(alpha_only), "| gamma_minus_alpha:", length(gamma_only),
        "| AxisE (minus core20):", length(axisE), "| AxisB_excl:", length(axisB_excl)),
  paste("REACTOME set:", ifelse(is.na(reactome_name), "NOT FOUND", reactome_name),
        "(n =", length(reactome), ")"),
  paste("KEGG lysosome set:", ifelse(is.na(kegg_name), "NOT FOUND", kegg_name),
        "(n =", length(kegg_lys), ")"),
  "Curated axes (A, A1-3, A-gamma, C1, C2, D, mito, lyso extras) are hand-specified here."
), "genesets/provenance.txt")

# ---- Console report for the plan --------------------------------------------
cat("\n==== C3 REPORT ====\n")
cat("alpha", length(a), "| gamma", length(g), "| shared", length(shared),
    "| alpha-only", length(alpha_only), "| gamma-minus-alpha", length(gamma_only),
    "| Axis E (final)", length(axisE), "| Axis B-excl", length(axisB_excl), "\n\n")
cat("(a) Core-20 status:\n"); print(table(core_status$status))
print(core_status[core_status$status != "shared alpha+gamma", c("gene","status")], row.names = FALSE)
cat("\n(b) Curated axes overlapping Axis B (Hallmark alpha):\n")
print(ov[ov$overlap_with_AxisB_hallmark_alpha != "", 1:2], row.names = FALSE)
cat("\n(c) REACTOME:", ifelse(is.na(reactome_name), "NOT FOUND", reactome_name),
    "| KEGG lysosome:", ifelse(is.na(kegg_name), "NOT FOUND", kegg_name), "\n")
cat("\nFiles written to genesets/. Send the C3 REPORT above.\n")
