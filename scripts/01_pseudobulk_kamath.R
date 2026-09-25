#!/usr/bin/env Rscript
# =============================================================================
# 01_pseudobulk_kamath.R
#
# Builds donor x cell-class pseudobulk count matrices from the Kamath SCP1768
# matrix WITHOUT loading it into memory.
#
# The .mtx file is a plain list of "row col value" triplets. We read it in
# chunks, keep only the nuclei we care about (DA neurons and microglia), and
# add each value into a small running total. Peak memory stays well under 2 GB.
#
# INPUT  (all in the same folder as this script):
#   Homo_matrix.mtx.gz      4.74 GB   the expression matrix
#   Homo_features.tsv.gz              gene names
#   Homo_bcd.tsv.gz                   nucleus barcodes (column order)
#   da_UMAP.tsv                       DA neuron annotations
#   mg_UMAP.tsv                       microglia annotations (incl. one macrophage cluster)
#   astro_UMAP.tsv, olig_UMAP.tsv,    OPTIONAL (plan check C4) - needed for S-6
#   opc_UMAP.tsv, nonda_UMAP.tsv      cell-type specificity. Skipped with a warning if absent.
#   METADATA_PD.tsv (.gz ok)          donor + diagnosis per nucleus
#
# OUTPUT (small - these are what you send back):
#   pseudobulk_counts.csv.gz          genes x (donor_cellclass) counts
#   pseudobulk_groups.csv             group definitions + nuclei counts
#   pseudobulk_log.txt                run log
#
# RUN:  Rscript 01_pseudobulk_kamath.R
# TIME: 20-60 minutes. Progress prints as it goes. Leave it running.
# =============================================================================

t_start <- Sys.time()
options(stringsAsFactors = FALSE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(msg, "\n"); cat(msg, "\n", file = "pseudobulk_log.txt", append = TRUE)
}

pick <- function(...) { for (f in c(...)) if (file.exists(f)) return(f); NA }

f_mtx  <- pick("Homo_matrix.mtx.gz", "Homo_matrix.mtx")
f_feat <- pick("Homo_features.tsv.gz", "Homo_features.tsv")
f_bcd  <- pick("Homo_bcd.tsv.gz", "Homo_bcd.tsv")
f_meta <- pick("METADATA_PD.tsv.gz", "METADATA_PD.tsv")
for (f in c(f_mtx, f_feat, f_bcd, f_meta, "da_UMAP.tsv", "mg_UMAP.tsv"))
  if (is.na(f) || !file.exists(f)) stop("Missing input file: ", f)

log_msg("=== Kamath pseudobulk build ===")

# ---- 1. Annotations -> which nuclei, and which group each belongs to --------
read_ann <- function(f) { x <- read.delim(f, check.names = FALSE); x[x$NAME != "TYPE", ] }

# Read only the 3 needed columns (the full file is ~470 MB uncompressed)
hdr  <- strsplit(readLines(gzfile(f_meta), n = 1), "\t")[[1]]
hdr  <- gsub('"', "", hdr)
need <- c("NAME", "donor_id", "Status")
cc   <- ifelse(hdr %in% need, "character", "NULL")
meta <- read.delim(gzfile(f_meta), check.names = FALSE, colClasses = cc, quote = '"')
meta <- meta[meta$NAME != "TYPE", ]
key  <- meta[, c("NAME", "donor_id", "Status")]

mk <- function(df, class) data.frame(NAME = df$NAME, donor = df$donor_id,
                                     status = df$Status, class = class)

da <- merge(read_ann("da_UMAP.tsv"), key, by = "NAME")
mg <- merge(read_ann("mg_UMAP.tsv"), key, by = "NAME")
da$family  <- sub("_.*$", "", da$Cell_Type)            # SOX6_AGTR1 -> SOX6
is_macro   <- grepl("^Macro", mg$Cell_Type)            # Macro_CD200R1

ann <- rbind(
  mk(da, "DA_all"),
  mk(da, paste0("DA_", da$family)),
  mk(mg[!is_macro, ], "MG_all"),                        # PRIMARY H1 class (plan 6.2)
  mk(mg, "MG_incl_macro"),                              # sensitivity
  mk(mg, paste0("MGstate_", mg$Cell_Type))              # E2
)

# Optional classes for S-6 (plan check C4)
optional <- c(ASTRO = "astro_UMAP.tsv", OLIG = "olig_UMAP.tsv",
              OPC = "opc_UMAP.tsv", NONDA = "nonda_UMAP.tsv")
for (cls in names(optional)) {
  f <- optional[[cls]]
  if (file.exists(f)) {
    x <- merge(read_ann(f), key, by = "NAME")
    ann <- rbind(ann, mk(x, paste0(cls, "_all")))
    log_msg("Loaded ", f, ": ", nrow(x), " nuclei")
  } else log_msg("WARNING: ", f, " not found - S-6 will lack ", cls)
}
ann$group <- paste(ann$donor, ann$class, sep = "|")

log_msg("Annotated nuclei: DA = ", nrow(da), ", MG = ", sum(!is_macro),
        ", macrophage = ", sum(is_macro))
log_msg("Groups defined: ", length(unique(ann$group)))

# ---- 2. Barcodes -> column indices -----------------------------------------
bcd <- readLines(gzfile(f_bcd)); bcd <- trimws(bcd)
bcd <- bcd[bcd != "" & bcd != "NAME"]
log_msg("Barcodes in matrix: ", length(bcd))

idx <- match(ann$NAME, bcd)
if (all(is.na(idx))) stop("No annotation barcodes matched the matrix barcode file.")
miss <- sum(is.na(idx))
if (miss > 0) log_msg("WARNING: ", miss, " annotated nuclei not found in barcodes; dropped.")
ann <- ann[!is.na(idx), ]; ann$col <- idx[!is.na(idx)]

groups   <- sort(unique(ann$group))
gidx     <- match(ann$group, groups)
n_groups <- length(groups)

# A nucleus belongs to SEVERAL groups at once (e.g. DA_all AND DA_SOX6; MG_all AND
# MG_incl_macro AND its MG state). Store every membership: one column per slot.
memb    <- split(gidx, ann$col)
memb_n  <- lengths(memb)
K       <- max(memb_n)
col2grp <- matrix(0L, nrow = length(bcd), ncol = K)
cols    <- as.integer(names(memb))
for (k in seq_len(K)) {
  has_k <- memb_n >= k
  col2grp[cols[has_k], k] <- vapply(memb[has_k], function(v) v[k], integer(1))
}
log_msg("Membership slots per nucleus: max ", K)

# ---- 3. Genes ---------------------------------------------------------------
feat <- read.delim(gzfile(f_feat), header = FALSE, check.names = FALSE)
genes <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
log_msg("Genes: ", length(genes))

# ---- 4. Stream the matrix ---------------------------------------------------
con <- gzfile(f_mtx, "r")
hdr <- readLines(con, n = 1)
if (!grepl("^%%MatrixMarket", hdr)) stop("Not a MatrixMarket file.")
repeat { ln <- readLines(con, n = 1); if (!grepl("^%", ln)) break }
dims <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
n_row <- dims[1]; n_col <- dims[2]; n_nz <- dims[3]
log_msg("Matrix: ", n_row, " x ", n_col, " with ", n_nz, " non-zero entries")

# Auto-detect orientation: genes may be rows or columns
if (n_row == length(genes) && n_col == length(bcd)) {
  genes_are_rows <- TRUE
} else if (n_col == length(genes) && n_row == length(bcd)) {
  genes_are_rows <- FALSE
  log_msg("NOTE: matrix is cells x genes; transposing on the fly.")
} else {
  stop("Dimensions do not match features/barcodes. Matrix ", n_row, "x", n_col,
       "; genes ", length(genes), "; barcodes ", length(bcd))
}

acc <- matrix(0, nrow = length(genes), ncol = n_groups)
CHUNK <- 5e6
done <- 0; kept <- 0

repeat {
  v <- scan(con, what = numeric(), n = CHUNK * 3, quiet = TRUE)
  if (length(v) == 0) break
  if (length(v) %% 3 != 0) stop("Malformed matrix: entries are not complete triplets.")
  m <- matrix(v, ncol = 3, byrow = TRUE)
  if (genes_are_rows) { g <- as.integer(m[, 1]); cl <- as.integer(m[, 2]) }
  else                { g <- as.integer(m[, 2]); cl <- as.integer(m[, 1]) }
  x <- m[, 3]

  kept <- kept + sum(col2grp[cl, 1] > 0L)
  for (k in seq_len(K)) {                      # add each entry to EVERY group
    grp  <- col2grp[cl, k]
    keep <- grp > 0L
    if (!any(keep)) next
    flat <- (grp[keep] - 1L) * length(genes) + g[keep]      # integer index
    add  <- tapply(x[keep], flat, sum)
    idx  <- as.integer(names(add))
    acc[idx] <- acc[idx] + add
  }
  done <- done + nrow(m)
  log_msg(sprintf("  %.1f%% (%d / %d entries; %d from annotated nuclei)",
                  100 * done / n_nz, done, n_nz, kept))
}
close(con)

# ---- 5. Write outputs -------------------------------------------------------
dimnames(acc) <- list(genes, groups)
acc <- acc[rowSums(acc) > 0, , drop = FALSE]     # drop all-zero genes
log_msg("Pseudobulk: ", nrow(acc), " genes x ", ncol(acc), " groups")

gz <- gzfile("pseudobulk_counts.csv.gz", "w")
write.csv(acc, gz); close(gz)

gi <- unique(ann[, c("group", "donor", "status", "class")])
nn <- as.data.frame(table(ann$group)); names(nn) <- c("group", "n_nuclei")
gi <- merge(gi, nn, by = "group")
write.csv(gi[order(gi$class, gi$status, gi$donor), ],
          "pseudobulk_groups.csv", row.names = FALSE)

log_msg("DONE in ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")
log_msg("Send back: pseudobulk_counts.csv.gz, pseudobulk_groups.csv, pseudobulk_log.txt")
