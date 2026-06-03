#DDS --> to Union table
# Differential expression integration across 4h + 24h IL-4
# (pulse vs continuous) with OR-gate DEG selection, UNION,
# and clean gene-symbol filtering. Saves:
#   1) Master expressed-gene log2FC+padj table (all genes kept)
#   2) Pulse-significant DEG table (OR-gate)
#   3) Continuous-significant DEG table (OR-gate)
#   4) UNION table (pulse wins ties) with uninformative symbols removed
#   5) Uncapped log2FC matrix used for downstream clustering (optional later)
#   6) Removed-uninformative symbol list

# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
})

# ----------------------------
# 0) User-defined paths
# ----------------------------

# Base analysis directory (all outputs will be saved under this folder)
base_dir <- "/Users/iyad-si/Desktop/KMJ lab/Mccleary data analysis/07_15_2025_4h_24h_bulkrnaseq_cVSp_ISI_analysis/07_16_2025_DEG analysis CVSP 4h_24h_try3"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# DESeq2 objects (RDS)
dds_4hrs_rds  <- file.path(base_dir, "../../07_15_2025_4hours_cvsp_bulkrnaseq_isi_analysis_batch_corrected_vstnorm/dds_4hrs.rds")
dds_24hrs_rds <- file.path(base_dir, "../../07_15_2025_24hours_cvsp_bulkrnaseq_isi_analysis_batch_corrected_vstnorm/dds_24hrs.rds")

if (!file.exists(dds_4hrs_rds))  stop("Missing file: ", dds_4hrs_rds)
if (!file.exists(dds_24hrs_rds)) stop("Missing file: ", dds_24hrs_rds)

# Output folders
or_gate_dir <- file.path(base_dir, "OR gate filtering And hierarchial clustering")
dir.create(or_gate_dir, recursive = TRUE, showWarnings = FALSE)

final_out_dir <- file.path(or_gate_dir, "Union_and_clean_symbols__publication_ready__2026_03_02")
dir.create(final_out_dir, recursive = TRUE, showWarnings = FALSE)

# Output files
master_outfile <- file.path(
  final_out_dir,
  "log2fc_padj_table__c100_vs_p100__4h_24h__expressed_genes_only.csv"
)

pulse_sig_file <- file.path(
  final_out_dir,
  "log2fc_padj_table__c100_vs_p100__4h_24h__pulse_significant_genes.csv"
)

cont_sig_file <- file.path(
  final_out_dir,
  "log2fc_padj_table__c100_vs_p100__4h_24h__continuous_significant_genes.csv"
)

removed_symbols_file <- file.path(
  final_out_dir,
  "removed_uninformative_gene_symbols__optionA.csv"
)

union_clean_wide_file <- file.path(
  final_out_dir,
  "log2fc_and_padj_table__union_pulse_wins__clean_symbols__optionA.csv"
)

log2fc_matrix_uncapped_file <- file.path(
  final_out_dir,
  "log2fc_matrix_uncapped__union_pulse_wins__clean_symbols__optionA.csv"
)

# OR-gate thresholds (match your existing logic)
lfc_cutoff  <- 1
padj_cutoff <- 0.05

# ----------------------------
# 1) Load DESeq2 objects
# ----------------------------
cat("\n[1/6] Loading DESeq2 objects...\n")
dds_4hrs  <- readRDS(dds_4hrs_rds)
dds_24hrs <- readRDS(dds_24hrs_rds)

# ----------------------------
# 2) Expression filter (keeps genes expressed in at least one
#    stimulated condition across either timepoint)
#
#    Rule:
#      Keep gene if summed raw counts > 0 across:
#        - 4h c100 OR p100 samples, OR
#        - 24h C OR P samples
# ----------------------------
cat("[2/6] Applying expressed-gene filter...\n")
counts_4h  <- counts(dds_4hrs)
counts_24h <- counts(dds_24hrs)

# Identify relevant stimulated samples by name pattern
samples_4h  <- colnames(counts_4h)[grepl("c100|p100", colnames(counts_4h))]
samples_24h <- colnames(counts_24h)[grepl("^C\\d+|^P\\d+", colnames(counts_24h))]

if (length(samples_4h) == 0) {
  stop("No 4h stimulated samples found with grepl('c100|p100'). Check colnames(counts_4h).")
}
if (length(samples_24h) == 0) {
  stop("No 24h stimulated samples found with grepl('^C\\d+|^P\\d+'). Check colnames(counts_24h).")
}

genes_to_keep <- rowSums(counts_4h[, samples_4h, drop = FALSE]) > 0 |
                 rowSums(counts_24h[, samples_24h, drop = FALSE]) > 0

dds_4hrs  <- dds_4hrs[genes_to_keep, ]
dds_24hrs <- dds_24hrs[genes_to_keep, ]

cat("  Genes retained after expressed-gene filter:", sum(genes_to_keep), "\n")

# ----------------------------
# 3) Extract DESeq2 results for the four pre-specified contrasts
#
#    Contrasts (each vs untreated control 'ct'):
#      - 4h:  p100 vs ct, c100 vs ct
#      - 24h: p100 vs ct, c100 vs ct
#
#    Output columns:
#      ensembl_gene_id, log2fc_<contrast>, padj_<contrast>
# ----------------------------
cat("[3/6] Extracting DESeq2 contrasts...\n")

get_clean_result <- function(dds_obj, contrast, label) {
  # results() returns a DESeqResults object; convert to data.frame and keep only
  # log2FC + padj, with standardized column names.
  as.data.frame(results(dds_obj, contrast = contrast)) %>%
    rownames_to_column("ensembl_gene_id") %>%
    as_tibble() %>%
    dplyr::select(ensembl_gene_id, log2FoldChange, padj) %>%
    dplyr::rename(
      !!paste0("log2fc_", label) := log2FoldChange,
      !!paste0("padj_",  label) := padj
    )
}

# IMPORTANT:
# This assumes your DESeq2 design uses `condition` and that levels include:
#   ct, p100, c100
res_p100_4h  <- get_clean_result(dds_4hrs,  c("condition", "p100", "ct"), "p100_4h")
res_c100_4h  <- get_clean_result(dds_4hrs,  c("condition", "c100", "ct"), "c100_4h")
res_p100_24h <- get_clean_result(dds_24hrs, c("condition", "p100", "ct"), "p100_24h")
res_c100_24h <- get_clean_result(dds_24hrs, c("condition", "c100", "ct"), "c100_24h")

# Merge all contrasts into one master table (outer joins preserve all genes)
all_res <- list(res_p100_4h, res_c100_4h, res_p100_24h, res_c100_24h) %>%
  purrr::reduce(full_join, by = "ensembl_gene_id")

# ----------------------------
# 4) Map Ensembl gene IDs to mouse gene symbols (MGI)
#    using biomaRt (Ensembl mmusculus dataset)
# ----------------------------
cat("[4/6] Mapping Ensembl IDs to MGI gene symbols...\n")

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(all_res$ensembl_gene_id),
  mart       = ensembl
) %>%
  as_tibble() %>%
  rename(gene_symbol = mgi_symbol) %>%
  mutate(across(everything(), as.character))

final_df <- all_res %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  relocate(gene_symbol, .after = ensembl_gene_id)

# Save the master expressed-gene table (all genes after expression filter)
write_csv(final_df, master_outfile)
cat("  Wrote master table:\n  ", master_outfile, "\n", sep = "")

# ----------------------------
# 5) OR-gate DEG selection (pulse significant vs continuous significant)
#
#    Pulse-significant gene if:
#      (|log2FC_p100_4h|  > lfc_cutoff AND padj_p100_4h  < padj_cutoff) OR
#      (|log2FC_p100_24h| > lfc_cutoff AND padj_p100_24h < padj_cutoff)
#
#    Continuous-significant gene if:
#      (|log2FC_c100_4h|  > lfc_cutoff AND padj_c100_4h  < padj_cutoff) OR
#      (|log2FC_c100_24h| > lfc_cutoff AND padj_c100_24h < padj_cutoff)
# ----------------------------
cat("[5/6] Applying OR-gate DEG filters...\n")

pulse_filtered <- final_df %>%
  filter(
    (abs(log2fc_p100_4h)  > lfc_cutoff & padj_p100_4h  < padj_cutoff) |
      (abs(log2fc_p100_24h) > lfc_cutoff & padj_p100_24h < padj_cutoff)
  )

cont_filtered <- final_df %>%
  filter(
    (abs(log2fc_c100_4h)  > lfc_cutoff & padj_c100_4h  < padj_cutoff) |
      (abs(log2fc_c100_24h) > lfc_cutoff & padj_c100_24h < padj_cutoff)
  )

write_csv(pulse_filtered, pulse_sig_file)
write_csv(cont_filtered,  cont_sig_file)

cat("  Pulse-significant genes:", nrow(pulse_filtered), "\n")
cat("  Continuous-significant genes:", nrow(cont_filtered), "\n")

# ----------------------------
# 6) UNION (pulse wins ties) + clean symbol filtering (Option A)
#
#    UNION:
#      bind_rows(pulse, continuous) then distinct(ensembl_gene_id)
#      => If a gene appears in both, the pulse row is retained (wins ties).
#
#    Option A removes uninformative symbols:
#      - NA or ""
#      - symbols that begin with ENSMUSG
#      - symbols that begin with Gm
#      - symbols ending in Rik / Rik2 / Riky / ... (Rik + trailing alphanumerics)
#      - symbols that begin with LOC
#
#    Outputs:
#      - removed gene list (transparency)
#      - union-wide table (log2FC + padj)
#      - uncapped log2FC matrix CSV (P4, P24, C4, C24)
# ----------------------------
cat("[6/6] UNION + clean-symbol filtering + outputs...\n")

pulse_df <- readr::read_csv(pulse_sig_file, show_col_types = FALSE)
cont_df  <- readr::read_csv(cont_sig_file,  show_col_types = FALSE)

combined_df <- bind_rows(pulse_df, cont_df) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(gene_symbol = as.character(gene_symbol))

is_unhelpful_symbol <- function(sym) {
  is.na(sym) |
    sym == "" |
    str_starts(sym, "ENSMUSG") |
    str_starts(sym, "Gm") |
    str_detect(sym, "Rik[[:alnum:]]*$") |
    str_starts(sym, "LOC")
}

removed_genes <- combined_df %>% filter(is_unhelpful_symbol(gene_symbol))
write_csv(removed_genes, removed_symbols_file)

combined_df_clean <- combined_df %>% filter(!is_unhelpful_symbol(gene_symbol))

# Wide table including both log2FC + padj, with standardized names
wide_df <- combined_df_clean %>%
  transmute(
    ensembl_gene_id,
    gene_symbol,

    log2fc_p100_4h  = log2fc_p100_4h,
    padj_p100_4h    = padj_p100_4h,

    log2fc_p100_24h = log2fc_p100_24h,
    padj_p100_24h   = padj_p100_24h,

    log2fc_c100_4h  = log2fc_c100_4h,
    padj_c100_4h    = padj_c100_4h,

    log2fc_c100_24h = log2fc_c100_24h,
    padj_c100_24h   = padj_c100_24h
  ) %>%
  mutate(
    # Replace missing log2FC with 0 for matrix construction only.
    # NOTE: padj remains NA if missing (do not coerce unless you want padj=1).
    across(starts_with("log2fc_"), ~replace_na(., 0)),
    label = make.unique(gene_symbol)
  )

# Save the UNION wide table (no "label" column in the saved file)
write_csv(wide_df %>% select(-label), union_clean_wide_file)

# Build & save uncapped log2FC matrix for downstream clustering (optional)
log2fc_matrix_uncapped <- wide_df %>%
  select(log2fc_p100_4h, log2fc_p100_24h, log2fc_c100_4h, log2fc_c100_24h) %>%
  as.matrix()
rownames(log2fc_matrix_uncapped) <- wide_df$label

write.csv(
  log2fc_matrix_uncapped,
  file = log2fc_matrix_uncapped_file,
  quote = FALSE,
  row.names = TRUE
)

# ----------------------------
# Final console summary
# ----------------------------
cat("\n=== Outputs written ===\n")
cat("Master expressed-gene table:\n  ", master_outfile, "\n", sep = "")
cat("Pulse-significant table:\n  ", pulse_sig_file, "\n", sep = "")
cat("Continuous-significant table:\n  ", cont_sig_file, "\n", sep = "")
cat("Removed uninformative symbols:\n  ", removed_symbols_file, "\n", sep = "")
cat("UNION (pulse wins) clean-symbol wide table:\n  ", union_clean_wide_file, "\n", sep = "")
cat("Uncapped log2FC matrix (for clustering later):\n  ", log2fc_matrix_uncapped_file, "\n", sep = "")
cat("\nCounts:\n")
cat("  Pulse sig genes:", nrow(pulse_filtered), "\n")
cat("  Continuous sig genes:", nrow(cont_filtered), "\n")
cat("  Total after UNION:", nrow(combined_df), "\n")
cat("  Removed symbols:", nrow(removed_genes), "\n")
cat("  Kept for downstream clustering:", nrow(combined_df_clean), "\n")