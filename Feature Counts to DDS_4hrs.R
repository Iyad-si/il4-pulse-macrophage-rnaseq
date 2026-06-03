
# Script: run_deseq2_normalization_4hrs.R                              # File name
# Purpose: Normalize and batch-correct 4-hour RNA-seq data             # High-level goal
# Author: Iyad Sayed Issa                                              # Author metadata
# Date: July 2025                                                      # Date metadata
# Input: featureCounts output (raw counts)                              # Expected input format/source
# Output:                                                              # Outputs produced by the script
#   - dds_4hrs.rds: DESeq2 object                                      # Saved DESeq2 object for downstream DE
#   - vst_normalized_batchcorrected_4hrs_with_symbols.csv: normalized matrix with gene symbols  # Batch-corrected VST matrix
############################################################

#===========================                                            # Section divider
# Load required libraries                                               # Load all dependencies
#===========================
library(DESeq2)    # Differential expression analysis                  # Core package for DE and normalization
library(tibble)    # Tidy data frames                                   # Provides tibble data structure
library(dplyr)     # Data manipulation                                  # Verbs for data wrangling (mutate/select/etc.)
library(readr)     # Reading TSV files                                  # Fast, friendly read_tsv/read_csv
library(stringr)   # String handling                                    # Regex/string utilities for name cleanup
library(limma)     # Batch correction                                   # removeBatchEffect for expression matrices
library(biomaRt)   # Gene annotation via Ensembl                        # Map Ensembl IDs to MGI symbols

#===========================                                            # Section divider
# Step 1: Define file paths                                             # Point to inputs and where to save outputs
#===========================
counts_file <- "/Users/iyad-si/Desktop/KMJ lab/Mccleary data analysis/07_15_2025_4h_24h_bulkrnaseq_cVSp_ISI_analysis/gene_counts_4hours_CvsP copy.txt"  # featureCounts raw counts file
output_dir  <- "/Users/iyad-si/Desktop/KMJ lab/Mccleary data analysis/07_15_2025_4h_24h_bulkrnaseq_cVSp_ISI_analysis"                                 # Output directory

vst_csv  <- file.path(output_dir, "vst_normalized_batchcorrected_4hrs_with_symbols.csv")   # Path for VST+batch-corrected CSV
dds_rds  <- file.path(output_dir, "dds_4hrs.rds")                                          # Path for serialized DESeq2 object

#===========================                                            # Section divider
# Step 2: Load raw count matrix                                         # Read counts table and prepare matrix
#===========================
# Skip header lines starting with "#"                                   # featureCounts can prepend comment header lines
raw_counts <- read_tsv(counts_file, comment = "#")                      # Load counts; "#" lines ignored

# Sample columns start at column 7 (first 6 columns are gene info)      # Assume first 6 cols are gene metadata
sample_cols <- 7:ncol(raw_counts)                                       # Indices of sample count columns

# Clean sample names from full BAM paths                                # Normalize sample names to concise labels
# e.g., "/.../c_100_1_Aligned.sortedByCoord.out.bam" → "c100_1"         # Example transformation target
colnames(raw_counts)[sample_cols] <- colnames(raw_counts)[sample_cols] %>%  # Take current sample column names
  basename() %>%                                # keep only file name     # Strip directories from paths
  str_remove("_Aligned.*") %>%                  # remove "_Aligned.sortedByCoord.out.bam"  # Drop aligner suffix
  str_replace("^c_", "c") %>%                   # fix c_100_1 → c100_1    # Normalize continuous prefix
  str_replace("^p_", "p") %>%                   # fix p_10_1 → p10_1      # Normalize pulse prefix
  str_replace("^Ct_", "ct_")                    # fix Ct_1 → ct_1         # Normalize control prefix/case

# Extract gene-level count matrix                                       # Build numeric matrix with genes as rows
gene_ids <- raw_counts$Geneid                                           # Store Ensembl IDs as row names source
count_matrix <- raw_counts[, sample_cols] %>% as.data.frame()           # Subset to sample columns and coerce to data.frame
rownames(count_matrix) <- gene_ids                                      # Set rownames to Ensembl gene IDs

#===========================                                            # Section divider
# Step 3: Create sample metadata (colData)                              # Derive condition/batch for each sample
#===========================
sample_names <- colnames(count_matrix)                                   # Vector of sample names (columns)

# Assign IL-4 treatment condition based on sample name prefix           # Parse condition from standardized names
# Order matters: match "c100" before "c10", "p100" before "p10"         # Ensure more specific patterns match first
condition <- case_when(                                                  # Map name patterns to condition labels
  str_detect(sample_names, "^c100") ~ "c100",                            # Continuous 100 ng/mL
  str_detect(sample_names, "^c10")  ~ "c10",                             # Continuous 10 ng/mL
  str_detect(sample_names, "^p100") ~ "p100",                            # Pulse 100 ng/mL
  str_detect(sample_names, "^p10")  ~ "p10",                             # Pulse 10 ng/mL
  str_detect(sample_names, "^ct")   ~ "ct",                              # Control
  TRUE ~ NA_character_                                                  # Fallback if no pattern matched
)
condition <- factor(condition, levels = c("ct", "c10", "c100", "p10", "p100"))  # Encode as factor with ct as reference

# Assign sequencing batch: replicate 1 is batch1; all others batch2     # Derive batch from suffix
batch <- ifelse(str_detect(sample_names, "_1$"), "batch1", "batch2")    # _1 → batch1; others → batch2
batch <- factor(batch)                                                  # Store as factor

# Construct sample metadata table                                       # Combine into colData frame
coldata <- data.frame(                                                  # Create data.frame for DESeq2 colData
  row.names = sample_names,                                             # Row names must match columns of count_matrix
  condition = condition,                                                # Experimental condition
  batch = batch                                                         # Batch indicator
)

#===========================                                            # Section divider
# Step 4: Create DESeq2 object and normalize                            # Build DESeq2 dataset and run DESeq pipeline
#===========================
dds <- DESeqDataSetFromMatrix(                                          # Construct DESeq2 object from counts
  countData = count_matrix,                                             # Raw counts matrix (integers)
  colData = coldata,                                                    # Sample metadata with condition/batch
  design = ~ batch + condition                                          # Model: adjust for batch, test condition
)

# Run DESeq2 pipeline                                                   # Estimate size factors/dispersion and fit GLM
dds <- DESeq(dds)                                                       # Executes normalization and model fitting

#DDS_4hrs ends here!
# ===========================             

# Variance Stabilizing Transformation (VST)                             # Transform counts to homoscedastic scale
# Set blind = FALSE to retain experimental design structure             # Uses design for dispersion trend (not blind)
vst <- vst(dds, blind = FALSE)                                          # Produce VST-transformed assay

#===========================                                            # Section divider
# Step 5: Batch correction (removeBatchEffect)                          # Remove batch trends from VST for visualization
#===========================
vst_mat <- assay(vst)  # Extract VST matrix                             # Get numeric matrix from DESeqTransform
vst_mat_corrected <- removeBatchEffect(vst_mat,                         # Apply limma batch correction
                                       batch = coldata$batch)           # Specify batch covariate

#===========================                                            # Section divider
# Step 6: Add gene symbols using biomaRt                                # Annotate Ensembl IDs with MGI symbols
#===========================
# Connect to Ensembl and select mouse gene dataset                      # Initialize biomaRt connection
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")       # Choose Mus musculus dataset

# Retrieve gene symbol mapping for all Ensembl gene IDs                 # Query Ensembl for ID→symbol mapping
gene_map <- getBM(                                                      # Download mapping table
  attributes = c("ensembl_gene_id", "mgi_symbol"),                      # Columns to retrieve
  filters = "ensembl_gene_id",                                          # Filter type is Ensembl gene IDs
  values = rownames(vst_mat_corrected),                                 # The IDs present in our matrix
  mart = ensembl                                                        # Mart connection object
)

# Remove duplicated mappings to keep 1-to-1 relationships               # Ensure unique Ensembl→symbol rows
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]           # Drop duplicate Ensembl IDs

# Prepare VST matrix as data frame                                      # Convert corrected matrix for joining
vst_df <- as.data.frame(vst_mat_corrected)                              # Coerce to data.frame
vst_df$ensembl_gene_id <- rownames(vst_df)                              # Preserve Ensembl IDs as explicit column

# Join gene symbols and re-order columns                                # Merge symbols and move to front
vst_df <- left_join(vst_df, gene_map, by = "ensembl_gene_id")           # Left-join to keep all rows
vst_df <- vst_df %>% relocate(mgi_symbol, .before = ensembl_gene_id)    # Place symbol column first

#===========================                                            # Section divider
# Step 7: Export outputs                                                # Write final artifacts to disk
#===========================
# Save normalized matrix with gene symbols                              # Export batch-corrected VST with symbols
write.csv(vst_df, file = vst_csv, row.names = FALSE)                    # Write CSV (no rownames column)

# Save DESeq2 object for downstream DE analysis                         # Persist DESeq2 object for contrasts
saveRDS(dds, file = dds_rds)                                            # Save RDS for later reuse
                        # Echo RDS path
