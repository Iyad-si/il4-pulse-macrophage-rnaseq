# Script: run_deseq2_normalization_24hrs.R                              # Script name
# Purpose: Normalize and annotate 24-hour IL-4 RNA-seq data             # Goal: normalization + gene annotation
# Input: featureCounts output with 9 samples                            # Input file type and sample count
# Output:                                                               # Outputs produced by this script
#   - dds_24hrs.rds: DESeq2 object                                      # DESeq2 object for downstream DE
#   - vst_normalized_24hrs_with_symbols.csv: VST matrix + gene symbols  # VST matrix with annotations
############################################################

#===========================                                            # Section divider
# Load required libraries                                               # Load dependencies
#===========================
library(DESeq2)    # Differential expression analysis                   # Core DESeq2 functionality
library(tibble)    # Tidy data frames                                   # Modern dataframe class
library(dplyr)     # Data manipulation                                  # Wrangling functions
library(readr)     # Reading TSV files                                  # Fast IO
library(stringr)   # String handling                                    # String manipulation helpers
library(biomaRt)   # Gene annotation via Ensembl                        # Fetch gene symbols

#===========================                                            # Section divider
# Step 1: Set file paths                                                # Input and output paths
#===========================
counts_file <- "/Users/iyad-si/Desktop/KMJ lab/Mccleary data analysis/07_15_2025_4h_24h_bulkrnaseq_cVSp_ISI_analysis/gene_counts_24hours_CvsP copy.txt"  # featureCounts counts file
output_dir  <- "/Users/iyad-si/Desktop/KMJ lab/Mccleary data analysis/07_15_2025_4h_24h_bulkrnaseq_cVSp_ISI_analysis"                                  # Output directory

vst_csv <- file.path(output_dir, "vst_normalized_24hrs_with_symbols.csv")  # VST CSV output path
dds_rds <- file.path(output_dir, "dds_24hrs.rds")                          # RDS DESeq2 object output path

#===========================                                            # Section divider
# Step 2: Load and clean raw counts                                     # Read and format count matrix
#===========================
raw_counts <- read_tsv(counts_file, comment = "#")                        # Read counts, skip header lines starting "#"
sample_cols <- 7:ncol(raw_counts)                                         # Sample columns start at col 7 (first 6 are metadata)

# Clean sample names from BAM paths                                      # Simplify to short labels
colnames(raw_counts)[sample_cols] <- colnames(raw_counts)[sample_cols] %>% 
  basename() %>%                                                          # Keep file name only
  str_remove("_Aligned.*") %>%                                            # Remove aligner suffix
  str_replace("Untreated", "untreated")                                   # Make untreated lowercase for consistency

# Extract count matrix                                                   # Build gene × sample counts
gene_ids <- raw_counts$Geneid                                             # Extract Ensembl gene IDs
count_matrix <- raw_counts[, sample_cols] %>% as.data.frame()             # Subset to samples and convert to data.frame
rownames(count_matrix) <- gene_ids                                        # Set rownames to Ensembl IDs
sample_names <- colnames(count_matrix)                                    # Vector of cleaned sample names

#===========================                                            # Section divider
# Step 3: Build metadata (colData)                                      # Annotate experimental conditions
#===========================
condition <- case_when(                                                   # Map sample names to condition
  str_detect(sample_names, "^C")         ~ "c100",                        # Continuous IL-4 100 ng/mL
  str_detect(sample_names, "^P")         ~ "p100",                        # Pulse IL-4 100 ng/mL
  str_detect(sample_names, "^untreated") ~ "ct",                          # Untreated control
  TRUE ~ NA_character_                                                     # Fallback if unmatched
)
condition <- factor(condition, levels = c("ct", "c100", "p100"))          # Factor with control as reference

coldata <- data.frame(                                                    # Construct colData
  row.names = sample_names,                                               # Sample names as rownames
  condition = condition                                                   # Condition column
)

#===========================                                            # Section divider
# Step 4: DESeq2 normalization                                          # Create DESeq2 object and normalize
#===========================
dds <- DESeqDataSetFromMatrix(                                            # Build DESeq2 dataset
  countData = count_matrix,                                               # Counts input
  colData = coldata,                                                      # Metadata
  design = ~ condition                                                    # Model with condition only
)
dds <- DESeq(dds)                                                         # Run DESeq2 pipeline (size factors, dispersion, GLM)

# Variance stabilizing transformation                                    # Normalize counts with VST
vst <- vst(dds, blind = FALSE)                                            # VST transform, keeping design info
vst_mat <- assay(vst)                                                     # Extract numeric VST matrix

#===========================                                            # Section divider
# Step 5: Gene symbol annotation                                        # Map Ensembl IDs to MGI symbols
#===========================
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")         # Connect to Ensembl, mouse dataset
gene_map <- getBM(                                                        # Retrieve ID ↔ symbol mapping
  attributes = c("ensembl_gene_id", "mgi_symbol"),                        # Columns to fetch
  filters = "ensembl_gene_id",                                            # Filter type
  values = rownames(vst_mat),                                             # Ensembl IDs in dataset
  mart = ensembl                                                          # biomaRt object
)
gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]             # Remove duplicate mappings

# Add gene symbols                                                       # Merge symbols into VST matrix
vst_df <- as.data.frame(vst_mat)                                          # Convert matrix to data.frame
vst_df$ensembl_gene_id <- rownames(vst_df)                                # Preserve Ensembl IDs as a column
vst_df <- left_join(vst_df, gene_map, by = "ensembl_gene_id")             # Join with symbol mapping
vst_df <- vst_df %>% relocate(mgi_symbol, .before = ensembl_gene_id)      # Move symbol column to front

#===========================                                            # Section divider
# Step 6: Save outputs                                                  # Export results
#===========================
write.csv(vst_df, file = vst_csv, row.names = FALSE)                      # Save VST + symbols as CSV
saveRDS(dds, file = dds_rds)                                              # Save DESeq2 object as RDS

