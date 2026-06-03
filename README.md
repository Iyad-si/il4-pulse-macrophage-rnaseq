# Transient IL-4 stimulation durably activates an alternative macrophage program and modulates the inflammatory response

---

## Overview

This repository contains the R code used to process and analyze bulk RNA-seq data for the manuscript submitted to *iScience*. The analysis integrates DESeq2 differential expression results across 4-hour and 24-hour timepoints comparing transient (15-minute pulse) versus continuous IL-4 stimulation in murine bone marrow-derived macrophages (BMDMs), and generates the union DEG table used for downstream clustering and figure generation.

---

## Data Availability

Raw RNA-seq data are deposited in the NCBI Gene Expression Omnibus (GEO) under accession number **GSE322520**. The DESeq2 objects (`.rds` files) required to run this script are derived from the raw counts in that repository.

---

## Repository Contents

```
il4-transient-macrophage-program/
├── README.md
└── scripts/
    └── DDS_to_Union_Table.R      # Main differential expression integration script
```

---

## Script Description: `DDS_to_Union_Table.R`

This script takes pre-built DESeq2 objects from the 4h and 24h IL-4 pulse vs. continuous experiments and performs the following steps:

1. **Expression filtering** — retains genes with non-zero counts in at least one stimulated condition (p100 or c100) at either timepoint
2. **Contrast extraction** — pulls log2FC and adjusted p-value for four contrasts: p100 vs. control and c100 vs. control at 4h and 24h
3. **Gene symbol mapping** — maps Ensembl IDs to MGI gene symbols via biomaRt
4. **OR-gate DEG selection** — a gene is significant if |log2FC| > 1 and padj < 0.05 at *either* 4h or 24h (applied separately for pulse and continuous)
5. **Union table construction** — merges pulse and continuous DEG lists with pulse rows winning ties for genes in both
6. **Symbol cleaning (Option A)** — removes uninformative symbols (ENSMUSG..., Gm..., ...Rik, LOC...)

**Outputs:**
- Master expressed-gene table (log2FC + padj, all genes)
- Pulse-significant DEG table
- Continuous-significant DEG table
- UNION table (pulse wins ties, clean symbols)
- Uncapped log2FC matrix for downstream clustering
- Removed uninformative symbol list

---

## Software Requirements

| Software | Version |
|----------|---------|
| R | 4.4.1 |
| DESeq2 | 1.42.0 |
| tidyverse | — |
| biomaRt | 2.58.0 |

---

## Usage

Before running the script, update the user-defined paths in **Section 0** at the top of `DDS_to_Union_Table.R`:

```r
base_dir     <- "/path/to/your/analysis/directory"
dds_4hrs_rds  <- "/path/to/dds_4hrs.rds"
dds_24hrs_rds <- "/path/to/dds_24hrs.rds"
```

All other parameters (thresholds, output filenames) are set in the same section and can be adjusted as needed.

---

## License

This code is provided for reproducibility purposes in association with the manuscript. Please cite the paper if you use or adapt this code.
