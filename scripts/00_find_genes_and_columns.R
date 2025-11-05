# ==============================================================================
# Find RASSF Genes and Survival Columns
# RASSF Gene Expression & Survival Analysis in TCGA Melanoma
# Author: Samson Kosemani
# Updated: November 2025
# ==============================================================================

cat("=================================================================\n")
cat("FINDING RASSF GENES AND SURVIVAL COLUMNS\n")
cat("Author: Samson Kosemani\n")
cat("Updated: November 2025\n")
cat("=================================================================\n\n")

# Load required library
library(SummarizedExperiment)

# ---- 1. EXPLORE GENE IDENTIFIERS ----
cat("1. EXPLORING GENE IDENTIFIERS IN YOUR DATASET\n")

expression_matrix <- assay(skcm_data)
clinical_data <- colData(skcm_data)

cat("Expression matrix dimensions:", dim(expression_matrix), "\n")
cat("First 20 gene identifiers (row names):\n")
print(head(rownames(expression_matrix), 20))

# Check what type of identifiers we have
gene_ids <- rownames(expression_matrix)
cat("\nGene ID patterns:\n")
cat("First gene ID:", gene_ids[1], "\n")
cat("Sample of gene IDs:\n")
print(sample(gene_ids, 10))

# Check if we have ENSEMBL IDs (common in TCGA)
is_ensembl <- grepl("^ENSG", gene_ids[1])
cat("Using ENSEMBL IDs:", is_ensembl, "\n")

# ---- 2. CHECK ROWDATA FOR GENE SYMBOLS ----
cat("\n2. CHECKING ROWDATA FOR GENE SYMBOLS\n")

if (ncol(rowData(skcm_data)) > 0) {
  rd <- rowData(skcm_data)
  cat("RowData columns available:\n")
  print(colnames(rd))
  
  # Look for gene symbol columns
  symbol_cols <- colnames(rd)[grepl("symbol|gene_name|hgnc|entrez", colnames(rd), ignore.case = TRUE)]
  if (length(symbol_cols) > 0) {
    cat("Potential gene symbol columns:", paste(symbol_cols, collapse = ", "), "\n")
    
    # Check the first symbol column
    symbol_col <- symbol_cols[1]
    cat("Using column:", symbol_col, "\n")
    cat("Sample values from", symbol_col, ":\n")
    print(head(rd[[symbol_col]]))
    
    # Check if RASSF genes are in this column
    rassf_genes <- c("RASSF1", "RASSF2", "RASSF3", "RASSF4", "RASSF5", 
                     "RASSF6", "RASSF7", "RASSF8", "RASSF9", "RASSF10")
    
    available_rassf <- rassf_genes[rassf_genes %in% rd[[symbol_col]]]
    cat("RASSF genes found in", symbol_col, ":", length(available_rassf), "/ 10\n")
    if (length(available_rassf) > 0) {
      cat("Available RASSF genes:", paste(available_rassf, collapse = ", "), "\n")
      
      # Show the corresponding row indices
      for (gene in available_rassf) {
        gene_idx <- which(rd[[symbol_col]] == gene)
        cat("  ", gene, "-> Row", gene_idx[1], "(ID:", rownames(expression_matrix)[gene_idx[1]], ")\n")
      }
    } else {
      cat("No RASSF genes found in", symbol_col, "\n")
    }
  } else {
    cat("No gene symbol columns found in rowData\n")
  }
} else {
  cat("No rowData available\n")
}

# ---- 3. SEARCH FOR RASSF GENES USING PATTERN MATCHING ----
cat("\n3. SEARCHING FOR RASSF GENES USING PATTERN MATCHING\n")

# Try to find RASSF genes by pattern matching in all rowData columns
if (ncol(rowData(skcm_data)) > 0) {
  rd <- rowData(skcm_data)
  
  for (col in colnames(rd)) {
    if (is.character(rd[[col]]) || is.factor(rd[[col]])) {
      # Look for RASSF pattern in this column
      rassf_matches <- grep("RASSF", rd[[col]], ignore.case = TRUE, value = TRUE)
      if (length(rassf_matches) > 0) {
        cat("Found RASSF patterns in column", col, ":\n")
        print(unique(rassf_matches))
      }
    }
  }
}

# ---- 4. CHECK SURVIVAL COLUMNS ----
cat("\n4. IDENTIFYING SURVIVAL COLUMNS\n")

clinical_cols <- colnames(clinical_data)

# Time columns
time_cols <- clinical_cols[grepl("time|day|follow", clinical_cols, ignore.case = TRUE)]
cat("Time-related columns:\n")
print(time_cols)

# Status columns
status_cols <- clinical_cols[grepl("status|vital|dead|alive", clinical_cols, ignore.case = TRUE)]
cat("\nStatus-related columns:\n")
print(status_cols)

# Check specific columns
cat("\nChecking specific survival columns:\n")
for (col in c("days_to_last_follow_up", "vital_status", "days_to_death")) {
  if (col %in% clinical_cols) {
    cat(col, ":", head(clinical_data[[col]]), "\n")
  }
}

# ---- 5. PROVIDE SOLUTIONS ----
cat("\n")
cat(rep("=", 60), "\n")
cat("SOLUTIONS FOR GENE IDENTIFICATION\n")
cat(rep("=", 60), "\n")
cat("\n")

if (exists("available_rassf") && length(available_rassf) > 0) {
  cat("✓ RASSF genes found! Use these in your analysis:\n")
  cat("  Genes:", paste(available_rassf, collapse = ", "), "\n")
  if (exists("symbol_col")) {
    cat("  Symbol column:", symbol_col, "\n")
  }
} else {
  cat("✗ RASSF genes not found by symbol. Alternative approaches:\n")
  cat("  1. Your dataset may use ENSEMBL IDs instead of gene symbols\n")
  cat("  2. We need to map ENSEMBL IDs to gene symbols\n")
  cat("  3. First gene ID in your dataset:", gene_ids[1], "\n")
}

cat("\nRecommended survival columns:\n")
cat("  Time: days_to_last_follow_up\n")
cat("  Status: vital_status\n")

cat("\n")
cat(rep("=", 60), "\n")
cat("NEXT STEPS\n")
cat(rep("=", 60), "\n")
cat("\n")

cat("1. Run this script to see what gene identifiers you have\n")
cat("2. Based on the output, we'll create the appropriate mapping\n")
cat("3. Then run the main analysis with the correct gene mapping\n")
