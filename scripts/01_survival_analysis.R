# ==============================================================================
# RASSF Gene Expression & Survival Analysis in TCGA Melanoma
# Using cBioPortal TCGA PanCancer Atlas Data
# Author: Samson Kosemani
# Updated: November 2025
# ==============================================================================

cat("=================================================================\n")
cat("RASSF GENE SURVIVAL ANALYSIS IN SKCM TCGA\n")
cat("Author: Samson Kosemani\n")
cat("Updated: November 2025\n")
cat("Data Source: cBioPortal TCGA PanCancer Atlas\n")
cat("=================================================================\n\n")

# Load required libraries
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

# ---- 1. DATA EXTRACTION BASED ON EXPLORATION FINDINGS ----
cat("1. EXTRACTING DATA FROM SUMMARIZEDEXPERIMENT\n")

# Extract components (based on exploration results)
expression_matrix <- assay(skcm_data)
clinical_data <- colData(skcm_data)

cat("Expression matrix dimensions:", dim(expression_matrix), "\n")
cat("Clinical data dimensions:", dim(clinical_data), "\n")

# ---- 2. MANUALLY SET THE CORRECT COLUMN NAMES ----
# REPLACE THESE WITH THE COLUMN NAMES FOUND IN YOUR EXPLORATION
# Check the output from 00_setup_and_explore.R and update these:

time_column <- "OS.time"    # CHANGE THIS to your actual time column
status_column <- "OS"       # CHANGE THIS to your actual status column

cat("Using time column:", time_column, "\n")
cat("Using status column:", status_column, "\n")

# Verify columns exist
if (!time_column %in% colnames(clinical_data)) {
  stop("Time column '", time_column, "' not found. Available columns: ", 
       paste(colnames(clinical_data), collapse = ", "))
}

if (!status_column %in% colnames(clinical_data)) {
  stop("Status column '", status_column, "' not found. Available columns: ", 
       paste(colnames(clinical_data), collapse = ", "))
}

# ---- 3. DEFINE RASSF GENES ----
# UPDATE THIS LIST BASED ON WHAT WAS FOUND IN EXPLORATION
rassf_genes <- c("RASSF1", "RASSF2", "RASSF3", "RASSF4", "RASSF5", 
                 "RASSF6", "RASSF7", "RASSF8", "RASSF9", "RASSF10")

# Find available RASSF genes
available_rassf <- rassf_genes[rassf_genes %in% rownames(expression_matrix)]
cat("RASSF genes to analyze:", length(available_rassf), "/ 10\n")
cat("Genes:", paste(available_rassf, collapse = ", "), "\n")

if (length(available_rassf) == 0) {
  stop("No RASSF genes found. Check gene identifiers in your dataset.")
}

# ---- 4. CREATE SURVIVAL DATASET ----
cat("\n2. CREATING SURVIVAL DATASET\n")

# Create survival dataset
survival_data <- data.frame(
  patient_id = colnames(expression_matrix),
  time = as.numeric(clinical_data[[time_column]]),
  status = as.numeric(clinical_data[[status_column]])
)

# Add RASSF gene expression data
for (gene in available_rassf) {
  gene_idx <- which(rownames(expression_matrix) == gene)
  survival_data[[gene]] <- as.numeric(expression_matrix[gene_idx, ])
}

# Data cleaning
initial_patients <- nrow(survival_data)
survival_data_clean <- survival_data %>%
  filter(complete.cases(time, status)) %>%
  filter(time > 0) %>%
  filter(status %in% c(0, 1))

final_patients <- nrow(survival_data_clean)
cat("Data quality filtering:\n")
cat("  Initial patients:", initial_patients, "\n")
cat("  Final patients:", final_patients, "\n")
cat("  Patients excluded:", initial_patients - final_patients, "\n")
cat("  Event rate:", round(100 * mean(survival_data_clean$status), 1), "%\n")

# Save processed data
saveRDS(survival_data_clean, "data/processed/survival_data_clean.rds")
write.csv(survival_data_clean, "data/processed/survival_data_clean.csv", row.names = FALSE)

# ---- 5. SURVIVAL ANALYSIS FUNCTIONS ----
cat("\n3. SETTING UP SURVIVAL ANALYSIS\n")

# Expression categorization function
categorize_expression <- function(gene_expression) {
  quantiles <- quantile(gene_expression, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  cut(gene_expression, breaks = quantiles, labels = c("Low", "Medium", "High"), 
      include.lowest = TRUE)
}

# KM plot function with author branding
create_km_plot <- function(gene, data, output_dir = "results/figures") {
  
  # Prepare data for this gene
  gene_data <- data.frame(
    time = data$time,
    status = data$status,
    expression = data[[gene]]
  )
  
  gene_data <- gene_data[complete.cases(gene_data), ]
  gene_data$expression_group <- categorize_expression(gene_data$expression)
  gene_data <- gene_data[complete.cases(gene_data), ]
  
  if (nrow(gene_data) < 10) {
    cat("  Insufficient data for", gene, "(n =", nrow(gene_data), ")\n")
    return(NULL)
  }
  
  # Survival analysis
  surv_obj <- Surv(time = gene_data$time, event = gene_data$status)
  km_fit <- survfit(surv_obj ~ expression_group, data = gene_data)
  logrank_test <- survdiff(surv_obj ~ expression_group, data = gene_data)
  p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # Create plot
  km_plot <- ggsurvplot(
    km_fit,
    data = gene_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    legend = "top",
    legend.title = paste(gene, "Expression"),
    legend.labs = levels(gene_data$expression_group),
    palette = c("Low" = "#E41A1C", "Medium" = "#377EB8", "High" = "#4DAF4A"),
    title = paste("Overall Survival by", gene, "Expression in SKCM"),
    subtitle = "Analysis: Samson Kosemani | Data: cBioPortal TCGA PanCancer Atlas",
    xlab = "Time (Days)",
    ylab = "Overall Survival Probability",
    ggtheme = theme_minimal()
  )
  
  km_plot$plot <- km_plot$plot +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  
  # Save plot
  png_file <- file.path(output_dir, paste0("KM_", gene, "_SKCM.png"))
  png(png_file, width = 1000, height = 800, res = 150)
  print(km_plot)
  dev.off()
  
  pdf_file <- file.path(output_dir, paste0("KM_", gene, "_SKCM.pdf"))
  pdf(pdf_file, width = 10, height = 8)
  print(km_plot)
  dev.off()
  
  return(list(
    gene = gene,
    p_value = p_value,
    n_patients = nrow(gene_data),
    plot = km_plot,
    group_sizes = table(gene_data$expression_group)
  ))
}

# ---- 6. EXECUTE SURVIVAL ANALYSIS ----
cat("\n4. EXECUTING SURVIVAL ANALYSIS FOR ALL RASSF GENES\n")

results <- list()
significant_genes <- c()

for (gene in available_rassf) {
  cat("Analyzing", gene, "... ")
  result <- create_km_plot(gene, survival_data_clean)
  
  if (!is.null(result)) {
    results[[gene]] <- result
    if (result$p_value < 0.05) {
      significant_genes <- c(significant_genes, gene)
      cat("SIGNIFICANT (p =", round(result$p_value, 4), ")\n")
    } else {
      cat("not significant (p =", round(result$p_value, 4), ")\n")
    }
  } else {
    cat("FAILED - insufficient data\n")
  }
}

# ---- 7. CREATE SUMMARY RESULTS ----
cat("\n5. CREATING SUMMARY TABLES AND VISUALIZATIONS\n")

# Create results table
results_table <- data.frame(
  Gene = names(results),
  P_Value = sapply(results, function(x) x$p_value),
  Patients = sapply(results, function(x) x$n_patients),
  Significant = sapply(results, function(x) ifelse(x$p_value < 0.05, "Yes", "No")),
  Low_Group = sapply(results, function(x) x$group_sizes["Low"]),
  Medium_Group = sapply(results, function(x) x$group_sizes["Medium"]),
  High_Group = sapply(results, function(x) x$group_sizes["High"])
)

# Format p-values
results_table$P_Value_Formatted <- ifelse(
  results_table$P_Value < 0.001, 
  "< 0.001", 
  sprintf("%.4f", results_table$P_Value)
)

# Save results
write.csv(results_table, "results/tables/rassf_survival_results.csv", row.names = FALSE)
saveRDS(results_table, "results/tables/rassf_survival_results.rds")

# Create summary visualization
pval_plot <- ggplot(results_table, aes(x = reorder(Gene, -P_Value), y = -log10(P_Value), 
                                       fill = Significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Yes" = "#E41A1C", "No" = "#377EB8")) +
  labs(
    title = "RASSF Genes Survival Analysis P-Values in SKCM",
    subtitle = "Analysis: Samson Kosemani | Data: cBioPortal TCGA PanCancer Atlas",
    x = "Gene", 
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
  )

ggsave("results/figures/rassf_pvalue_summary.png", pval_plot, width = 10, height = 6)
ggsave("results/figures/rassf_pvalue_summary.pdf", pval_plot, width = 10, height = 6)

# ---- 8. FINAL SUMMARY ----
cat("\n" + rep("=", 60) + "\n")
cat("ANALYSIS COMPLETE - SUMMARY RESULTS\n")
cat(rep("=", 60) + "\n\n")

cat("Significant genes found:", length(significant_genes), "/", length(results), "\n")
if (length(significant_genes) > 0) {
  cat("Significant genes:", paste(significant_genes, collapse = ", "), "\n")
}

cat("\nDetailed Results:\n")
print(results_table)

cat("\nOutput Files Created:\n")
cat("• results/figures/KM_[GENE]_SKCM.png/pdf - Individual KM plots\n")
cat("• results/figures/rassf_pvalue_summary.png/pdf - Summary plot\n")
cat("• results/tables/rassf_survival_results.csv/rds - Results table\n")
cat("• data/processed/survival_data_clean.csv/rds - Processed data\n")

cat("\n" + rep("=", 60) + "\n")
cat("RASSF SURVIVAL ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Author: Samson Kosemani\n")
cat("Date:", Sys.Date(), "\n")
cat(rep("=", 60) + "\n")




# Quick exploration to find survival columns
cat("=== FINDING SURVIVAL COLUMNS ===\n")

# Look at all column names that might be related to survival
clinical_cols <- colnames(clinical_data)
cat("All clinical columns (first 30):\n")
print(clinical_cols[1:30])

# Look for time-related columns
time_cols <- clinical_cols[grepl("time|day|survival|follow|death", clinical_cols, ignore.case = TRUE)]
cat("\nTime-related columns:\n")
print(time_cols)

# Look for status-related columns  
status_cols <- clinical_cols[grepl("status|event|death|vital|alive|dead", clinical_cols, ignore.case = TRUE)]
cat("\nStatus-related columns:\n")
print(status_cols)

# Check specific columns that are commonly used
common_time_cols <- c("days_to_death", "days_to_last_followup", "OS.time")
common_status_cols <- c("vital_status", "OS", "OS_STATUS")

cat("\nChecking common TCGA survival columns:\n")
for (col in common_time_cols) {
  if (col %in% clinical_cols) {
    cat("Found time column:", col, "\n")
    cat("Sample values:", head(clinical_data[[col]]), "\n\n")
  }
}

for (col in common_status_cols) {
  if (col %in% clinical_cols) {
    cat("Found status column:", col, "\n")
    cat("Sample values:", head(clinical_data[[col]]), "\n\n")
  }
}
# ==============================================================================
# RASSF Gene Expression & Survival Analysis in TCGA Melanoma
# Using cBioPortal TCGA PanCancer Atlas Data
# Author: Samson Kosemani
# Updated: November 2025
# ==============================================================================

cat("=================================================================\n")
cat("RASSF GENE SURVIVAL ANALYSIS IN SKCM TCGA\n")
cat("Author: Samson Kosemani\n")
cat("Updated: November 2025\n")
cat("Data Source: cBioPortal TCGA PanCancer Atlas\n")
cat("=================================================================\n\n")

# Load required libraries
library(survival)
library(survminer)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

# ---- 1. DATA EXTRACTION ----
cat("1. EXTRACTING DATA FROM SUMMARIZEDEXPERIMENT\n")

expression_matrix <- assay(skcm_data)
clinical_data <- colData(skcm_data)

cat("Expression matrix dimensions:", dim(expression_matrix), "\n")
cat("Clinical data dimensions:", dim(clinical_data), "\n")

# ---- 2. USE CORRECT SURVIVAL COLUMNS ----
# Based on your dataset structure, use these columns:
time_column <- "days_to_last_follow_up"    # Time to last follow-up
status_column <- "vital_status"            # Vital status

cat("Using time column:", time_column, "\n")
cat("Using status column:", status_column, "\n")

# Verify columns exist
if (!time_column %in% colnames(clinical_data)) {
  stop("Time column '", time_column, "' not found.")
}

if (!status_column %in% colnames(clinical_data)) {
  stop("Status column '", status_column, "' not found.")
}

# Check the values in these columns
cat("\nTime column sample values:\n")
print(head(clinical_data[[time_column]]))
cat("\nStatus column sample values:\n")
print(head(clinical_data[[status_column]]))

# ---- 3. DEFINE RASSF GENES ----
rassf_genes <- c("RASSF1", "RASSF2", "RASSF3", "RASSF4", "RASSF5", 
                 "RASSF6", "RASSF7", "RASSF8", "RASSF9", "RASSF10")

# Find available RASSF genes
available_rassf <- rassf_genes[rassf_genes %in% rownames(expression_matrix)]
cat("\nRASSF genes to analyze:", length(available_rassf), "/ 10\n")
cat("Genes:", paste(available_rassf, collapse = ", "), "\n")

if (length(available_rassf) == 0) {
  stop("No RASSF genes found. Check gene identifiers in your dataset.")
}

# ---- 4. CREATE SURVIVAL DATASET WITH PROPER CONVERSION ----
cat("\n2. CREATING SURVIVAL DATASET\n")

# Function to convert vital status to numeric
convert_vital_status <- function(status_vector) {
  # Convert "Alive"/"Dead" to 0/1
  if (is.character(status_vector) || is.factor(status_vector)) {
    status_numeric <- ifelse(grepl("dead|deceased", status_vector, ignore.case = TRUE), 1, 0)
    return(status_numeric)
  }
  return(as.numeric(status_vector))
}

# Create survival dataset
survival_data <- data.frame(
  patient_id = colnames(expression_matrix),
  time = as.numeric(clinical_data[[time_column]]),
  status = convert_vital_status(clinical_data[[status_column]])
)

# Add RASSF gene expression data
for (gene in available_rassf) {
  gene_idx <- which(rownames(expression_matrix) == gene)
  survival_data[[gene]] <- as.numeric(expression_matrix[gene_idx, ])
}

# Data cleaning
initial_patients <- nrow(survival_data)
survival_data_clean <- survival_data %>%
  filter(complete.cases(time, status)) %>%
  filter(time > 0) %>%
  filter(status %in% c(0, 1))

final_patients <- nrow(survival_data_clean)
cat("Data quality filtering:\n")
cat("  Initial patients:", initial_patients, "\n")
cat("  Final patients:", final_patients, "\n")
cat("  Patients excluded:", initial_patients - final_patients, "\n")
cat("  Event rate:", round(100 * mean(survival_data_clean$status), 1), "%\n")

# Check if we have enough data
if (final_patients < 10) {
  stop("Insufficient data after filtering. Only", final_patients, "patients available.")
}

# Save processed data
saveRDS(survival_data_clean, "data/processed/survival_data_clean.rds")
write.csv(survival_data_clean, "data/processed/survival_data_clean.csv", row.names = FALSE)

# ---- 5. SURVIVAL ANALYSIS FUNCTIONS ----
cat("\n3. SETTING UP SURVIVAL ANALYSIS\n")

# Expression categorization function
categorize_expression <- function(gene_expression) {
  quantiles <- quantile(gene_expression, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  cut(gene_expression, breaks = quantiles, labels = c("Low", "Medium", "High"), 
      include.lowest = TRUE)
}

# KM plot function with author branding
create_km_plot <- function(gene, data, output_dir = "results/figures") {
  
  # Prepare data for this gene
  gene_data <- data.frame(
    time = data$time,
    status = data$status,
    expression = data[[gene]]
  )
  
  gene_data <- gene_data[complete.cases(gene_data), ]
  gene_data$expression_group <- categorize_expression(gene_data$expression)
  gene_data <- gene_data[complete.cases(gene_data), ]
  
  if (nrow(gene_data) < 10) {
    cat("  Insufficient data for", gene, "(n =", nrow(gene_data), ")\n")
    return(NULL)
  }
  
  # Survival analysis
  surv_obj <- Surv(time = gene_data$time, event = gene_data$status)
  km_fit <- survfit(surv_obj ~ expression_group, data = gene_data)
  logrank_test <- survdiff(surv_obj ~ expression_group, data = gene_data)
  p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # Create plot
  km_plot <- ggsurvplot(
    km_fit,
    data = gene_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    legend = "top",
    legend.title = paste(gene, "Expression"),
    legend.labs = levels(gene_data$expression_group),
    palette = c("Low" = "#E41A1C", "Medium" = "#377EB8", "High" = "#4DAF4A"),
    title = paste("Overall Survival by", gene, "Expression in SKCM"),
    subtitle = "Analysis: Samson Kosemani | Data: cBioPortal TCGA PanCancer Atlas",
    xlab = "Time to Last Follow-up (Days)",
    ylab = "Overall Survival Probability",
    ggtheme = theme_minimal()
  )
  
  km_plot$plot <- km_plot$plot +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  
  # Save plot
  png_file <- file.path(output_dir, paste0("KM_", gene, "_SKCM.png"))
  png(png_file, width = 1000, height = 800, res = 150)
  print(km_plot)
  dev.off()
  
  return(list(
    gene = gene,
    p_value = p_value,
    n_patients = nrow(gene_data),
    plot = km_plot,
    group_sizes = table(gene_data$expression_group)
  ))
}

# ---- 6. EXECUTE SURVIVAL ANALYSIS ----
cat("\n4. EXECUTING SURVIVAL ANALYSIS FOR ALL RASSF GENES\n")

results <- list()
significant_genes <- c()

for (gene in available_rassf) {
  cat("Analyzing", gene, "... ")
  result <- create_km_plot(gene, survival_data_clean)
  
  if (!is.null(result)) {
    results[[gene]] <- result
    if (result$p_value < 0.05) {
      significant_genes <- c(significant_genes, gene)
      cat("SIGNIFICANT (p =", round(result$p_value, 4), ")\n")
    } else {
      cat("not significant (p =", round(result$p_value, 4), ")\n")
    }
  } else {
    cat("FAILED - insufficient data\n")
  }
}

# ---- 7. CREATE SUMMARY RESULTS ----
cat("\n5. CREATING SUMMARY TABLES AND VISUALIZATIONS\n")

# Create results table
results_table <- data.frame(
  Gene = names(results),
  P_Value = sapply(results, function(x) x$p_value),
  Patients = sapply(results, function(x) x$n_patients),
  Significant = sapply(results, function(x) ifelse(x$p_value < 0.05, "Yes", "No")),
  Low_Group = sapply(results, function(x) x$group_sizes["Low"]),
  Medium_Group = sapply(results, function(x) x$group_sizes["Medium"]),
  High_Group = sapply(results, function(x) x$group_sizes["High"])
)

# Save results
write.csv(results_table, "results/tables/rassf_survival_results.csv", row.names = FALSE)

cat("\n" + rep("=", 60) + "\n")
cat("ANALYSIS COMPLETE - SUMMARY RESULTS\n")
cat(rep("=", 60) + "\n\n")

cat("Significant genes found:", length(significant_genes), "/", length(results), "\n")
if (length(significant_genes) > 0) {
  cat("Significant genes:", paste(significant_genes, collapse = ", "), "\n")
}

cat("\nDetailed Results:\n")
print(results_table)

cat("\n" + rep("=", 60) + "\n")
cat("RASSF SURVIVAL ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Author: Samson Kosemani\n")
cat("Date:", Sys.Date(), "\n")
cat(rep("=", 60) + "\n")