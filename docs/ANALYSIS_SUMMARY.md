# RASSF Gene Survival Analysis - Final Summary

## Project Overview
- **Analysis**: Association between RASSF gene expression and overall survival in cutaneous melanoma
- **Author**: Samson Kosemani
- **Date**: 2024-2025
- **Data Source**: cBioPortal TCGA PanCancer Atlas (SKCM)

## Methods

### Data Source
- TCGA Skin Cutaneous Melanoma (SKCM) dataset
- 473 initial patients, 149 after quality filtering
- RNA-seq expression data for 10 RASSF genes

### Analysis Approach
1. **Data Extraction**: Expression and clinical data from SummarizedExperiment object
2. **Quality Control**: Removal of incomplete survival data
3. **Gene Mapping**: ENSEMBL IDs mapped to gene symbols via `gene_name` column
4. **Stratification**: Tertile-based grouping (Low/Medium/High expression)
5. **Survival Analysis**: Kaplan-Meier curves with log-rank tests

### Statistical Methods
- Kaplan-Meier survival estimation
- Log-rank test for group comparisons  
- Significance threshold: p < 0.05
- Multiple testing: Unadjusted (exploratory analysis)

## Key Findings

### Significant Genes
- RASSF2

### All Results
```
    Gene    P_Value Patients Significant Low_Group Medium_Group High_Group
  RASSF1 0.78853163      149          No        50           49         50
  RASSF7 0.49497674      149          No        50           49         50
  RASSF2 0.02720452      149         Yes        50           49         50
  RASSF4 0.54368931      149          No        50           49         50
  RASSF8 0.11611840      149          No        50           49         50
  RASSF3 0.21897124      149          No        50           49         50
  RASSF6 0.93316833      149          No        52           49         48
 RASSF10 0.34442788      149          No        56           46         47
  RASSF9 0.44066880      149          No        50           49         50
  RASSF5 0.42062647      149          No        50           49         50
 P_Value_Formatted
            0.7885
            0.4950
            0.0272
            0.5437
            0.1161
            0.2190
            0.9332
            0.3444
            0.4407
            0.4206
```

## Files Generated

### Analysis Scripts
- `00_setup_and_explore.R` - Data exploration and setup
- `01_survival_analysis_simple.R` - Main survival analysis

### Results
- `results/figures/KM_[GENE]_SKCM.png/pdf` - Individual Kaplan-Meier plots
- `results/figures/rassf_pvalue_summary.png` - P-value summary plot
- `results/tables/rassf_survival_results.csv` - Complete results table

### Data
- `data/processed/survival_data_clean.csv` - Processed survival dataset
- `data/processed/rassf_gene_mapping.csv` - Gene identifier mapping

## Interpretation

The analysis identified 1 RASSF genes with significant associations with overall survival in cutaneous melanoma. These findings suggest potential roles for these genes in melanoma progression and patient outcomes.

## Limitations
- Exploratory analysis without multiple testing correction
- TCGA dataset limitations (retrospective, selection bias)
- Sample size reduced due to data quality filtering

## Reproducibility
All code and processed data are available for complete reproducibility of the analysis.

---
*Analysis conducted by Samson Kosemani using R and TCGA data from cBioPortal*
