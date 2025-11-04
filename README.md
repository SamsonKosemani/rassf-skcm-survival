cat > README.md << 'EOF'
# RASSF Gene Expression & Survival Analysis in TCGA Melanoma

## Project Overview
Analysis of RASSF tumor suppressor gene family expression and overall survival in cutaneous melanoma using TCGA SKCM data from cBioPortal.

**Author**: Samson Kosemani  
**Data Source**: cBioPortal TCGA PanCancer Atlas (SKCM)  
**Analysis Date**: 2024-2025

## Project Structure
- `scripts/` - Analysis R scripts
- `data/processed/` - Processed datasets and results
- `results/figures/` - Kaplan-Meier plots and visualizations
- `results/tables/` - Statistical results
- `docs/` - Documentation and analysis summaries

## Analysis Description
Survival analysis of 10 RASSF genes (RASSF1-RASSF10) in TCGA SKCM cohort using Kaplan-Meier curves and log-rank tests.

## Key Files
- `scripts/01_survival_analysis_simple.R` - Main analysis script
- `results/figures/KM_*_SKCM.png` - Individual Kaplan-Meier plots
- `data/processed/rassf_survival_results.csv` - Statistical results
- `docs/ANALYSIS_SUMMARY.md` - Comprehensive analysis summary

## Results
Analysis completed for all 10 RASSF genes with individual survival curves and statistical significance testing.
