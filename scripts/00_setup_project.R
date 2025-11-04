# Create project directory structure with author metadata
create_project_structure <- function(author = "Samson Kosemani", update_date = "November 2025") {
  # Create directories
  dirs <- c(
    "data/raw",
    "data/processed", 
    "scripts",
    "results/figures",
    "results/tables",
    "docs",
    "reports"
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created:", dir, "\n")
    }
  }
  
  # Create author metadata
  author_info <- list(
    author = author,
    update_date = update_date,
    project = "RASSF Gene Expression & Survival Analysis in TCGA Melanoma",
    data_source = "cBioPortal TCGA PanCancer Atlas",
    created_date = Sys.Date()
  )
  
  saveRDS(author_info, "docs/author_metadata.rds")
  cat("Author metadata saved\n")
}

create_project_structure("Samson Kosemani", "November 2025")