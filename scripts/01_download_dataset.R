#!/usr/bin/env Rscript

# 01_download_dataset.R
# Purpose: Load the Bioconductor airway dataset and save standardized inputs
# for downstream RNA-seq differential expression analysis.

suppressPackageStartupMessages({
  library(airway)
  library(SummarizedExperiment)
  library(yaml)
})

# Load config
config <- yaml::read_yaml("config/config.yaml")

# Output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading Bioconductor airway dataset...")

data("airway")

# Save full SummarizedExperiment object
saveRDS(
  airway,
  file = "data/processed/airway.rds"
)

# Extract count matrix and sample metadata
count_matrix <- assay(airway)
sample_metadata <- as.data.frame(colData(airway))

# Save processed inputs for transparency
write.csv(
  count_matrix,
  file = "data/processed/airway_counts.csv"
)

write.csv(
  sample_metadata,
  file = "data/processed/airway_sample_metadata.csv",
  row.names = TRUE
)

# Save dataset summary
dataset_summary <- data.frame(
  dataset = "airway",
  source = "Bioconductor airway package",
  genes = nrow(count_matrix),
  samples = ncol(count_matrix),
  design = "paired treated vs untreated",
  model = "design = ~ cell + dex",
  padj_cutoff = config$padj_cutoff,
  log2fc_cutoff = config$log2fc_cutoff
)

write.csv(
  dataset_summary,
  file = "results/tables/dataset_summary.csv",
  row.names = FALSE
)

message("Dataset loading completed.")
message("Genes: ", nrow(count_matrix))
message("Samples: ", ncol(count_matrix))
message("Saved: data/processed/airway.rds")
message("Saved: data/processed/airway_counts.csv")
message("Saved: data/processed/airway_sample_metadata.csv")
message("Saved: results/tables/dataset_summary.csv")
