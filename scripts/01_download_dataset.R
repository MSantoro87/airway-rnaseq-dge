# 01_download_dataset.R
# Goal: Download and inspect GEO dataset GSE183947

.libPaths(c("/home/mariano/R/library", .libPaths()))

# Install GEOquery if missing
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", lib = "/home/mariano/R/library")
}

library(GEOquery)

# Download GEO metadata
gse <- getGEO("GSE183947", GSEMatrix = TRUE)

# Inspect object
print(gse)

# Save metadata object
saveRDS(gse, file = "data/raw/GSE183947_geo_metadata.rds")
