#!/usr/bin/env Rscript

# 03_enrichment_analysis.R
# Purpose: Run GO and KEGG enrichment analysis using config-driven databases.

suppressPackageStartupMessages({
  library(enrichR)
  library(DESeq2)
  library(airway)
  library(yaml)
})

# Load config
config <- yaml::read_yaml("config/config.yaml")

enrichment_databases <- config$enrichment_databases

# Output directories
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading significant DESeq2 results...")

sig_genes <- read.csv("results/tables/deseq2_results_significant.csv")

if (nrow(sig_genes) == 0) {
  stop("No significant genes found in results/tables/deseq2_results_significant.csv")
}

message("Loading gene annotation from airway dataset...")

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

gene_info <- as.data.frame(rowData(dds))[, c("gene_id", "symbol")]

sig_genes <- merge(
  sig_genes,
  gene_info,
  by = "gene_id"
)

gene_symbols <- unique(sig_genes$symbol)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]

if (length(gene_symbols) == 0) {
  stop("No valid gene symbols found for enrichment analysis.")
}

message("Running enrichR analysis using databases:")
message(paste(enrichment_databases, collapse = ", "))

enrich_results <- enrichr(gene_symbols, enrichment_databases)

go_results <- enrich_results[["GO_Biological_Process_2021"]]
kegg_results <- enrich_results[["KEGG_2021_Human"]]

write.csv(
  go_results,
  "results/tables/go_enrichment.csv",
  row.names = FALSE
)

write.csv(
  kegg_results,
  "results/tables/kegg_enrichment.csv",
  row.names = FALSE
)

message("Enrichment analysis completed.")
message("Significant genes used: ", nrow(sig_genes))
message("Unique gene symbols used: ", length(gene_symbols))
message("GO terms returned: ", nrow(go_results))
message("KEGG terms returned: ", nrow(kegg_results))
