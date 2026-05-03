#!/usr/bin/env Rscript

# 02_deseq2_analysis.R
# Purpose: Perform differential expression analysis on the airway RNA-seq dataset.

suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
})

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading airway dataset...")

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

dds$dex <- relevel(dds$dex, ref = "untrt")

message("Running DESeq2 model: design = ~ cell + dex")

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("dex", "trt", "untrt"),
  alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df[order(res_df$padj), ]

sig_df <- subset(
  res_df,
  !is.na(padj) & padj < 0.05
)

up_df <- subset(
  sig_df,
  log2FoldChange > 1
)

down_df <- subset(
  sig_df,
  log2FoldChange < -1
)

saveRDS(dds, "data/processed/deseq2_dds.rds")
saveRDS(res, "data/processed/deseq2_results.rds")

write.csv(res_df, "results/tables/deseq2_results_all.csv", row.names = FALSE)
write.csv(sig_df, "results/tables/deseq2_results_significant.csv", row.names = FALSE)
write.csv(up_df, "results/tables/deseq2_results_upregulated.csv", row.names = FALSE)
write.csv(down_df, "results/tables/deseq2_results_downregulated.csv", row.names = FALSE)

message("DESeq2 analysis completed.")
message("Total genes tested: ", nrow(res_df))
message("Significant genes, padj < 0.05: ", nrow(sig_df))
message("Upregulated genes, padj < 0.05 and log2FC > 1: ", nrow(up_df))
message("Downregulated genes, padj < 0.05 and log2FC < -1: ", nrow(down_df))
