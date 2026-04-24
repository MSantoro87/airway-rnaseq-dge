# 02_deseq2_analysis.R
# Differential expression analysis using DESeq2 on airway dataset

# Set library path
.libPaths(c("/home/mariano/R/library", .libPaths()))

# Load libraries
library(DESeq2)
library(airway)

# Load dataset
data("airway")

# Create DESeq2 dataset
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

# Set reference level
dds$dex <- relevel(dds$dex, ref = "untrt")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Save results
write.csv(
  res_df,
  "results/tables/deseq2_results.csv",
  row.names = FALSE
)

# Volcano plot
library(EnhancedVolcano)

pdf("results/figures/volcano_plot.pdf", width = 8, height = 7)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Differential Expression: Dexamethasone Treatment",
  subtitle = "Airway RNA-seq dataset"
)

dev.off()

# Heatmap (top 50 genes)
library(pheatmap)

resOrdered <- res[order(res$pvalue), ]
top_genes <- rownames(resOrdered)[1:50]

normalized_counts <- counts(dds, normalized = TRUE)
mat <- normalized_counts[top_genes, ]

pdf("results/figures/heatmap_top50.pdf", width = 8, height = 10)

pheatmap(
  mat,
  scale = "row",
  show_rownames = FALSE,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()
