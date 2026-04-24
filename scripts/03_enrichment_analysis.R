# 03_enrichment_analysis.R
# Functional enrichment analysis using enrichR
# Input: DESeq2 results table
# Output: GO and KEGG enrichment tables + GO barplot

.libPaths(c("/home/mariano/R/library", .libPaths()))

# Load libraries
library(enrichR)
library(ggplot2)

# Load DESeq2 results
res_df <- read.csv("results/tables/deseq2_results.csv")

# Filter significant genes
sig_genes <- res_df[
  !is.na(res_df$padj) & res_df$padj < 0.05,
]

# Load airway annotation again to map Ensembl IDs to gene symbols
library(DESeq2)
library(airway)

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

gene_info <- as.data.frame(rowData(dds))[, c("gene_id", "symbol")]

# Merge DESeq2 results with gene symbols
sig_genes <- merge(
  sig_genes,
  gene_info,
  by = "gene_id"
)

# Clean gene symbol list
gene_symbols <- unique(sig_genes$symbol)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]

# Run enrichment analysis
selected_dbs <- c(
  "GO_Biological_Process_2021",
  "KEGG_2021_Human"
)

enrich_results <- enrichr(gene_symbols, selected_dbs)

# Save enrichment tables
write.csv(
  enrich_results[["GO_Biological_Process_2021"]],
  "results/tables/go_enrichment.csv",
  row.names = FALSE
)

write.csv(
  enrich_results[["KEGG_2021_Human"]],
  "results/tables/kegg_enrichment.csv",
  row.names = FALSE
)

# Create GO barplot
go <- enrich_results[["GO_Biological_Process_2021"]]
go_sorted <- go[order(go$Adjusted.P.value), ]
top_go <- go_sorted[1:10, ]

go_plot <- ggplot(
  top_go,
  aes(x = reorder(Term, Combined.Score), y = Combined.Score)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top GO Biological Processes",
    x = "GO Term",
    y = "Combined Score"
  ) +
  theme_minimal()

ggsave(
  "results/figures/go_barplot.pdf",
  plot = go_plot,
  width = 8,
  height = 6
)

# Print short summary
cat("Enrichment analysis completed.\n")
cat("Number of significant genes:", nrow(sig_genes), "\n")
cat("Number of unique gene symbols:", length(gene_symbols), "\n")
