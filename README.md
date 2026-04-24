# RNA-seq Differential Expression Analysis (Airway Dataset)

## 📌 Overview
This project performs a complete RNA-seq differential expression analysis workflow using the **airway dataset** from Bioconductor. The goal is to identify gene expression changes in human airway smooth muscle cells following **dexamethasone treatment**.

The project demonstrates a reproducible computational biology pipeline including:
- Differential expression analysis (DESeq2)
- Data visualization (volcano plot, heatmap)
- Functional enrichment analysis (GO, KEGG)

---

## 🎯 Objective
To analyze transcriptional changes induced by dexamethasone and identify biologically relevant pathways associated with treatment response.

---

## 🧬 Dataset
- Source: Bioconductor `airway` package  
- Samples: 8 (treated vs untreated)  
- Data type: RNA-seq count data  

---

## ⚙️ Environment
- OS: Ubuntu 20.04  
- R version: 4.5.2  

### Key packages
- DESeq2  
- ggplot2  
- pheatmap  
- EnhancedVolcano  
- enrichR  

---

## 🧪 Workflow

### 1. Data Loading & Exploration
- Loaded dataset and inspected metadata
- Identified experimental variables (`cell`, `dex`)

### 2. Differential Expression Analysis
- Built DESeq2 model: `~ cell + dex`
- Controlled for donor variability
- Performed statistical testing

### 3. Visualization
- Volcano plot  
  → `results/figures/volcano_plot.pdf`  
- Heatmap (top 50 genes)  
  → `results/figures/heatmap_top50.pdf`  

### 4. Gene Filtering
- Selected significant genes:  
  → adjusted p-value (`padj < 0.05`)

### 5. Functional Enrichment
- Gene Ontology (GO) analysis  
- KEGG pathway analysis  
- GO barplot visualization  
  → `results/figures/go_barplot.pdf`  

---

## 📊 Results

### Differential Expression
Approximately **4000 genes** were identified as significantly differentially expressed (adjusted p-value < 0.05) between treated and untreated samples. This indicates a broad transcriptional response to dexamethasone.

### Visualization
- The **volcano plot** shows a substantial number of both upregulated and downregulated genes.
- The **heatmap** demonstrates clear clustering of samples by treatment condition, confirming a strong treatment-specific transcriptional signature.

### Functional Enrichment
- GO analysis revealed enrichment in:
  - immune response
  - inflammatory signaling
  - cellular stress processes  

- KEGG pathways highlighted:
  - signaling and regulatory mechanisms
  - pathways associated with glucocorticoid response  

---

## 🧠 Interpretation
Dexamethasone induces coordinated transcriptional changes affecting immune and inflammatory pathways. This is consistent with its known role as a **glucocorticoid with anti-inflammatory and immunomodulatory effects**.

The results validate the computational workflow and demonstrate how RNA-seq analysis can uncover biologically meaningful patterns.

---

## 📁 Project Structure

---

## ⚠️ Notes / Issues
- Fixed plotting device issue using `dev.off()`
- Simplified dependency stack (removed tidyverse/clusterProfiler)

---

## 🚀 Next Steps
- Refine biological interpretation
- Improve visualization quality
- Add **Python-based machine learning analysis**
- Apply pipeline to a **cancer RNA-seq dataset**

---

## 🔁 Reproducibility

Work in progress
