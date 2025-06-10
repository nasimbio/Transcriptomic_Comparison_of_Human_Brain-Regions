# Bulk RNA-seq Analysis — Brain Region Comparison

This project focuses on **bulk RNA-seq analysis** comparing two distinct brain regions: **Hippocampus (HP)** and **Motor Cortex (MC)**. Each region includes two experimental groups, labeled **A** and **B**. The objective is to identify differentially expressed genes (DEGs) and enriched biological processes that distinguish conditions both *within* and *across* brain regions.

---

## 📘 Project Overview

The goal of this project is to investigate transcriptional differences within and between the hippocampus (HP) and motor cortex (MC) brain regions by analyzing bulk RNA-seq data. Each brain region is subdivided into groups A and B. Differential expression and Gene Ontology enrichment analyses were conducted to explore region-specific biological variation and experimental effects.

**Data Preprocessing**
- Raw FASTQ files were quality checked and aligned to the reference human genome
- Gene count matrices were generated from aligned read

---

## 📊 Analysis Tasks


### **Task: Differential Gene Expression (DGE) and GO Enrichment**
- **Within-region comparisons:**
  - DGE and GO analysis of **HP A vs. HP B**
  - DGE and GO analysis of **MC A vs. MC B**
- **Across-region comparisons (same group):**
  - DGE and GO analysis of **HP A vs. MC A**
  - DGE and GO analysis of **HP B vs. MC B**
- DEG analysis was performed using standard pipelines (DESeq2)
- Significant DEGs were filtered by log2 fold change and adjusted p-value thresholds
- Gene Ontology (GO) enrichment was performed using the `clusterProfiler` package

---

# 💡 Notes

- Raw data is not publicly available due to client ownership and confidentiality.
- Some example outputs plots are organized by task in the `output/` folder.
- This project is designed for both reproducibility and clarity.

---

## 📬 Contact

*Author:* Nasim Rahmatpour 
*Email:* nasimrahmatpour1@gmail.com 
*GitHub:* (https://github.com/nasimbio)

