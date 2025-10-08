# Bulk RNA‑seq Expression Analysis Mini Project 
This repository presents a mini‑project in bulk RNA‑seq analysis, structured as a reproducible workflow. It illustrates the core steps of exploratory RNA‑seq data analysis using a small demonstration dataset, while providing a framework that can be readily scaled to real‑world datasets.

## What is Bulk RNA‑seq? 🧬
Bulk RNA sequencing (bulk RNA‑seq) is a high‑throughput technique used to measure gene expression across an entire population of cells.
- RNA molecules are extracted from a sample (e.g., tissue, blood, or cultured cells), converted into cDNA, and sequenced.
- The resulting reads are aligned to a reference genome or transcriptome and counted per gene.
- Unlike single‑cell RNA‑seq, which captures expression at the level of individual cells, bulk RNA‑seq provides an average expression profile across all cells in the sample.
- It is widely used to identify differentially expressed genes, study biological pathways, and compare conditions (e.g., Healthy vs Disease).

*This project uses a toy dataset to illustrate the logic of bulk RNA‑seq analysis in a simplified, reproducible way.*

## Overview
- Type: Bulk RNA‑seq (synthetic dataset)
- Samples: 5 (3 Healthy, 2 TB)
- Genes: 6 synthetic gene IDs
- Goal: Showcase the essentials of bulk RNA‑seq analysis in a modular, annotated, and reproducible format

## Repository Structure


| Folder / File                          | Description                          |
|----------------------------------------|--------------------------------------|
| `script_and_report/`                   | Contains R scripts and reports       |
| ├── `Workflow_RScript.R`               | Main R script for the analysis       |
| ├── `RNA-seq_Mini_Analysis_Report.Rmd` | R Markdown report (source)           |
| └── `RNA-seq_Mini_Analysis_Report.html`| Knitted HTML report                  |
|                                        |                                      |
| `results/`                             | Exported results and figures         |
| ├── `gene_summary.csv`                 | Unified summary table of results     |
| └── `library_sizes_barplot.png`        | Barplot of library sizes             |
|                                        |                                      |
| `README.md`                            | Project documentation                |

## Workflow Steps
- Dataset creation – build a mini dataset
- Inspection – check dimensions, IDs, and metadata integrity
- Per‑gene statistics – compute mean, SD, variability
- Per‑sample totals – library size assessment
- Condition means – Healthy vs TB averages
- Fold‑change proxy – TB/Healthy ratios with safe division
- Upregulated genes – identify TB‑responsive genes
- Filtering & pruning – remove low/noisy genes
- Reordering & renaming – tidy sample and gene IDs
- Visualization – barplot of library sizes
- Unified summary table – consolidate all metrics
- Export – save results into results/

## Usage
- Run the analysis script in R:
source("script_and_report/Workflow_RScript.R")

## Outputs:
- Console summaries at each step
- Barplot of library sizes (saved in results/)
- Exported gene_summary.csv (in results/)
- Full report in script_and_report/

## Significance of the Analysis
This mini-project illustrates the core logic of bulk RNA-sequencing analysis in a modular and reproducible manner.
- Bridges raw data to insight: from raw counts to interpretable summaries
- Highlights QC principles: library size checks, filtering, pruning
- Condition‑specific insights: identify TB‑responsive genes
- Reproducible structure: unified summary tables and exports
- Scalable logic: same workflow scales to real RNA‑seq datasets

## License
The project is released under the MIT License, permitting reuse, modification, and distribution with proper attribution.
