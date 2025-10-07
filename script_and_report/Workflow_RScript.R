#-------------------------------------------------------------------------------
# Exploratory RNA-seq Mini Analysis Workflow
# Purpose: Demonstrate stepwise RNA-seq style analysis
# Dataset: 6 genes × 5 samples (3 Healthy, 2 TB)
#-------------------------------------------------------------------------------

set.seed(123)

# 0) Create dataset 
# -----------------------------
counts <- matrix(
  c(
    50,  55,  52, 120, 130,
    300, 280, 310, 200, 190,
    15,  14,  16,  17,  15,
    400, 420, 410, 800, 820,
    150, 160, 155, 120, 110,
    8,   7,   9,   8,   7
  ),
  nrow = 6, byrow = TRUE
)

rownames(counts) <- paste0("G", 1:6)
colnames(counts) <- paste0("S", 1:5)

samples <- data.frame(
  SampleID  = paste0("S", 1:5),
  Condition = c("Healthy", "Healthy", "Healthy", "TB", "TB"),
  stringsAsFactors = FALSE
)

# 1) Inspect dataset structure 
# -----------------------------
dim(counts)         # dimensions (genes x samples)
rownames(counts)    # gene IDs
colnames(counts)    # sample IDs
str(samples)        # metadata structure
samples             # full metadata table

# 2) Per-gene mean & SD; top 2 most variable  
# -----------------------------
# Calculate mean and standard deviation for each gene across all samples
gene_mean <- apply(counts, 1, mean)
gene_sd   <- apply(counts, 1, sd)

# Combine results into a summary table
gene_stats <- data.frame(
  Gene = rownames(counts),
  Mean = round(gene_mean, 2),
  SD   = round(gene_sd, 2)
)
gene_stats   # view per-gene mean and SD

# Identify the top 2 most variable genes (highest SD)
top2_var <- head(gene_stats[order(-gene_stats$SD), ], 2)
top2_var

# 3) Per-sample totals (library size)  
# -----------------------------
# Compute total counts per sample (column sums)
sample_totals <- colSums(counts)
sample_totals

# Identify which sample has the highest total (largest library size)
names(which.max(sample_totals))

# 4) Gene x Condition mean  
# -----------------------------
# Compute mean expression per gene within each condition (Healthy vs TB)
gene_cond_means <- tapply(
  counts,
  list(
    rep(rownames(counts), ncol(counts)),
    rep(samples$Condition, each = nrow(counts))
  ),
  mean
)
gene_cond_means

# 5) Genes upregulated in TB  
# -----------------------------
# Select genes where TB mean > Healthy mean
up_in_TB <- rownames(counts)[gene_cond_means[, "TB"] > gene_cond_means[, "Healthy"]]
up_in_TB

# 6) Largest fold change proxy (TB/Healthy)  
# -----------------------------
# Safe division to avoid division by zero
safe_div <- function(tb, healthy) {
  if (healthy == 0) return(if (tb > 0) Inf else NA)
  tb / healthy
}

# Compute fold-change proxy for each gene
safe_fc <- mapply(
  safe_div,
  tb = gene_cond_means[, "TB"],
  healthy = gene_cond_means[, "Healthy"]
)

# Preserve gene names
names(safe_fc) <- rownames(gene_cond_means)

# Round for readability
safe_fc <- round(safe_fc, 3)

# Subset to only upregulated genes
fc_up <- safe_fc[up_in_TB]

# Identify gene with largest fold-change
if (length(fc_up)) {
  top_fc_gene  <- names(which.max(fc_up))
  top_fc_value <- max(fc_up)
} else {
  top_fc_gene  <- NA
  top_fc_value <- NA
}

safe_fc
fc_up
top_fc_gene
top_fc_value

# 7) Label samples HighLoad vs Normal  
# -----------------------------
# Define HighLoad as samples with total counts above the median
median_total <- median(sample_totals)
LoadLabel <- ifelse(sample_totals > median_total, "HighLoad", "Normal")

# Add total counts and load label to the sample metadata
samples$Total     <- sample_totals[samples$SampleID]
samples$LoadLabel <- LoadLabel[samples$SampleID]
samples

# 8) Table(condition, HighLoad)  
# -----------------------------
# Cross-tabulate condition vs load label
tab_cond_load <- table(samples$Condition, samples$LoadLabel)
tab_cond_load

# 9) Row-standardized counts (subtract row mean)  
# -----------------------------
# Center each gene’s expression by subtracting its mean
# This highlights relative up/down shifts across samples
row_centered <- counts - gene_mean
row_centered

# 10) Build list rna_seq; extract TB mean for G3  
# -----------------------------
# Store counts, samples, and condition means in a list
rna_seq <- list(
  counts = counts,
  samples = samples,
  gene_cond_means = gene_cond_means
)

# Extract the TB mean for gene G3 using nested indexing
rna_seq$gene_cond_means["G3", "TB"]

# 11) Gene with highest single value 
# -----------------------------
max_val <- max(counts)
pos <- which(counts == max_val, arr.ind = TRUE)

rownames(counts)[pos[1, "row"]]   # gene with max
colnames(counts)[pos[1, "col"]]   # sample with max
max_val

# 12) Per-condition total counts 
# -----------------------------
# Sum total counts across all samples within each condition (Healthy vs TB)
cond_totals <- tapply(colSums(counts), samples$Condition, sum)

# View the total counts per condition
cond_totals

# Identify which condition has the higher overall total expression
names(which.max(cond_totals))

# 13) Low-count filtering (<5 -> 0)
# -----------------------------

# Replace all counts < 5 with 0 to simulate low-count filtering.
# Count how many entries were affected.
counts_filtered <- counts
changed <- sum(counts_filtered < 5)
counts_filtered[counts_filtered < 5] <- 0

changed          # number of entries changed
counts_filtered  # filtered matrix

# 14) Recompute per-gene mean after filtering
# -----------------------------
# Compare top-mean gene before vs after filtering.
gene_mean_after <- apply(counts_filtered, 1, mean)
round(gene_mean_after, 3)

names(which.max(gene_mean))        # top mean gene before filtering
names(which.max(gene_mean_after))  # top mean gene after filtering

# 15) Reorder columns by condition (Healthy first, then TB)
# -----------------------------
# Create an ordering factor for samples by condition.
order_by_cond <- order(factor(samples$Condition, levels = c("Healthy", "TB")))
ordered_samples <- samples$SampleID[order_by_cond]

# Reorder the counts matrix accordingly.
counts_reordered <- counts[, ordered_samples, drop = FALSE]
colnames(counts_reordered)

# 16) Gene max–min range; most dynamic gene
# -----------------------------
# Compute the range (max - min) of expression for each gene.
gene_range <- apply(counts, 1, function(x) max(x) - min(x))
gene_range

# Identify the gene with the largest dynamic range.
names(which.max(gene_range))

# 17) Barplot of per-sample totals
# -----------------------------
# Visualize library sizes per sample.
op <- par(mar = c(5,4,2,1))
barplot(sample_totals,
        main = "Library sizes per sample",
        ylab = "Total counts",
        col = c(rep("steelblue", 3), rep("tomato", 2)))
par(op)

# 18) Rename gene G2 to G2A
# -----------------------------
# Update rownames to reflect new gene ID.
rownames(counts)[rownames(counts) == "G2"] <- "G2A"
rownames(counts)

# 19) Remove lowest-mean gene
# -----------------------------
# Identify the lowest-mean gene and remove it from the counts matrix.
lowest_mean_gene <- names(which.min(gene_mean))
counts_pruned <- counts[setdiff(rownames(counts), lowest_mean_gene), , drop = FALSE]

lowest_mean_gene  # removed gene
counts_pruned     # pruned matrix

# Justification:
# Very low-expression genes often add noise, reduce statistical power,
# and inflate multiple-testing burden. Removing them can improve
# downstream normalization and differential expression analysis.

# 20) Conclusion
# -----------------------------
# Build unified summary table
gene_summary <- data.frame(
  Gene = rownames(counts),
  Mean = round(gene_mean, 2),
  SD = round(gene_sd, 2),
  HealthyMean = round(gene_cond_means[, "Healthy"], 2),
  TBMean = round(gene_cond_means[, "TB"], 2),
  FoldChange = round(safe_fc, 3),
  Range = gene_range,
  UpregulatedInTB = rownames(counts) %in% up_in_TB,
  row.names = NULL   # <- prevents genes from appearing twice
)

# View the full table
gene_summary

# Save the summary table to CSV
write.csv(gene_summary, file = "gene_summary.csv", row.names = FALSE)


