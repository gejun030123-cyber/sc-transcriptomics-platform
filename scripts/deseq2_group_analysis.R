#!/usr/bin/env Rscript

# Eight-group DESeq2 workflow
#
# DEG inference is fitted once across all treatment x compartment groups.
# B and En PCA plots are transformed separately for within-compartment viewing
# only; this does not alter the DEG model or its contrasts.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

counts_path <- "counts/raw_counts_clean.tsv"
metadata_path <- "metadata.tsv"
output_dir <- "deg"
plot_dir <- file.path(output_dir, "plots")

treatment_levels <- c("Ctr", "NH4Cl", "PEA", "TMAO")
compartment_levels <- c("B", "En")
group_levels <- as.vector(outer(treatment_levels, compartment_levels,
  FUN = function(treatment, compartment) paste(treatment, compartment, sep = "_")
))

# =========================
# 1. Read and validate input
# =========================

counts <- read.delim(
  counts_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
metadata <- read.delim(
  metadata_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (anyDuplicated(colnames(counts))) {
  stop("Duplicate sample IDs in count-matrix columns.")
}
if (anyDuplicated(rownames(metadata))) {
  stop("Duplicate sample IDs in metadata row names.")
}
if (!setequal(colnames(counts), rownames(metadata))) {
  missing_metadata <- setdiff(colnames(counts), rownames(metadata))
  missing_counts <- setdiff(rownames(metadata), colnames(counts))
  stop(
    "Count matrix and metadata do not contain the same sample IDs. ",
    "Missing metadata: ", paste(missing_metadata, collapse = ", "), "; ",
    "missing count columns: ", paste(missing_counts, collapse = ", ")
  )
}
if (!all(c("treatment", "compartment") %in% colnames(metadata))) {
  stop("metadata.tsv must contain treatment and compartment columns.")
}
if (!all(vapply(counts, is.numeric, logical(1)))) {
  stop("Count matrix contains non-numeric values.")
}

counts <- as.matrix(counts)
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
  stop("DESeq2 requires a complete matrix of non-negative integer raw counts.")
}

# Reorder only after checking the two sample-ID sets are identical.
metadata <- metadata[colnames(counts), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))

metadata$treatment <- factor(as.character(metadata$treatment), levels = treatment_levels)
metadata$compartment <- factor(as.character(metadata$compartment), levels = compartment_levels)
if (anyNA(metadata$treatment) || anyNA(metadata$compartment)) {
  stop(
    "Unexpected treatment or compartment value. Expected treatments: ",
    paste(treatment_levels, collapse = ", "), "; expected compartments: ",
    paste(compartment_levels, collapse = ", ")
  )
}
metadata$group <- factor(
  paste(metadata$treatment, metadata$compartment, sep = "_"),
  levels = group_levels
)

group_sizes <- table(metadata$group)
if (any(group_sizes == 0)) {
  stop("Every group needed for the six requested contrasts must have samples: ",
       paste(names(group_sizes)[group_sizes == 0], collapse = ", "))
}
if (any(group_sizes < 2)) {
  warning("At least one group has fewer than two biological replicates: ",
          paste(names(group_sizes)[group_sizes < 2], collapse = ", "))
}
print(group_sizes)

# =========================
# 2. Global low-expression filter
# =========================

keep <- rowSums(counts >= 10) >= 3
counts_filt <- counts[keep, , drop = FALSE]
if (nrow(counts_filt) == 0) {
  stop("No genes remain after requiring count >= 10 in at least 3 samples.")
}
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering :", nrow(counts_filt), "\n")

# =========================
# 3. One global DEG fit
# =========================

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = metadata,
  design = ~ group
)
dds <- DESeq(dds)
vsd_all <- vst(dds, blind = FALSE)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 4. Within-compartment PCA
# =========================

save_compartment_pca <- function(counts, metadata, compartment, output_file) {
  sample_keep <- metadata$compartment == compartment
  metadata_sub <- droplevels(metadata[sample_keep, , drop = FALSE])
  counts_sub <- counts[, sample_keep, drop = FALSE]

  # This filter and fit are intentionally local to the PCA visualization.
  # The DEG model above remains the single global fit.
  gene_keep <- rowSums(counts_sub >= 10) >= 3
  counts_sub <- counts_sub[gene_keep, , drop = FALSE]
  if (nrow(counts_sub) == 0) {
    warning("No genes remain for ", compartment, " PCA; skipping plot.")
    return(invisible(NULL))
  }

  dds_sub <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = metadata_sub,
    design = ~ treatment
  )
  dds_sub <- DESeq(dds_sub, quiet = TRUE)
  vsd_sub <- vst(dds_sub, blind = FALSE)

  pca_data <- plotPCA(
    vsd_sub,
    intgroup = "treatment",
    ntop = min(5000, nrow(vsd_sub)),
    returnData = TRUE
  )
  percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

  p <- ggplot(pca_data, aes(PC1, PC2, color = treatment, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.7, size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(paste(compartment, "samples")) +
    theme_bw()

  ggsave(output_file, p, width = 8, height = 6)
}

save_compartment_pca(counts_filt, metadata, "B", file.path(plot_dir, "PCA_B_only.pdf"))
save_compartment_pca(counts_filt, metadata, "En", file.path(plot_dir, "PCA_En_only.pdf"))

# =========================
# 5. Pairwise DEG function
# =========================

run_deg <- function(dds, vsd, numerator, denominator, prefix) {
  res <- results(
    dds,
    contrast = c("group", numerator, denominator),
    alpha = 0.05
  )

  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]

  res_df$status <- "NS"
  is_significant <- !is.na(res_df$padj) & res_df$padj < 0.05
  res_df$status[is_significant & res_df$log2FoldChange >= 1] <- "Up"
  res_df$status[is_significant & res_df$log2FoldChange <= -1] <- "Down"

  write.table(
    res_df,
    file = file.path(output_dir, paste0(prefix, "_all.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  sig <- res_df[res_df$status %in% c("Up", "Down"), , drop = FALSE]
  write.table(
    sig,
    file = file.path(output_dir, paste0(prefix, "_DEG.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  up_n <- sum(res_df$status == "Up")
  down_n <- sum(res_df$status == "Down")
  cat(prefix, "Up:", up_n, "Down:", down_n, "\n")

  plot_df <- res_df
  plot_df$minuslog10padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = minuslog10padj)) +
    geom_point(aes(shape = status), alpha = 0.6, size = 1.5, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted P value") +
    ggtitle(prefix) +
    theme_bw()
  ggsave(file.path(plot_dir, paste0(prefix, "_volcano.pdf")), p_volcano, width = 7, height = 6)

  pdf(file.path(plot_dir, paste0(prefix, "_MA.pdf")), width = 7, height = 6)
  plotMA(res, alpha = 0.05, ylim = c(-6, 6))
  dev.off()

  # Show only the two groups in this contrast. Including all eight groups here
  # can make a pairwise DEG heatmap look inconsistent with its contrast.
  if (nrow(sig) >= 2) {
    top_genes <- head(sig$Geneid, 50)
    sample_keep <- colData(dds)$group %in% c(numerator, denominator)
    mat <- assay(vsd)[top_genes, sample_keep, drop = FALSE]
    mat <- mat - rowMeans(mat)
    annotation_col <- as.data.frame(
      colData(dds)[sample_keep, c("treatment", "compartment", "group"), drop = FALSE]
    )
    annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

    pdf(file.path(plot_dir, paste0(prefix, "_top50_heatmap.pdf")), width = 10, height = 10)
    pheatmap(
      mat,
      annotation_col = annotation_col,
      show_rownames = TRUE,
      fontsize_row = 6
    )
    dev.off()
  }

  data.frame(
    comparison = prefix,
    up = up_n,
    down = down_n,
    total = up_n + down_n
  )
}

# =========================
# 6. Six pre-specified contrasts
# Positive log2FC means higher expression in the numerator group.
# =========================

comparisons <- data.frame(
  numerator = c("NH4Cl_B", "PEA_B", "TMAO_B", "NH4Cl_En", "PEA_En", "TMAO_En"),
  denominator = c("Ctr_B", "Ctr_B", "Ctr_B", "Ctr_En", "Ctr_En", "Ctr_En"),
  prefix = c(
    "NH4Cl_B_vs_Ctr_B", "PEA_B_vs_Ctr_B", "TMAO_B_vs_Ctr_B",
    "NH4Cl_En_vs_Ctr_En", "PEA_En_vs_Ctr_En", "TMAO_En_vs_Ctr_En"
  ),
  stringsAsFactors = FALSE
)

summary_list <- lapply(seq_len(nrow(comparisons)), function(i) {
  run_deg(
    dds = dds,
    vsd = vsd_all,
    numerator = comparisons$numerator[i],
    denominator = comparisons$denominator[i],
    prefix = comparisons$prefix[i]
  )
})
deg_summary <- do.call(rbind, summary_list)
write.table(
  deg_summary,
  file.path(output_dir, "DEG_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
print(deg_summary)
cat("\nAll DEG analyses finished.\n")
