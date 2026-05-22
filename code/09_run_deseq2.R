suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

COUNT_MATRIX <- "results/htseq_combined/combined.tsv"
OUTDIR <- "results/deseq2"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "plots"), recursive = TRUE, showWarnings = FALSE)

raw <- read.table(COUNT_MATRIX,
                  header = TRUE,
                  sep = "\t",
                  row.names = 1,
                  check.names = FALSE,
                  comment.char = "")

gene_names <- data.frame(
  GeneID = rownames(raw),
  Name = raw[, 1],
  stringsAsFactors = FALSE
)
counts_raw <- raw[, -1]

counts <- counts_raw[!grepl("^__", rownames(counts_raw)),]

counts_mat <- as.matrix(counts)
mode(counts_mat) <- "integer"

cat("Count matrix dimensions:", nrow(counts_mat), "genes x",
    ncol(counts_mat), "samples\n")
cat("Samples:", paste(colnames(counts_mat), collapse = ", "), "\n")

sample_names <- colnames(counts_mat)

col_data <- data.frame(
  row.names = sample_names,
  condition = factor(
    ifelse(grepl("_BH$", sample_names), "BH", "Serum"),
    levels = c("BH", "Serum")
  ),
  sample = sub("_(BH|Serum)$", "", sample_names)
)

cat("Sample metadata:\n")
print(col_data)

keep <- rowSums(counts_mat) >= 10
cat("Genes before filter:", nrow(counts_mat), "\n")
counts_mat <- counts_mat[keep,]
gene_names <- gene_names[keep,]
cat("Genes after filter (rowSum >= 10):", nrow(counts_mat), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_data,
  design = ~condition
)

dds <- DESeq(dds)

print(sizeFactors(dds))

res <- results(dds,
               contrast = c("condition", "Serum", "BH"),
               alpha = 0.05)

res_df <- as.data.frame(res)
res_df <- cbind(
  GeneID = rownames(res_df),
  Name = gene_names$Name[match(rownames(res_df), gene_names$GeneID)],
  res_df
)
res_df <- res_df[order(res_df$padj, na.last = TRUE),]

res_df$significant <- "Not significant"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in Serum"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Up in BH"

sig_up <- sum(res_df$significant == "Up in Serum", na.rm = TRUE)
sig_down <- sum(res_df$significant == "Up in BH", na.rm = TRUE)
cat("DEGs (padj < 0.05, |log2FC| > 1):",
    sig_up, "up in Serum,", sig_down, "up in BH\n")

write.table(res_df,
            file = file.path(OUTDIR, "deseq2_results_all.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

sig_df <- res_df[!is.na(res_df$padj) &
                   res_df$padj < 0.05 &
                   abs(res_df$log2FoldChange) > 1,]
write.table(sig_df,
            file = file.path(OUTDIR, "deseq2_results_significant.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Significant DEGs saved:", nrow(sig_df), "genes\n")

vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("BH" = "#2166ac", "Serum" = "#d6604d")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA - BH vs Serum") +
  theme_bw(base_size = 13)

ggsave(file.path(OUTDIR, "pca_plot.pdf"),
       pca_plot, width = 7, height = 5)
ggsave(file.path(OUTDIR, "pca_plot.png"),
       pca_plot, width = 7, height = 5, dpi = 300)

plot_df <- res_df[!is.na(res_df$padj),]
plot_df$neg_log10_padj <- -log10(plot_df$padj)

top_genes <- head(plot_df[plot_df$significant != "Not significant",], 20)

volcano <- ggplot(plot_df,
                  aes(x = log2FoldChange,
                      y = neg_log10_padj,
                      color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text_repel(data = top_genes,
                  aes(label = ifelse(Name != GeneID, Name, GeneID)),
                  size = 2.8, max.overlaps = 30) +
  scale_color_manual(values = c(
    "Not significant" = "grey70",
    "Up in Serum" = "#d6604d",
    "Up in BH" = "#2166ac"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  xlab("log2 Fold Change (Serum / BH)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot - Serum vs BH",
          subtitle = paste0(sig_up, " up in Serum | ", sig_down, " up in BH")) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank())

ggsave(file.path(OUTDIR, "volcano_plot.pdf"),
       volcano, width = 8, height = 6)
ggsave(file.path(OUTDIR, "volcano_plot.png"),
       volcano, width = 8, height = 6, dpi = 300)

pdf(file.path(OUTDIR, "ma_plot.pdf"), width = 7, height = 5)
plotMA(res, ylim = c(-6, 6), main = "MA plot - Serum vs BH",
       alpha = 0.05, colSig = "#d6604d")
dev.off()

sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
anno_col <- data.frame(Condition = col_data$condition,
                       row.names = rownames(col_data))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(OUTDIR, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = anno_col,
         col = colors,
         main = "Sample-to-sample distances (VST)")
dev.off()

top50 <- head(sig_df, 50)

if (nrow(top50) >= 2) {
  mat <- assay(vsd)[top50$GeneID,]
  mat <- mat - rowMeans(mat)

  row_labels <- ifelse(top50$Name != top50$GeneID,
                       paste0(top50$Name, " (", top50$GeneID, ")"),
                       top50$GeneID)
  rownames(mat) <- row_labels

  pdf(file.path(OUTDIR, "top50_deg_heatmap.pdf"), width = 9, height = 12)
  pheatmap(mat,
           annotation_col = anno_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           fontsize_row = 7,
           color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
           main = "Top 50 DEGs centered VST counts")
  dev.off()
} else {
  cat("Fewer than 2 significant DEGs, skipping top50 heatmap\n")
}

pdf(file.path(OUTDIR, "dispersion_plot.pdf"), width = 7, height = 5)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

top20 <- head(sig_df, 20)

if (nrow(top20) >= 2) {
  mat <- assay(vsd)[top20$GeneID,]
  mat <- mat - rowMeans(mat)

  row_labels <- ifelse(top20$Name != top20$GeneID,
                       paste0(top20$Name, " (", top20$GeneID, ")"),
                       top20$GeneID)
  rownames(mat) <- row_labels

  pdf(file.path(OUTDIR, "top20_deg_heatmap.pdf"), width = 7, height = 6)
  pheatmap(mat,
           annotation_col = anno_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           fontsize_row = 8,
           color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
           main = "Top 20 DEGs centered VST counts")
  dev.off()

  write.csv(top20, file.path(OUTDIR, "top20.csv"), row.names = FALSE)
}