
# airway_rnaseq_pipeline.R
# End-to-end RNA-seq case study using Bioconductor `airway`
# Outputs: results/*.csv, results/*.txt, figures/*.png
R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('DESeq2','airway'), ask=FALSE, update=TRUE)"
suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
  library(ggplot2)
  library(pROC)
  library(randomForest)
  library(matrixStats)
})

set.seed(123)

# Ensure output directories
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ------------------------------------------------------------------
# Load data, set design, filter
# ------------------------------------------------------------------
data(airway)  # SummarizedExperiment with counts + colData

# Recode condition and set reference level
colData(airway)$dex <- factor(colData(airway)$dex, levels = c("untrt","trt"),
                              labels = c("untreated","treated"))

dds <- DESeqDataSet(airway, design = ~ dex)

# Filter low-count genes to reduce noise
dds <- dds[rowSums(counts(dds)) >= 10, ]

# ------------------------------------------------------------------
# Differential expression
# ------------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(res_df, "results/DESeq2_results.csv", row.names = FALSE)

# ------------------------------------------------------------------
# Volcano
# ------------------------------------------------------------------
res_df$neglog10padj <- -log10(res_df$padj)
res_df$signif <- res_df$padj < 0.05 & !is.na(res_df$padj)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = signif), alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  theme_minimal() +
  labs(title = "Volcano: treated vs untreated",
       x = "Log2 Fold Change", y = "-log10 adjusted p-value")
ggsave("figures/volcano.png", p_volcano, width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------------
# PCA on VST-normalised counts
# ------------------------------------------------------------------
rld <- vst(dds, blind = FALSE)
pca_data <- plotPCA(rld, intgroup = "dex", returnData = TRUE)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = dex)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA (VST)")
ggsave("figures/pca.png", p_pca, width = 7, height = 5, dpi = 300)

# ------------------------------------------------------------------
# Toy ML: RF predicts sample condition from top-variable VST genes
# (Demonstration only — tiny sample size)
# ------------------------------------------------------------------
mat <- assay(rld)  # genes x samples (VST)
vars <- rowVars(mat)
top_idx <- order(vars, decreasing = TRUE)[1:min(500, nrow(mat))]
X <- t(mat[top_idx, ])  # samples x genes
y <- factor(colData(dds)$dex, levels = c("untreated","treated"))

n <- nrow(X)

# guard: if n < 4, skip ML
do_ml <- n >= 4

acc <- NA_real_; auc_val <- NA_real_; cm <- NULL

if (do_ml) {
  # re-sample until both classes exist in train and test (limit tries)
  tries <- 0
  ok <- False <- FALSE
  while (tries < 50 && !ok) {
    idx <- sample(seq_len(n), round(0.7 * n))
    trainX <- X[idx, , drop = FALSE]
    testX  <- X[-idx, , drop = FALSE]
    trainY <- y[idx]
    testY  <- y[-idx]
    ok <- length(unique(trainY)) == 2 && length(unique(testY)) == 2
    tries <- tries + 1
  }

  if (ok) {
    rf <- randomForest(x = trainX, y = trainY,
                       ntree = 500,
                       mtry = max(1, round(sqrt(ncol(trainX)))),
                       importance = TRUE)

    pred_prob <- predict(rf, testX, type = "prob")[, "treated"]
    pred_class <- ifelse(pred_prob >= 0.5, "treated", "untreated")
    pred_class <- factor(pred_class, levels = levels(testY))

    acc <- mean(pred_class == testY)
    cm <- table(Predicted = pred_class, Actual = testY)

    # ROC only if both classes present
    if (length(unique(testY)) == 2) {
      roc_obj <- roc(response = testY, predictor = pred_prob, quiet = TRUE)
      auc_val <- as.numeric(auc(roc_obj))
      png("figures/roc.png", width = 800, height = 600)
      plot(roc_obj, main = sprintf("ROC (AUC = %.3f)", auc_val))
      dev.off()
    }
  }
}

# Save metrics + confusion matrix (if available)
con <- file("results/metrics.txt", open = "wt")
writeLines(sprintf("Accuracy: %s", ifelse(is.na(acc), "NA", sprintf("%.3f", acc))), con)
writeLines(sprintf("AUC: %s", ifelse(is.na(auc_val), "NA", sprintf("%.3f", auc_val))), con)
if (!is.null(cm)) {
  writeLines("\nConfusion Matrix:", con)
  capture.output(cm, file = con, append = TRUE)
}
close(con)

# Save session info for reproducibility
sink("results/sessionInfo.txt"); sessionInfo(); sink()
