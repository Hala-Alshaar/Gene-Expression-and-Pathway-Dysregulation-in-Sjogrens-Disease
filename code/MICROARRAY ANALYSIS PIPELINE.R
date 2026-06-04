################################################################
# COMPLETE MICROARRAY ANALYSIS PIPELINE
# GSE51092
################################################################

# =========================
# 1. LOAD LIBRARIES
# =========================
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

library(ggplot2)
library(ggrepel)

# =========================
# 2. DOWNLOAD GEO DATA
# =========================

gset <- getGEO("GSE51092",
               GSEMatrix = TRUE,
               AnnotGPL = TRUE)

if (length(gset) > 1) {
  idx <- grep("GPL6884", attr(gset, "names"))
} else {
  idx <- 1
}

gset <- gset[[idx]]

# clean feature names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# =========================
# 3. DEFINE GROUPS
# =========================

gsms <- paste0(
  "00000000000000000000000000000000111111111111111111",
  "11111111111111111111111111111111111111111111111111",
  "11111111111111111111111111111111111111111111111111",
  "11111111111111111111111111111111111111111111111111",
  "1111111111111111111111"
)

sml <- strsplit(gsms, split="")[[1]]

gs <- factor(sml)

groups <- c("control", "case")

levels(gs) <- groups

gset$group <- gs

# =========================
# 4. EXTRACT EXPRESSION
# =========================

ex <- exprs(gset)

# =========================
# 5. LOG2 TRANSFORMATION
# =========================

qx <- as.numeric(quantile(ex,
                          c(0,0.25,0.5,0.75,0.99,1),
                          na.rm=TRUE))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogC) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

exprs(gset) <- ex

# =========================
# 6. NORMALIZATION
# =========================

ex <- normalizeBetweenArrays(ex)

exprs(gset) <- ex

# =========================
# 7. REMOVE MISSING VALUES
# =========================

gset <- gset[complete.cases(exprs(gset)), ]

ex <- exprs(gset)

# =========================
# 8. QUALITY CONTROL
# =========================

# Boxplot
boxplot(ex,
        outline=FALSE,
        las=2,
        col=as.numeric(gs),
        main="Normalized Expression")

legend("topright",
       legend=levels(gs),
       fill=1:length(levels(gs)))

# Density plot
plotDensities(ex,
              group=gs,
              main="Expression Density")

# =========================
# 9. DESIGN MATRIX
# =========================

design <- model.matrix(~0 + group, gset)

colnames(design) <- levels(gs)

design

# =========================
# 10. LIMMA DIFFERENTIAL EXPRESSION
# =========================

fit <- lmFit(gset, design)

# CASE - CONTROL
cont.matrix <- makeContrasts(
  case-control,
  levels=design
)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

# =========================
# 11. RESULTS TABLE
# =========================

results <- topTable(fit2,
                    adjust="fdr",
                    number=Inf)

write.csv(results,
          "Differential_Expression_Results.csv",
          row.names=FALSE)

# =========================
# 12. SIGNIFICANT GENES
# =========================

sig <- subset(results,
              adj.P.Val < 0.05)

# Upregulated in CASE
upregulated <- subset(sig,
                      logFC > 1)

# Downregulated in CASE
downregulated <- subset(sig,
                        logFC < -1)

write.csv(upregulated,
          "Upregulated_Case.csv",
          row.names=FALSE)

write.csv(downregulated,
          "Downregulated_Case.csv",
          row.names=FALSE)

# =========================
# 13. VOLCANO PLOT
# =========================

results$threshold <- as.factor(
  abs(results$logFC) > 1 &
    results$adj.P.Val < 0.05
)

ggplot(results,
       aes(x=logFC,
           y=-log10(adj.P.Val),
           color=threshold)) +
  geom_point(alpha=0.7) +
  geom_text_repel(
    data=head(results[order(results$adj.P.Val), ], 15),
    aes(label=Gene.symbol),
    size=3
  ) +
  theme_minimal() +
  labs(title="Volcano Plot")
#another way to plot it 
library(ggrepel)

# Select top 15 most significant genes
topgenes <- head(
  results[order(abs(results$logFC), decreasing = TRUE), ],
  20
)

ggplot(results,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = status)) +
  
  geom_point(alpha = 0.8,
             size = 2) +
  
  scale_color_manual(values = c(
    "Downregulated" = "blue",
    "Not Significant" = "grey70",
    "Upregulated" = "red"
  )) +
  
  # Add gene labels
  geom_text_repel(
    data = topgenes,
    aes(label = Gene.symbol),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20
  ) +
  
  theme_minimal() +
  
  labs(title = "Differential Expression Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
##
library(ggrepel)

# Top 10 upregulated genes
top_up <- results[
  results$logFC > 1 &
    results$adj.P.Val < 0.05,
]

top_up <- head(
  top_up[order(top_up$logFC, decreasing = TRUE), ],
  43
)

# Top 10 downregulated genes
top_down <- results[
  results$logFC < -1 &
    results$adj.P.Val < 0.05,
]

top_down <- head(
  top_down[order(top_down$logFC), ],
  10
)

# Combine both
topgenes <- rbind(top_up, top_down)

# Volcano plot
ggplot(results,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = status)) +
  
  geom_point(alpha = 0.8,
             size = 2) +
  
  scale_color_manual(values = c(
    "Downregulated" = "blue",
    "Not Significant" = "grey70",
    "Upregulated" = "red"
  )) +
  
  # Threshold lines
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  # Labels
  geom_text_repel(
    data = topgenes,
    aes(label = Gene.symbol),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 30
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Differential Expression Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
# =========================
# 15. PCA
# =========================

pca <- prcomp(t(ex),
              scale.=TRUE)

plot(pca$x[,1],
     pca$x[,2],
     col=as.numeric(gs),
     pch=20,
     cex=1.5,
     xlab="PC1",
     ylab="PC2",
     main="PCA")

legend("topright",
       legend=levels(gs),
       col=1:length(levels(gs)),
       pch=20)

# =========================
# 16. UMAP
# =========================

ex2 <- na.omit(ex)

ex2 <- ex2[!duplicated(ex2), ]

ump <- umap(t(ex2),
            n_neighbors=15,
            random_state=123)

plot(ump$layout,
     col=as.numeric(gs),
     pch=20,
     cex=1.5,
     main="UMAP")

legend("topright",
       legend=levels(gs),
       col=1:length(levels(gs)),
       pch=20)

# =========================
# 17. K-MEANS CLUSTERING
# =========================

k <- 2

kmeans_result <- kmeans(t(ex2),
                        centers=k)

table(kmeans_result$cluster)

plot(pca$x[,1],
     pca$x[,2],
     col=kmeans_result$cluster,
     pch=20,
     cex=1.5,
     main="PCA + Kmeans")

legend("topright",
       legend=paste("Cluster", 1:k),
       col=1:k,
       pch=20)

# =========================
# 18. GO ENRICHMENT
# =========================

genes_up <- unique(upregulated$Gene.symbol)

gene.df <- bitr(
  genes_up,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)

gene_ids <- gene.df$ENTREZID

ego <- enrichGO(
  gene          = gene_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(ego),
          "GO_Enrichment.csv",
          row.names=FALSE)

dotplot(ego,
        showCategory=10)

# =========================
# 19. KEGG ENRICHMENT
# =========================

ekegg <- enrichKEGG(
  gene         = gene_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(ekegg),
          "KEGG_Enrichment.csv",
          row.names=FALSE)

dotplot(ekegg,
        showCategory=20)

# =========================
# 20. REACTOME ENRICHMENT
# =========================

react <- enrichPathway(
  gene=gene_ids,
  organism="human"
)

write.csv(as.data.frame(react),
          "Reactome_Enrichment.csv",
          row.names=FALSE)

dotplot(react,
        showCategory=15)

# =========================
# 21. GSEA
# =========================

gene_ranks <- results$logFC

names(gene_ranks) <- results$Gene.symbol

gene_ranks <- sort(gene_ranks,
                   decreasing=TRUE)

gsea <- gseGO(
  geneList     = gene_ranks,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "SYMBOL",
  verbose      = FALSE
)

ridgeplot(gsea,  showCategory=10)

gseaplot2(gsea,
          geneSetID=1)

# =========================
# 22. MEAN VARIANCE TREND
# =========================

plotSA(fit2,
       main="Mean-Variance Trend")

# =========================
# 23. SAVE COMPLETE SESSION
# =========================

save.image("Complete_GSE51092_Analysis.RData")

################################################################
# END OF PIPELINE
################################################################
