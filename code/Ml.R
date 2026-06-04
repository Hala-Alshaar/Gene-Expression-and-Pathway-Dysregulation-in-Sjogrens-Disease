################################################################
# ADVANCED MACHINE LEARNING + PATHWAY ANALYSIS PIPELINE
# Compatible with GSE51092 LIMMA PIPELINE
################################################################

# ===============================================================
# 24. ADDITIONAL LIBRARIES
# ===============================================================
install.packages("caret", dependencies = c("Depends", "Suggests"))
library(randomForest)
library(MASS)
library(glmnet)

library(gprofiler2)
library(pathfindR)

library(caret)

# ===============================================================
# 25. PREPARE DATA FOR MACHINE LEARNING
# ===============================================================

# transpose expression matrix
# rows = samples
# columns = genes

ml_data <- as.data.frame(t(ex))

# add phenotype
ml_data$group <- gs

# remove problematic columns
ml_data <- ml_data[, colSums(is.na(ml_data)) == 0]

# ===============================================================
# 26. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ===============================================================

pca_ml <- prcomp(
  ml_data[, -ncol(ml_data)],
  scale. = TRUE
)

# variance explained
pca_var <- (pca_ml$sdev^2 / sum(pca_ml$sdev^2)) * 100

# PCA dataframe
pca_df <- data.frame(
  PC1 = pca_ml$x[,1],
  PC2 = pca_ml$x[,2],
  Group = ml_data$group
)

# PCA plot
ggplot(pca_df,
       aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1: ", round(pca_var[1],2), "%"),
    y = paste0("PC2: ", round(pca_var[2],2), "%")
  )

# scree plot
plot(
  pca_var[1:20],
  type="b",
  pch=19,
  xlab="Principal Component",
  ylab="Variance Explained (%)",
  main="Scree Plot"
)

# ===============================================================
# 27. K-MEANS CLUSTERING
# ===============================================================

set.seed(123)

# use PCA coordinates
kmeans_input <- pca_ml$x[,1:2]

# choose number of clusters
k <- 2

kmeans_res <- kmeans(
  kmeans_input,
  centers = k,
  nstart = 25
)

# cluster assignments
clusters <- as.factor(kmeans_res$cluster)

# compare clusters with true groups
table(Cluster=clusters,
      Group=gs)

# plotting
kmeans_df <- data.frame(
  PC1 = kmeans_input[,1],
  PC2 = kmeans_input[,2],
  Cluster = clusters,
  Group = gs
)

ggplot(kmeans_df,
       aes(PC1, PC2,
           color=Cluster,
           shape=Group)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="K-Means Clustering on PCA Space")

# ===============================================================
# 28. RANDOM FOREST CLASSIFICATION
# ===============================================================

set.seed(123)

# training data
x_rf <- ml_data[, -ncol(ml_data)]
y_rf <- ml_data$group

# train random forest
rf_model <- randomForest(
  x = x_rf,
  y = y_rf,
  ntree = 500,
  importance = TRUE
)

print(rf_model)

# predictions
rf_pred <- predict(rf_model)

# confusion matrix
confusionMatrix(rf_pred, y_rf)

# variable importance
importance_df <- importance(rf_model)

importance_df <- data.frame(
  Gene = rownames(importance_df),
  MeanDecreaseGini = importance_df[, "MeanDecreaseGini"]
)

importance_df <- importance_df[
  order(importance_df$MeanDecreaseGini,
        decreasing=TRUE),
]

# top genes
top_rf_genes <- head(importance_df, 20)

# plot importance
ggplot(top_rf_genes,
       aes(reorder(Gene, MeanDecreaseGini),
           MeanDecreaseGini)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title="Top Random Forest Important Genes",
    x="Gene",
    y="Importance"
  )

write.csv(
  importance_df,
  "RandomForest_FeatureImportance.csv",
  row.names=FALSE
)
##
# Variable importance
importance_df <- importance(rf_model)

importance_df <- data.frame(
  ID = rownames(importance_df),
  MeanDecreaseGini = importance_df[, "MeanDecreaseGini"]
)

# annotation table
annot <- fData(gset)

# map IDs to gene symbols
id_to_symbol <- annot[, c("ID", "Gene.symbol")]

importance_df <- merge(
  importance_df,
  id_to_symbol,
  by = "ID",
  all.x = TRUE
)

# remove missing symbols
importance_df <- importance_df[
  importance_df$Gene.symbol != "",
]

# sort
importance_df <- importance_df[
  order(importance_df$MeanDecreaseGini,
        decreasing = TRUE),
]

# top genes
top_rf_genes <- head(importance_df, 20)

# plot
ggplot(top_rf_genes,
       aes(reorder(Gene.symbol,
                   MeanDecreaseGini),
           MeanDecreaseGini)) +
  
  geom_bar(stat = "identity",
           fill = "steelblue") +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Top Random Forest Important Genes",
    x = "Gene Symbol",
    y = "Importance"
  )

#######  RANDOM FOREST CLASSIFICATION
library(randomForest)
library(caret)
library(ggplot2)
library(pheatmap)

set.seed(123)

# Training data
x_rf <- ml_data[, -ncol(ml_data)]
y_rf <- as.factor(ml_data$group)

# Random Forest with 1800 trees
rf_model <- randomForest(
  x = x_rf,
  y = y_rf,
  ntree = 1800,
  importance = TRUE,
  keep.forest = TRUE
)

print(rf_model)

# OOB prediction
rf_pred <- predict(rf_model)

# Confusion matrix
confusionMatrix(rf_pred, y_rf)
# Extract OOB error
oob_df <- data.frame(
  Trees = 1:nrow(rf_model$err.rate),
  OOB_Error = rf_model$err.rate[, "OOB"] * 100
)

ggplot(oob_df,
       aes(x = Trees,
           y = OOB_Error)) +
  geom_line(color = "steelblue",
            linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Out-of-Bag (OOB) Error Rate Across Trees",
    x = "Number of Trees",
    y = "OOB Error (%)"
  )
final_oob <- tail(oob_df$OOB_Error, 1)
final_oob

# ===============================================================
# 29. LINEAR DISCRIMINANT ANALYSIS (LDA)
# ===============================================================

# use top variable genes for stability
top_var_genes <- names(
  sort(apply(x_rf, 2, var),
       decreasing=TRUE)
)[1:100]

lda_data <- ml_data[, c(top_var_genes, "group")]

# fit LDA
lda_model <- lda(group ~ ., data=lda_data)

print(lda_model)

# predictions
lda_pred <- predict(lda_model)

# confusion matrix
confusionMatrix(
  as.factor(lda_pred$class),
  lda_data$group
)

# LDA plot
lda_df <- data.frame(
  LD1 = lda_pred$x[,1],
  Group = lda_data$group
)

ggplot(lda_df,
       aes(LD1,
           fill=Group)) +
  geom_density(alpha=0.5) +
  theme_minimal() +
  labs(title="LDA Discriminant Function")
# LDA coefficients
lda_coef <- lda_model$scaling[,1]

# select top contributing genes (adjust number if needed)
top_genes <- names(sort(abs(lda_coef), decreasing = TRUE))[1:30]

# subset expression data
heat_data <- ml_data[, top_genes]

# add group label
annotation <- data.frame(Group = lda_data$group)
rownames(annotation) <- rownames(lda_data)
library(pheatmap)

library(pheatmap)

# ===============================
# 1. Gene annotation (probe → symbol)
# ===============================
annot <- fData(gset)

id_to_symbol <- annot[, c("ID", "Gene.symbol")]
id_to_symbol <- id_to_symbol[
  !is.na(id_to_symbol$Gene.symbol) &
    id_to_symbol$Gene.symbol != "",
]

gene_map <- id_to_symbol$Gene.symbol
names(gene_map) <- id_to_symbol$ID

# ===============================
# 2. Select top LDA genes
# ===============================
lda_coef <- lda_model$scaling[,1]

top_genes_ids <- names(sort(abs(lda_coef), decreasing = TRUE))[1:30]

# map to gene symbols
top_genes_symbols <- gene_map[top_genes_ids]

# remove missing + make unique
keep <- !is.na(top_genes_symbols)
top_genes_ids <- top_genes_ids[keep]
top_genes_symbols <- make.unique(top_genes_symbols[keep])

# ===============================
# 3. Build expression matrix
# ===============================
heat_data <- ml_data[, top_genes_ids]

colnames(heat_data) <- top_genes_symbols

# ===============================
# 4. Sample annotation
# ===============================
annotation_col <- data.frame(
  Group = ml_data$group
)
rownames(annotation_col) <- rownames(ml_data)

# ===============================
# 5. Heatmap (Z-score scaled)
# ===============================
pheatmap(
  t(scale(heat_data)),
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  clustering_method = "complete",
  scale = "row",
  main = "Top LDA Genes Heatmap"
)
# order samples by group
ord <- order(ml_data$group)

heat_data_ordered <- heat_data[ord, ]

annotation_col <- data.frame(Group = ml_data$group[ord])
rownames(annotation_col) <- rownames(ml_data)[ord]

pheatmap(
  t(scale(heat_data_ordered)),
  annotation_col = annotation_col,
  cluster_cols = FALSE,   # IMPORTANT: disable clustering
  cluster_rows = TRUE,
  show_colnames = FALSE,
  scale = "row",
  main = "Top Genes (Samples Ordered by Group)"
)
# ===============================================================
# 30. LASSO LOGISTIC REGRESSION
# ===============================================================
set.seed(123)

train_index <- sample(seq_len(nrow(ml_data)), 0.7 * nrow(ml_data))

train_data <- ml_data[train_index, ]
test_data  <- ml_data[-train_index, ]

x_train <- as.matrix(train_data[, -ncol(train_data)])
y_train <- ifelse(train_data$group == "case", 1, 0)

x_test <- as.matrix(test_data[, -ncol(test_data)])
y_test <- ifelse(test_data$group == "case", 1, 0)
library(glmnet)

set.seed(123)

cvfit <- cv.glmnet(
  x_train,
  y_train,
  family = "binomial",
  alpha = 1,
  nfolds = 10
)

best_lambda <- cvfit$lambda.min
best_lambda
prob_test <- predict(
  cvfit,
  s = best_lambda,
  newx = x_test,
  type = "response"
)
library(pROC)

roc_obj <- roc(y_test, as.vector(prob_test))

plot(roc_obj, col="blue", main="LASSO Logistic Regression ROC")

auc_value <- auc(roc_obj)
auc_value
pred_class <- ifelse(prob_test > 0.5, 1, 0)

table(Predicted = pred_class, Actual = y_test)

mean(pred_class == y_test)  # accuracy
##### compare the modeks 
# ===============================================================
# MODEL COMPARISON PIPELINE (RF vs LDA vs LASSO)
# ===============================================================

library(caret)
library(randomForest)
library(glmnet)
library(MASS)
library(ggplot2)

set.seed(123)

# ---------------------------
# Data
# ---------------------------
ml_data$group <- as.factor(ml_data$group)

# remove NA if any
ml_data <- na.omit(ml_data)

# predictor matrix + label
x <- ml_data[, -ncol(ml_data)]
y <- ml_data$group

# ===============================================================
# 1. CONTROL SETTINGS
# ===============================================================

ctrl_5fold <- trainControl(
  method = "cv",
  number = 5
)

ctrl_repeated <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5
)

# ===============================================================
# 2. 5-FOLD CROSS VALIDATION MODELS
# ===============================================================

# Random Forest
set.seed(123)
rf_5fold <- train(
  x, y,
  method = "rf",
  trControl = ctrl_5fold,
  ntree = 500
)

# LDA
set.seed(123)
lda_5fold <- train(
  x, y,
  method = "lda",
  trControl = ctrl_5fold
)

# LASSO Logistic Regression
set.seed(123)
lasso_5fold <- train(
  x, y,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl_5fold,
  tuneLength = 10
)

# ===============================================================
# 3. REPEATED 5x5 CROSS VALIDATION MODELS
# ===============================================================

set.seed(123)

rf_rep <- train(
  x, y,
  method = "rf",
  trControl = ctrl_repeated,
  ntree = 500
)

lda_rep <- train(
  x, y,
  method = "lda",
  trControl = ctrl_repeated
)

lasso_rep <- train(
  x, y,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl_repeated,
  tuneLength = 10
)

# ===============================================================
# 4. EXTRACT RESULTS (MEAN ACCURACY)
# ===============================================================

results_5fold <- data.frame(
  Model = c("Random Forest", "LDA", "LASSO"),
  Accuracy = c(
    max(rf_5fold$results$Accuracy),
    max(lda_5fold$results$Accuracy),
    max(lasso_5fold$results$Accuracy)
  ),
  Type = "5-Fold CV"
)

results_repeated <- data.frame(
  Model = c("Random Forest", "LDA", "LASSO"),
  Accuracy = c(
    mean(rf_rep$resample$Accuracy),
    mean(lda_rep$resample$Accuracy),
    mean(lasso_rep$resample$Accuracy)
  ),
  Type = "Repeated 5x5 CV"
)

# combine
all_results <- rbind(results_5fold, results_repeated)

print(all_results)

# ===============================================================
# 5. PLOT COMPARISON (FIGURE 8 STYLE)
# ===============================================================

ggplot(results_5fold, aes(x = Model, y = Accuracy)) +
  geom_point(size = 4) +
  geom_line(group = 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    title = "Model Comparison (5-Fold Cross-Validation)",
    y = "Accuracy"
  )

# ===============================================================
# 6. PLOT REPEATED CV (FIGURE 9 STYLE)
# ===============================================================

ggplot(results_repeated, aes(x = Model, y = Accuracy)) +
  geom_point(size = 4) +
  geom_line(group = 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    title = "Model Comparison (Repeated 5×5 Cross-Validation)",
    y = "Accuracy"
  )

# ===============================================================
# 7. FINAL SUMMARY TABLE
# ===============================================================

aggregate_results <- data.frame(
  Model = c("Random Forest", "LDA", "LASSO"),
  CV_5Fold = results_5fold$Accuracy,
  CV_Repeated = results_repeated$Accuracy
)

print(aggregate_results)
