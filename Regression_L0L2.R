
library(bgsmtr)
library(L0Learn)
source("2-ADNI/bgsmtr/R/all_functions.R")
load("adni_sample_snv_80000.Rdata")
load("all_sample_label.Rdata")

# data_other_gene<- sample_gene_matrix2
data_other_gene <- adni_sample_snv_80000
data_other_gene <- data_other_gene[row.names(all_sample_label), ]
data_other_gene <- data.frame(data_other_gene, check.names = TRUE)
data_other_gene$label <- 0
data_other_gene[row.names(all_sample_label), "label"] <- all_sample_label[, "label"]
data_other_gene[which(data_other_gene$label == "CN"), "label"] <- 0
data_other_gene[which(data_other_gene$label == "MCI"), "label"] <- 1
data_other_gene[which(data_other_gene$label == "AD"), "label"] <- 2
print("label data done")
data_other_gene$label <- as.numeric(data_other_gene$label)
########
SNV_no <- ncol(data_other_gene) - 1
X_SNV <- as.matrix(data_other_gene[, 1:SNV_no])
row.names(X_SNV) <- row.names(data_other_gene)
Y_label <- as.vector(data_other_gene[, SNV_no + 1])
names(Y_label) <- row.names(data_other_gene)
maxSNVSize <- 100
#########
print("L0 learn")
cvfit <- L0Learn.cvfit(X_SNV, Y_label,
  penalty = "L0L2", nGamma = 5, gammaMin = 0.0001,
  gammaMax = 10, maxSuppSize = maxSNVSize
)
print("L0 done")
j <- which.min(lapply(cvfit$cvMeans, min))
optimalGammaIndex <- j # index of the optimal gamma identified previously
optimalLambdaIndex <- which.min(cvfit$cvMeans[[optimalGammaIndex]])
optimalLambda <- cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
print("predicting...")
y_cat <- predict(cvfit,
  newx = X_SNV,
  lambda = optimalLambda, gamma = cvfit$fit$gamma[j]
)
print("predict done")
y_hat <- as.vector(y_cat)

library(Metrics)
res_rmse <- rmse(Y_label, y_hat)
res_rse <- rse(Y_label, y_hat)
r_square <- 1 - res_rse
print("rmse done")
############
temp <- coef(cvfit, lambda = optimalLambda, gamma = cvfit$fit$gamma[j])
temp <- as.vector(temp)
temp <- temp[-1]
temp <- which(temp != 0)
temp <- colnames(X_SNV)[temp]
X_Y <- cbind(X_SNV[, temp], Y_label)
X_Y_frame <- as.data.frame(X_Y)
lmfit <- lm(Y_label ~ ., data = X_Y_frame)
fit_temp <- summary(lmfit)
########
write.csv(fit_temp$coefficients, file = "feature_selected_byL0.csv")

feature_SNV_set <- SNP_top50$V1
feature_SNV_set <- temp
sample_gene_test <- X_SNV[, feature_SNV_set]
sample_gene_test[which(sample_gene_test == 0)] <- " "
sample_gene_test[which(sample_gene_test == 1)] <- "Heterozygote"
sample_gene_test[which(sample_gene_test == 2)] <- "Mutanthomozygote"
mat <- t(as.matrix(sample_gene_test))
col <- c(Heterozygote = "blue", Mutanthomozygote = "red")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = "#CCCCCC", col = NA)
    )
  },
  # big blue
  Heterozygote = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = col["Heterozygote"], col = NA)
    )
  },
  # bug red
  Mutanthomozygote = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = col["Mutanthomozygote"], col = NA)
    )
  }
)
column_title <- "Feature SNPs selected OncoPrint"
column_title <- "Top 50 SNPs OncoPrint"
heatmap_legend_param <- list(
  title = "Alternations", at = c("Heterozygote", "Mutanthomozygote"),
  labels = c("Heterozygote", "Mutanthomozygote")
)
annotation_col <- data.frame(Class = factor(all_sample_label[, "label"])) # 27 AD samples，73 CN samples
rownames(annotation_col) <- row.names(all_sample_label)

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing = TRUE, index.return = TRUE)$ix
  scoreCol <- function(x) {
    score <- 0
    for (i in 1:length(x)) {
      if (x[i]) {
        score <- score + 2^(length(x) - i)
      }
    }
    return(score)
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol)
  sampleOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
  return(M[geneOrder, sampleOrder])
}

mat_encode <- t(X_SNV[, feature_SNV_set])
mat_encode_CN <- memoSort(mat_encode[, which(all_sample_label$label == "CN")])
mat_encode_MCI <- memoSort(mat_encode[, which(all_sample_label$label == "MCI")])
mat_encode_AD <- memoSort(mat_encode[, which(all_sample_label$label == "AD")])
mat_encode_all <- cbind(cbind(as.matrix(mat_encode_CN), as.matrix(mat_encode_MCI)), as.matrix(mat_encode_AD))

annotation_col_order <- data.frame(c(
  all_sample_label[which(all_sample_label$label == "CN"), "label"],
  all_sample_label[which(all_sample_label$label == "MCI"), "label"],
  all_sample_label[which(all_sample_label$label == "AD"), "label"]
))
row.names(annotation_col_order) <- c(
  all_sample_label[which(all_sample_label$label == "CN"), "sampleID"],
  all_sample_label[which(all_sample_label$label == "MCI"), "sampleID"],
  all_sample_label[which(all_sample_label$label == "AD"), "sampleID"]
)
colnames(annotation_col_order) <- c("Class")

mat_encode_all[which(mat_encode_all == 0)] <- " "
mat_encode_all[which(mat_encode_all == 1)] <- "Heterozygote"
mat_encode_all[which(mat_encode_all == 2)] <- "Mutanthomozygote"

oncoPrint(mat_encode_all,
  alter_fun = alter_fun, col = col,
  remove_empty_columns = TRUE, remove_empty_rows = TRUE,
  top_annotation = HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(),
    df = annotation_col_order,
    col = list(Class = c("CN" = "green", "MCI" = "yellow", "AD" = "red"))
  ),
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param,
  alter_fun_is_vectorized = FALSE, column_order = colnames(mat_encode_all)
)


feature_SNV_set <- temp
sample_gene_test <- X_SNV[, feature_SNV_set]
# 0 as “ ”; 1 as "MODERATE"; 2 as "HIGH"
sample_gene_test[which(sample_gene_test == 0)] <- "Wildhomozygote"
sample_gene_test[which(sample_gene_test == 1)] <- "Heterozygote"
sample_gene_test[which(sample_gene_test == 2)] <- "Mutanthomozygote"

mat <- t(as.matrix(sample_gene_test))
col <- c(Wildhomozygote = "green", Heterozygote = "blue", Mutanthomozygote = "red")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = "#CCCCCC", col = NA)
    )
  },
  Wildhomozygote = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = col["Wildhomozygote"], col = NA)
    )
  },
  # big blue
  Heterozygote = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = col["Heterozygote"], col = NA)
    )
  },
  # bug red
  Mutanthomozygote = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
      gp = gpar(fill = col["Mutanthomozygote"], col = NA)
    )
  }
)
###### column_title, heatmap_legend_param
column_title <- "OncoPrint for ADNI, genes in geneset"
heatmap_legend_param <- list(
  title = "Alternations", at = c("Wildhomozygote", "Heterozygote", "Mutanthomozygote"),
  labels = c("Wildhomozygote", "Heterozygote", "Mutanthomozygote")
)
annotation_col <- data.frame(Class = factor(all_sample_label[, "label"])) # 27 AD samples，73 CN samples
rownames(annotation_col) <- row.names(all_sample_label)
library(ComplexHeatmap)
oncoPrint(mat,
  alter_fun = alter_fun, col = col,
  remove_empty_columns = TRUE, remove_empty_rows = TRUE,
  top_annotation = HeatmapAnnotation(
    cbar = anno_oncoprint_barplot(),
    df = annotation_col
  ),
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param,
  alter_fun_is_vectorized = FALSE
)
SNP_top50 <- read.table("figure 4.txt", head = FALSE)
