

source("Function.R")
Explanation_r2(X2, Y2, 10)
Explanation_r2(X, Y, 20)
feature_list <- FeatureSummarize(X, Y, 20)

j <- 15
SNV_ROIlabel <- matrix(0, nrow = ncol(X), ncol = nrow(X) + 1)
SNV_ROIlabel[, 1:(ncol(SNV_ROIlabel) - 1)] <- t(X)
colnames(SNV_ROIlabel) <- c(row.names(X), row.names(Y)[j])
SNV_ROIlabel[, ncol(SNV_ROIlabel)] <- Y[j, ]
# SNV_ROIlabel[1:5,480:487]
library(psych)
x_temp <- SNV_ROIlabel[, "rs3783526"]
y_temp <- Y[j, ]
cor_results <- corr.test(x_temp, y_temp, method = "spearman", adjust = "none")
