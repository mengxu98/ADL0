

FeatureSummarize <- function(X2, Y2, maxSNVSize) {
  fit_L0L2 <- fit_singleROI_L0(X2, Y2, maxSNVSize)
  ROI_feature_list <- list()
  for (j in 1:nrow(Y2))
  {
    temp <- coef(fit_L0L2[[j]]$cvfit, lambda = fit_L0L2[[j]]$lambda, gamma = fit_L0L2[[j]]$gamma)
    temp <- as.vector(temp)
    temp <- temp[-1]
    temp <- which(temp != 0)
    temp <- row.names(X2)[temp]
    ROI_name <- gsub(".adj", "", row.names(Y2)[j])
    ROI_feature_list[[ROI_name]] <- temp
  }
  return(ROI_feature_list)
}

Explanation_r2 <- function(X2, Y2, maxSNVSize) {
  fit_L0L2 <- fit_singleROI_L0(X2, Y2, maxSNVSize)
  ROI_feature_list <- list()
  for (j in 1:nrow(Y2)) {
    temp <- coef(fit_L0L2[[j]]$cvfit, lambda = fit_L0L2[[j]]$lambda, gamma = fit_L0L2[[j]]$gamma)
    temp <- as.vector(temp)
    temp <- temp[-1]
    temp <- which(temp != 0)
    temp <- row.names(X2)[temp]
    ROI_name <- gsub(".adj", "", row.names(Y2)[j])
    ROI_feature_list[[ROI_name]] <- temp
  }
  
  for (j in 1:nrow(Y2)) {
    if (length(ROI_feature_list[[j]] > 0)) {
      feature_SNV <- ROI_feature_list[[j]]
      feature_ROI <- matrix(0, nrow = ncol(X2), ncol = length(feature_SNV) + 1)
      colnames(feature_ROI) <- c(feature_SNV, row.names(Y2)[j]) #<U+7B2C><U+7B2C>j<U+4E2A>ROI<U+4F5C><U+6D4B><U+8BD5>
      feature_ROI[, feature_SNV] <- t(X2[feature_SNV, ])
      feature_ROI[, ncol(feature_ROI)] <- Y2[j, ]
      # plot(feature_ROI[,feature_SNV],feature_ROI[,ncol(feature_ROI)])
      lmfit <- lm(feature_ROI[, ncol(feature_ROI)] ~ feature_ROI[, feature_SNV])
      # lmfit = lm(Left_AmygVol.adj~rs596577,data=feature_ROI)
      # abline(lmfit,col= "red")
      temp <- summary(lmfit) # Multiple R-squared:  0.004505,	Adjusted R-squared:  0.002925
      print(temp$r.squared)
      if (temp$r.squared > 0.5) {
        print(paste0(j, ">0.5"))
      }
    }
  }
}

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
