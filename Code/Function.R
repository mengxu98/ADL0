#' Title
#'
#' @param X 
#' @param Y 
#' @param maxSNVSize 
#'
#' @return
#' @export
#'
#' @examples
fit_singleROI_L0 <- function(X, Y, maxSNVSize) {
  W_singleROI <- list()
  for (i in (1:nrow(Y))) { # the i-th ROI
    cvfit <- L0Learn.cvfit(t(X), t(Y[i, ]), penalty = "L0L2", nGamma = 5, gammaMin = 0.0001, gammaMax = 10, maxSuppSize = maxSNVSize)
    j <- which.min(lapply(cvfit$cvMeans, min))
    # plot(cvfit, gamma=cvfit$fit$gamma[j])
    optimalGammaIndex <- j # index of the optimal gamma identified previously
    optimalLambdaIndex <- which.min(cvfit$cvMeans[[optimalGammaIndex]])
    optimalLambda <- cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
    coef(cvfit, lambda = optimalLambda, gamma = cvfit$fit$gamma[j])
    W_singleROI[[i]] <- list("cvfit" = cvfit, "lambda" = optimalLambda, "gamma" = cvfit$fit$gamma[j])
    print(paste(i, ": optimalLambda", optimalLambda, "; gamma:", cvfit$fit$gamma[j]))
  }
  return(W_singleROI)
}


#' Title
#'
#' @param X 
#' @param fit 
#'
#' @return
#' @export
#'
#' @examples
predict_singleROI <- function(X, fit) {
  Y_hat <- matrix(0, length(fit), ncol(X))
  for (i in (1:length(fit))) {
    y_cat <- predict(fit[[i]]$cvfit,
                     newx = t(X),
                     lambda = fit[[i]]$lambda, gamma = fit[[i]]$gamma
    )
    Y_hat[i, ] <- as.vector(y_cat)
  }
  return(Y_hat)
}

#' Title
#'
#' @param sample_no 
#' @param folder_no 
#'
#' @return
#' @export
#'
#' @examples
my_createFolds <- function(sample_no, folder_no = 10) {
  folders <- list()
  sample_index <- 1:sample_no # sample_no=632
  for (i in 1:folder_no) {
    temp <- sample(x = sample_index, size = sample_no / folder_no)
    folders[[i]] <- temp
  }
  return(folders)
}

#' Title
#'
#' @param W 
#' @param X 
#'
#' @return
#' @export
#'
#' @examples
predictMTR <- function(W,X) {
  # W*X=Y 
  y_prediction=t(W)%*%X
  return(y_prediction)
}

#' Title
#'
#' @param Y 
#' @param Y_predict 
#'
#' @return
#' @export
#'
#' @examples
RMSE_MTR <- function(Y,Y_predict) {
  RMSE_res = c()
  for(i in (1:nrow(Y))) {
    temp=rmse(Y[i,],Y_predict[i,])
    RMSE_res = c(RMSE_res,temp)
  }
  return(RMSE_res)
}

#' Title
#'
#' @param X 
#' @param folders 
#' @param Y 
#' @param Group 
#'
#' @return
#' @export
#'
#' @examples
test_MTR <- function(X, folders, Y, Group) {
  ptm <- proc.time()
  rmse_res <- list()
  folder_no <- length(folders)
  allsample_set <- 1:ncol(X)
  for (i in 1:folder_no) {
    # i=1
    print(c("test_MTR, folder: ", i))
    test_idx <- folders[[i]]
    train_idx <- setdiff(allsample_set, test_idx)
    fit <- Wang_CV_tuning_values_FUN(X[, train_idx], Y[, train_idx],
                                     group = Group
    )
    W_estimate <- fit$W_Wang_from_tuning_CV
    y_prediction <- predictMTR(W_estimate, X[, test_idx])
    temp <- RMSE_MTR(Y[, test_idx], y_prediction)
    rmse_res[[i]] <- temp
  }
  running_time <- proc.time() - ptm
  print(running_time)
  return(rmse_res)
}

#' Title
#'
#' @param X 
#' @param folders 
#' @param Y 
#' @param Group 
#'
#' @return
#' @export
#'
#' @examples
test_singleROI <- function(X, folders, Y, Group) {
  rmse_res <- list()
  global_list_fit <- list()
  folder_no <- length(folders)
  allsample_set <- 1:ncol(X)
  for (i in 1:folder_no) {
    # i=1
    print(c("test_singleROI, folder: ", i))
    test_idx <- folders[[i]]
    train_idx <- setdiff(allsample_set, test_idx)
    fit_L0L2 <- fit_singleROI_L0(X[, train_idx], Y[, train_idx], maxSNVSize)
    global_list_fit[[i]] <- fit_L0L2
    y_pre_test <- predict_singleROI(X = X[, test_idx], fit = fit_L0L2)
    temp <- RMSE_MTR(Y[, test_idx], y_pre_test)
    rmse_res[[i]] <- temp
  }
  list_to_return <- list("rmse_res" = rmse_res, "fit_list" = global_list_fit)
  return(rmse_res)
}

#' Title
#'
#' @param L0L2 
#' @param G_SMuRFS 
#'
#' @return
#' @export
#'
#' @examples
comparison_res_rmse <- function(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0) {
  library(ggplot2)
  rmse_pic_list <- list()
  x_axis <- 1:length(rmse_res_Wang)
  x_axis <- c(x_axis, x_axis)
  method_temp <- c(rep("L0L2", 10), rep("G_SMuRFS", 10))
  for (j in 1:15) {
    y1 <- c()
    y2 <- c()
    for (i in 1:length(rmse_res_Wang)) {
      y1 <- c(y1, rmse_res_l0[[i]][j])
      y2 <- c(y2, rmse_res_Wang[[i]][j])
    }
    y_axis <- c(y1, y2)
    mydata <- data.frame(x_axis, y_axis, method_temp)
    # p <- ggplot(mydata,aes(x=x_axis,y=y_axis,colour=method_temp,group=method_temp,fill=method_temp)) +
    #  geom_line(size =0.8)
    p <- ggplot(mydata, aes(x_axis, y_axis, fill = method_temp)) +
      geom_bar(stat = "identity", position = "dodge") +
      # theme_economist(base_size=14)+
      # scale_fill_economist()+
      theme(axis.ticks.length = unit(0.5, "cm")) +
      guides(fill = guide_legend(title = NULL)) +
      ggtitle(paste0("10-folds of RMSE on ", row.names(Y)[j])) +
      theme(axis.title = element_blank())
    rmse_pic_list[[j]] <- p
  }
  return(rmse_pic_list)
}

#' Title
#'
#' @param L0L2 
#' @param G_SMuRFS 
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
comparison_bar_rmse <- function(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0, file = NULL) {
  if (is.null(file)) {
    file <- ""
  } else {
    if (!dir.exists(file)) {
      dir.create(file, recursive = TRUE)
    }
  }
  x_axis <- 1:length(rmse_res_Wang)
  method_temp <- c("L0L2", "G_SMuRFS")
  mydata <- matrix(0, length(method_temp), length(rmse_res_Wang))
  row.names(mydata) <- method_temp
  colnames(mydata) <- x_axis
  for (j in 1:length(rmse_res_Wang[[1]])) {
    y1 <- c()
    y2 <- c()
    if (j %% 3 == 1) {
      png(
        file = paste(file, "Data_S1_", floor(j / 3) + 1, ".png", sep = ""),
        width = 1500, height = 950, res = 105
      )
      par(mfrow = c(1, 3))
    }
    
    for (i in 1:length(rmse_res_Wang)) {
      y1 <- c(y1, rmse_res_l0[[i]][j])
      y2 <- c(y2, rmse_res_Wang[[i]][j])
    }
    mydata["L0L2", ] <- y1
    mydata["G_SMuRFS", ] <- y2
    
    n <- nrow(mydata)
    cols <- rev(gray(0:(n + 1) / (n + 1)))[1:n]
    
    barplot(mydata,
            col = cols, beside = TRUE, ylab = "RMSE", xlab = "Folds",
            main = paste0("", gsub(".adj", "", row.names(Y)[j])), axisnames = T,
            ylim = c(0, max(mydata) * 1.1), cex.axis = 1.5, cex.name = 1.5, cex.lab = 1.5, cex.main = 1.8
    )
    legend("topleft", legend = method_temp, fill = cols, box.col = NA, cex = 1.5)
    # legend('top', inset = .01,
    #        legend = method_temp, fill = cols, box.col = NA,
    #        x.intersp=1,y.intersp=0.5,bty="n")
    box()
    if (j %% 3 == 0) {
      dev.off()
    }
  }
}


#' Title
#'
#' @param X 
#' @param folders 
#' @param Y 
#'
#' @return
#' @export
#'
#' @examples
test_singleROI2 <- function(X, folders, Y) {
  i <- 1
  X <- X2
  Y <- Y2
  Group <- Group2
  ###
  rmse_res <- list()
  folder_no <- length(folders)
  allsample_set <- 1:ncol(X)
  for (i in 1:folder_no) {
    test_idx <- folders[[i]]
    train_idx <- setdiff(allsample_set, test_idx)
    fit_L0L2 <- fit_singleROI_L0(X[, train_idx], Y[, train_idx], maxSNVSize)
    print(fit_L0L2)
    y_pre_test <- predict_singleROI2(X = X[, test_idx], fit = fit_L0L2)
    temp <- RMSE_MTR2(Group2[, test_idx], y_pre_test)
    rmse_res[[i]] <- temp
  }
  return(rmse_res)
}

#' Title
#'
#' @param X 
#' @param fit 
#'
#' @return
#' @export
#'
#' @examples
predict_singleROI2 <- function(X, fit) {
  Y_hat <- matrix(0, length(fit), ncol(X))
  for (i in (1:length(fit)))
  {
    y_cat <- predict(fit[[i]]$cvfit,
                     newx = t(X),
                     lambda = fit[[i]]$lambda, gamma = fit[[i]]$gamma
    )
    Y_hat[i, ] <- as.vector(y_cat)
  }
  return(Y_hat)
}

#' Title
#'
#' @param Y 
#' @param Y_predict 
#'
#' @return
#' @export
#'
#' @examples
RMSE_MTR2 <- function(Y, Y_predict) {
  library(Metrics)
  Y <- Group2[, test_idx]
  Y_predict <- y_pre_test
  RMSE_res <- c()
  for (i in (1:nrow(Y)))
  {
    temp <- rmse(Y[i, ], Y_predict[i, ])
    RMSE_res <- c(RMSE_res, temp)
  }
  return(RMSE_res)
}

#' Title
#'
#' @param L0L2 
#' @param G_SMuRFS 
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
comparison_bar_rmse2 <- function(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0, file = NULL) {
  if (is.null(file)) {
    file <- ""
  } else {
    if (!dir.exists(file)) {
      dir.create(file, recursive = TRUE)
    }
  }
  x_axis <- 1:length(rmse_res_Wang)
  method_temp <- c("L0L2", "G_SMuRFS")
  mydata <- matrix(0, length(method_temp), length(rmse_res_Wang))
  row.names(mydata) <- method_temp
  colnames(mydata) <- x_axis
  # par(mfrow=c(1,3))
  for (j in 1:length(rmse_res_Wang[[1]])) {
    y1 <- c()
    y2 <- c()
    if (j %% 3 == 1) {
      png(
        file = paste(file, "Data_S2_", floor(j / 3) + 1, ".png", sep = ""),
        width = 1500, height = 950, res = 105
      )
      par(mfrow = c(1, 3))
    }
    for (i in 1:length(rmse_res_Wang)) {
      y1 <- c(y1, rmse_res_l0[[i]][j])
      y2 <- c(y2, rmse_res_Wang[[i]][j])
    }
    mydata["L0L2", ] <- y1
    mydata["G_SMuRFS", ] <- y2
    ########
    n <- nrow(mydata)
    cols <- rev(gray(0:(n + 1) / (n + 1)))[1:n]
    barplot(mydata,
            col = cols, beside = TRUE, ylab = "RMSE", xlab = "Folds",
            main = paste0("", gsub(".adj", "", row.names(Y2)[j])), axisnames = T,
            ylim = c(0, max(mydata) * 1.1), cex.axis = 1.5, cex.name = 1.5, cex.lab = 1.5, cex.main = 1.8
    )
    legend("topleft", legend = method_temp, fill = cols, box.col = NA, cex = 1.5)
    # legend('top', inset = .01,
    #        legend = method_temp, fill = cols, box.col = NA,
    #        x.intersp=1,y.intersp=0.5,bty="n")
    box()
    if (j %% 3 == 0) {
      dev.off()
    }
  }
}

#' Title
#'
#' @param X2 
#' @param Y2 
#' @param maxSNVSize 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X2 
#' @param Y2 
#' @param maxSNVSize 
#'
#' @return
#' @export
#'
#' @examples
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