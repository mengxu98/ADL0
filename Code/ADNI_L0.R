

library(caret)
library(Metrics)
library(bgsmtr)
library(L0Learn)
source("bgsmtr/R/all_functions.R")
data("bgsmtr_example_data") # Simulated data with 632 subjects, 486 SNPs from 33 genes, 15 structural neuroimaging measures.
X <- bgsmtr_example_data$SNP_data
Y <- bgsmtr_example_data$BrainMeasures
maxSNVSize <- 20

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

my_createFolds <- function(sample_no, folder_no = 10) {
  folders <- list()
  sample_index <- 1:sample_no # sample_no=632
  for (i in 1:folder_no) {
    temp <- sample(x = sample_index, size = sample_no / folder_no)
    folders[[i]] <- temp
  }
  return(folders)
}

predictMTR <- function(W,X)
{
  ###W*X=Y 
  y_prediction=t(W)%*%X
  return(y_prediction)
}

RMSE_MTR <- function(Y,Y_predict)
{
  RMSE_res = c()
  for(i in (1:nrow(Y)))
  {
    temp=rmse(Y[i,],Y_predict[i,])
    RMSE_res = c(RMSE_res,temp)
  }
  return(RMSE_res)
}

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

folders <- my_createFolds(ncol(X), 10)
ptm <- proc.time()
rmse_res_l0 <- test_singleROI(X, folders, Y, bgsmtr_example_data$SNP_groups)
running_time <- proc.time() - ptm
running_time

ptm <- proc.time()
rmse_res_Wang <- test_MTR(X, folders, Y, bgsmtr_example_data$SNP_groups)
running_time <- proc.time() - ptm
running_time
# user   system  elapsed
# 6177.190  292.968 8884.787
###########
library(ggplot2)
comparison_res_rmse <- function(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0) {
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
rmse_pic_list <- comparison_res_rmse(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0)
#
comparison_bar_rmse <- function(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0) {
  x_axis <- 1:length(rmse_res_Wang)
  method_temp <- c("L0L2", "G_SMuRFS")
  mydata <- matrix(0, length(method_temp), length(rmse_res_Wang))
  row.names(mydata) <- method_temp
  colnames(mydata) <- x_axis
  for (j in 1:length(rmse_res_Wang[[1]])) {
    y1 <- c()
    y2 <- c()
    if (j %% 3 == 1) {
      print(paste("../Results/Supplementary/Data_S1_", floor(j / 3) + 1, ".png", sep = ""))
      png(
        file = paste("../Results/Supplementary/Data_S1_", floor(j / 3) + 1, ".png", sep = ""),
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
comparison_bar_rmse(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0)

xx <- c(105.521, 21.899)
xx <- rbind(xx, c(6177.190, 36.826))
row.names(xx) <- c("L0L2", "G_SMuRFS")
n <- nrow(xx)
cols <- rev(gray(0:(n + 1) / (n + 1)))[1:n]
barplot(xx,
  col = cols, beside = TRUE,
  main = "Running time", axisnames = T, ylab = "Seconds(s)",
  ylim = c(0, max(xx) * 1.1)
)
box()
############
xx <- c(105.521, 6177.190)
names(xx) <- c("L0L2", "G_SMuRFS")
n <- length(xx)
cols <- rev(gray(0:(n + 1) / (n + 1)))[1:n]
barplot(xx,
  col = cols, beside = TRUE,
  main = "Running time", axisnames = T, ylab = "Seconds(s)",
  ylim = c(0, max(xx) * 1.1)
)
box()

save(rmse_res_Wang, rmse_res_l0, file = "../Data/bgsmtr_example_data_results.Rdata")
# load("../Data/bgsmtr_example_data_results.Rdata")
