

library(caret)
library(Metrics)
library(bgsmtr)
library(L0Learn)
source("bgsmtr/R/all_functions.R")
source("Function.R")
# Simulated data with 632 subjects, 486 SNPs from 33 genes, 15 structural neuroimaging measures.
data("bgsmtr_example_data")
X <- bgsmtr_example_data$SNP_data
Y <- bgsmtr_example_data$BrainMeasures
maxSNVSize <- 20

folders <- my_createFolds(ncol(X), 10)
ptm <- proc.time()
rmse_res_l0 <- test_singleROI(X, folders, Y, bgsmtr_example_data$SNP_groups)
running_time <- proc.time() - ptm
running_time
# user  system elapsed 
# 94.45    0.78  103.86

ptm <- proc.time()
rmse_res_Wang <- test_MTR(X, folders, Y, bgsmtr_example_data$SNP_groups)
running_time <- proc.time() - ptm
running_time
# user  system elapsed 
# 3962.19   57.14 4475.31
###########

rmse_pic_list <- comparison_res_rmse(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0)
comparison_bar_rmse(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0, file = "../Results/Supplementary/")

xx <- c(103.86, 0.78)
xx <- rbind(xx, c(4475.30, 57.14))
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
xx <- c(103.86, 4475.30)
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
