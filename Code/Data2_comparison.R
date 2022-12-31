

source("Function.R")
Data2_path <- "../Data/Dataset2-HDBIG-SCCA-v1.0.0/"
X2 <- t(read.csv(paste0(Data2_path, "SNPdata.csv"), header = T, row.names = 1))
temp <- read.csv(paste0(Data2_path, "BrainMeasure.csv"), header = T, row.names = 1)
Y2 <- t(temp[, 1:50])
Group2 <- read.csv(paste0(Data2_path, "SNPGroup.csv"), header = T)
Group2 <- Group2[1:nrow(X2), "group_NO"]
Group2 <- as.character(Group2)
maxSNVSize <- 10

folders <- my_createFolds(ncol(X2), 10)
ptm <- proc.time()
rmse_res_l0 <- test_singleROI(X2, folders, Y2, Group2)
running_time <- proc.time() - ptm

ptm <- proc.time()
rmse_res_Wang <- test_MTR(X2, folders, Y2, Group2)
running_time <- proc.time() - ptm

comparison_bar_rmse2(L0L2 = rmse_res_Wang, G_SMuRFS = rmse_res_l0, file = "../Results/Supplementary/")
