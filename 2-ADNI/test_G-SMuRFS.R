predictMTR <- function(W,X)
{
 ###W*X=Y 
  y_prediction=t(W)%*%X
  return(y_prediction)
}
library(Metrics)
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
fit = Wang_CV_tuning_values_FUN(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures,
             group = bgsmtr_example_data$SNP_groups)
W_estimate = fit$W_Wang_from_tuning_CV
y_prediction=predictMTR(W_estimate,bgsmtr_example_data$SNP_data[,1:20])#前20个样本的所有SNP
RMSE_MTR(bgsmtr_example_data$BrainMeasures[,1:20],y_prediction)#测试前20个样本在15个ROI上的精度
#[1] 2.913117e+02 1.602795e+04 2.241469e+04 6.950813e+02 4.367939e+02 7.506675e+03 6.679157e-01 2.010401e-01
#[9] 2.014468e-01 2.191167e-01 2.107629e-01 4.393526e-01 1.706704e-01 1.551464e-01 2.023908e-01
######################
test_MTR <- function(X,folders,Y,Group)
{
  ptm = proc.time()
  rmse_res = list()
  folder_no = length(folders)
  allsample_set=1:ncol(X)
  for(i in 1:folder_no)
  {
    #i=1
    test_idx = folders[[i]]
    train_idx = setdiff(allsample_set,test_idx)
    fit = Wang_CV_tuning_values_FUN(X[,train_idx], Y[,train_idx],
                                    group = Group)
    W_estimate = fit$W_Wang_from_tuning_CV
    y_prediction=predictMTR(W_estimate,X[,test_idx])
    temp=RMSE_MTR(Y[,test_idx],y_prediction)
    rmse_res[[i]]=temp
  }
  running_time=proc.time() - ptm
  print(running_time)
  return(rmse_res)
}

