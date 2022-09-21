#Simulated data with 632 subjects, 486 SNPs from 33 genes, 15 structural neuroimaging measures.
library(bgsmtr)
library(L0Learn)
data("bgsmtr_example_data")
X = bgsmtr_example_data$SNP_data
Y = bgsmtr_example_data$BrainMeasures
maxSNVSize=20
#########
fit_singleROI_L0<- function(X,Y,maxSNVSize)
{
  W_singleROI = list()
  for(i in (1:nrow(Y)))#the i-th ROI
  {
    cvfit = L0Learn.cvfit(t(X), t(Y[i,]), penalty="L0L2", nGamma = 5, gammaMin = 0.0001, gammaMax = 10, maxSuppSize=maxSNVSize)
    j=which.min(lapply(cvfit$cvMeans, min))
    optimalGammaIndex = j # index of the optimal gamma identified previously
    optimalLambdaIndex = which.min(cvfit$cvMeans[[optimalGammaIndex]])
    optimalLambda = cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
    coef(cvfit, lambda=optimalLambda, gamma=cvfit$fit$gamma[j])
    W_singleROI[[i]]=list('cvfit'=cvfit,'lambda'=optimalLambda, 'gamma'=cvfit$fit$gamma[j])
  }
  return(W_singleROI)
}
predict_singleROI <- function(X,fit)
{
  Y_hat =matrix(0,length(fit),ncol(X))
  for(i in (1:length(fit)))
  {
    y_cat =predict(fit[[i]]$cvfit, 
                   newx=t(X), 
                   lambda=fit[[i]]$lambda, gamma=fit[[i]]$gamma)
    Y_hat[i,]=as.vector(y_cat)
  }
  return(Y_hat)
}
my_createFolds <- function(sample_no,folder_no=10)
{
  folders= list()
  #sample_no=632
  sample_index=1:sample_no
  for(i in 1:folder_no)
  {
    temp=sample(x=sample_index,size=sample_no/folder_no)
    folders[[i]]=temp
  }
  return(folders)
}

test_singleROI <- function(X,folders,Y,Group)
{
  rmse_res = list()
  global_list_fit = list()
  folder_no = length(folders)
  allsample_set=1:ncol(X)
  for(i in 1:folder_no)
  {
    test_idx = folders[[i]]
    train_idx = setdiff(allsample_set,test_idx)
    fit_L0L2=fit_singleROI_L0(X[,train_idx],Y[,train_idx],maxSNVSize)
    global_list_fit[[i]] = fit_L0L2
    y_pre_test =predict_singleROI(X=X[,test_idx],fit=fit_L0L2)
    temp=RMSE_MTR(Y[,test_idx],y_pre_test)
    rmse_res[[i]]=temp
  }
  list_to_return = list('rmse_res'=rmse_res,'fit_list'=global_list_fit)
  return(rmse_res)
}
########
folders=my_createFolds(ncol(X),10)
ptm = proc.time()
rmse_res_l0=test_singleROI(X,folders,Y,bgsmtr_example_data$SNP_groups)
running_time=proc.time() - ptm

