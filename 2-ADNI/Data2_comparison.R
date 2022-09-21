Data2_path = '2-ADNI\\Dataset2-HDBIG-SCCA-v1.0.0\\'
X2 = t(read.csv(paste0(Data2_path,'SNPdata.csv'),header = T,row.names = 1))
temp=read.csv(paste0(Data2_path,'BrainMeasure.csv'),header = T,row.names = 1)
Y2=t(temp[,1:50])
Group2=read.csv(paste0(Data2_path,'SNPGroup.csv'),header = T)
Group2=Group2[1:nrow(X2),'group_NO']
Group2=as.character(Group2)
maxSNVSize=10
##################
test_singleROI2 <- function(X,folders,Y)
{
  i=1
  X=X2
  Y=Y2
  Group=Group2
  ###
  rmse_res = list()
  folder_no = length(folders)
  allsample_set=1:ncol(X)
  for(i in 1:folder_no)
  {
    test_idx = folders[[i]]
    train_idx = setdiff(allsample_set,test_idx)
    fit_L0L2=fit_singleROI_L0(X[,train_idx],Y[,train_idx],maxSNVSize)
    print(fit_L0L2)
    y_pre_test =predict_singleROI2(X=X[,test_idx],fit=fit_L0L2)
    temp=RMSE_MTR2(Group2[,test_idx],y_pre_test)
    rmse_res[[i]]=temp
  }
  return(rmse_res)
}
predict_singleROI2 <- function(X,fit)
{
  Y_hat =matrix(0,length(fit),ncol(X))#行为ROI数量，列为测试样本数量
  for(i in (1:length(fit)))
  {
    y_cat =predict(fit[[i]]$cvfit, 
                   newx=t(X), #predict函数要求行为样本。所以需转置
                   lambda=fit[[i]]$lambda, gamma=fit[[i]]$gamma)
    Y_hat[i,]=as.vector(y_cat)
  }
  return(Y_hat)
}
library(Metrics)
RMSE_MTR2 <- function(Y,Y_predict)
{
  Y=Group2[,test_idx]
  Y_predict=y_pre_test
  RMSE_res = c()
  for(i in (1:nrow(Y)))
  {
    temp=rmse(Y[i,],Y_predict[i,])
    RMSE_res = c(RMSE_res,temp)
  }
  return(RMSE_res)
}
####################
folders=my_createFolds(ncol(X2), 10)
ptm = proc.time()
rmse_res_l0=test_singleROI(X2,folders,Y2,Group2)
running_time=proc.time() - ptm
# 用户   系统   流逝 
# 21.899  0.034 21.936 
ptm = proc.time()
rmse_res_Wang=test_MTR(X2,folders,Y2,Group2)
running_time=proc.time() - ptm
# user  system elapsed 
# 36.826   4.396  41.228 
###########
comparison_bar_rmse2 <- function(L0L2=rmse_res_Wang,G_SMuRFS=rmse_res_l0)
{
  x_axis=1:length(rmse_res_Wang)
  method_temp=c('L0L2','G_SMuRFS')
  mydata = matrix(0,length(method_temp),length(rmse_res_Wang))
  row.names(mydata)=method_temp
  colnames(mydata)=x_axis
  # par(mfrow=c(1,3))
  for(j in 1:length(rmse_res_Wang[[1]]))#ROI数量
  {
    y1=c()
    y2=c()
    if(j %% 3 == 1) {
      print(paste("D:\\test\\AD\\2-ADNI\\supplementary\\Data_S2_", floor(j / 3) + 1, ".png", sep=""))
      png(file=paste("D:\\test\\AD\\2-ADNI\\supplementary\\Data_S2_", floor(j / 3) + 1, ".png", sep=""), 
          width=1500, height=950, res=105)
      par(mfrow=c(1,3))
    }
    for(i in 1:length(rmse_res_Wang))
    {
      y1=c(y1,rmse_res_l0[[i]][j])
      y2=c(y2,rmse_res_Wang[[i]][j])
    }
    mydata['L0L2',] = y1
    mydata['G_SMuRFS',] = y2
    ########
    n = nrow(mydata)
    cols = rev(gray(0:(n + 1)/(n + 1)))[1:n]
    barplot(mydata, col = cols, beside = TRUE,ylab='RMSE',xlab='Folds',
            main = paste0("",gsub('.adj','',row.names(Y2)[j])),axisnames = T,
            ylim = c(0, max(mydata) * 1.1), cex.axis=1.5, cex.name=1.5, cex.lab=1.5, cex.main=1.8)
    legend("topleft", legend = method_temp, fill = cols, box.col = NA, cex=1.5)
    # legend('top', inset = .01,
    #        legend = method_temp, fill = cols, box.col = NA,
    #        x.intersp=1,y.intersp=0.5,bty="n")
    box()
    if(j %% 3 == 0) {
      dev.off()
    }
    ########
  }
}
comparison_bar_rmse2(L0L2=rmse_res_Wang,G_SMuRFS=rmse_res_l0)

