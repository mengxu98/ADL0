
rmse_res_l0_matrix = matrix(unlist(rmse_res_l0), ncol = length(rmse_res_l0[[1]]), byrow = TRUE)
rmse_res_Wang_matrix = matrix(unlist(rmse_res_Wang), ncol = length(rmse_res_Wang[[1]]), byrow = TRUE)

# 平均RMSE进行paired t-test
rmse_res_l0_mean = colMeans(rmse_res_l0_matrix)
rmse_res_Wang_mean = colMeans(rmse_res_Wang_matrix)
rmse_paired_ttest <- t.test(x=rmse_res_l0_mean, y=rmse_res_Wang_mean, paired=TRUE,alternative="two.sided")
print(rmse_paired_ttest$statistic)
print(rmse_paired_ttest$p.value)

# Data S1 t=1.692417, p=0.05634634,
# Data S2 t=61.11133, p=3.217035e-48

# two side
# S1 t=-1.692417, p=0.1126927 -240.6586
# S2 t=-62.52564, p=6.434069e-48

#每个folder进行paired t-test
for(i in 1:10){
  print(paste("fold", i, ":", t.test(rmse_res_Wang_matrix[, i], rmse_res_l0_matrix[, i], paired=TRUE,alternative="greater")$p.value))
}

plot(rmse_res_Wang_mean, type='l', xlab = "ROIs", ylab = "Average RMSE", col = "red")
lines(rmse_res_l0_mean, type = "l", col = "blue")

#箱线图
library(ggpubr)
ggpaired(data.frame(rmse_res_Wang_mean,rmse_res_l0_mean), cond1="rmse_res_Wang_mean", cond2="rmse_res_l0_mean")

#参考文献paired t-test
paper_train_KG_mean = c(0.461,0.464,0.4588,0.4638,0.4644)
paper_train_PMA_mean = c(0.4318,0.4306,0.4302,0.4314,0.4308)
t.test(paper_train_KG_mean, paper_train_PMA_mean, paired=TRUE,alternative="greater")

#箱线图
rmse_res_group = data.frame(c(rmse_res_l0_mean,rmse_res_Wang_mean),c(matrix(1, nrow=1, ncol=15), matrix(2, nrow=1, ncol=15)))
names(rmse_res_group) <- c("mean", "group")
boxplot(mean ~ group, data=rmse_res_group)

write.csv(data.frame(rmse_res_l0_mean,rmse_res_Wang_mean), file="/Users/linyangkai/Desktop/d.csv")

#Shapiro-Wilk正态性检验差值是否符合正态分布
shapiro_test = shapiro.test(rmse_res_l0_mean-rmse_res_Wang_mean)
print(shapiro_test)
wilcox_test = wilcox.test(rmse_res_l0_mean, rmse_res_Wang_mean, data=paired_data,  paired=TRUE, alternative = "two.sided")
print(wilcox_test)
print(median(rmse_res_l0_mean-rmse_res_Wang_mean))

