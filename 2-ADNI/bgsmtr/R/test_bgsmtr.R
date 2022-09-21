library(sparseMVN)
library(statmod)
library(EDISON)
library(coda)
library(mnormt)
#load('/Users/xionglee/Downloads/研究生论文指导/2-ADNI/bgsmtr/data/bgsmtr_example_data.RData')
library(bgsmtr)
data(bgsmtr_example_data) #这样使用的前提是先安装对应的package
names(bgsmtr_example_data)
# Not run:
## test run the sampler for 100 iterations with fixed tunning parameters and compute WAIC
## we recomend at least 5,000 iterations for actual use
fit = bgsmtr(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures,
             group = bgsmtr_example_data$SNP_groups, tuning = 'WAIC', lam_1_fixed = 2, lam_2_fixed = 2, iter_num = 100, burn_in = 50)
## posterior mean for regression parameter relating 100th SNP to 14th phenotype 
fit$Gibbs_W_summaries$W_post_mean[100,14]
## posterior mode for regression parameter relating 100th SNP to 14th phenotype 
fit$Gibbs_W_summaries$W_post_mode[100,14]
## posterior standard deviation for regression parameter relating 100th SNP to 14th phenotype 
fit$Gibbs_W_summaries$W_post_sd[100,14]
## 95% equal-tail credible interval for regression parameter relating 100th SNP to 14th phenotype
c(fit$Gibbs_W_summaries$W_2.5_quantile[100,14],fit$Gibbs_W_summaries$W_97.5_quantile[100,14])
###############
fold_predict <- predict(fit,
                        type='response',
                        newdata=bgsmtr_example_data$SNP_data[1:10])








## End(Not run)