library(L0Learn)
set.seed(1) # fix the seed to get a reproducible result
X <- matrix(rnorm(500 * 1000), nrow = 500, ncol = 1000)
B <- c(rep(1, 10), rep(0, 990))
e <- rnorm(500)
y <- X %*% B + e
fit <- L0Learn.fit(X, y, penalty = "L0", maxSuppSize = 20)
print(fit)
y_cat <- predict(fit, newx = X, lambda = 0.0325142, gamma = 0)
rmse <- sum(((y_cat - y)**2)) / 500
# L0L2
fit <- L0Learn.fit(X, y, penalty = "L0L2", nGamma = 5, gammaMin = 0.0001, gammaMax = 10, maxSuppSize = 20)
y_cat <- predict(fit, newx = X, lambda = 0.0011539, gamma = 10)
rmse <- sum(((y_cat - y)**2)) / 500
# Use K-fold cross-validation (CV) to select the optimal values of the tuning parameters
# To set lambda and gamma
cvfit <- L0Learn.cvfit(X, y, nFolds = 5, seed = 1, penalty = "L0L2", nGamma = 5, gammaMin = 0.0001, gammaMax = 0.1, maxSuppSize = 50)
lapply(cvfit$cvMeans, min)
plot(cvfit, gamma = cvfit$fit$gamma[4])
optimalGammaIndex <- 4 # index of the optimal gamma identified previously
optimalLambdaIndex <- which.min(cvfit$cvMeans[[optimalGammaIndex]])
optimalLambda <- cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
optimalLambda
coef(cvfit, lambda = optimalLambda, gamma = cvfit$fit$gamma[4])
y_cat <- predict(cvfit, newx = X, lambda = optimalLambda, gamma = cvfit$fit$gamma[4])
rmse <- sum(((y_cat - y)**2)) / 500
cvfit <- L0Learn.cvfit(X, y, nFolds = 5, seed = 1, penalty = "L0", maxSuppSize = 50)
