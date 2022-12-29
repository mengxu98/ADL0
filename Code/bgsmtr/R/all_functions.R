

Wang_algorithm_FUN = function(X, Y, group, r1, r2){

  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))

  eps=2.2204e-16

  d = dim(X)[1]
  p = dim(Y)[1]
  group_set=unique(group)
  group_num=length(group_set)


  ob=rep(NA,group_num)
  ob2=rep(NA,d)

  W_Wang=matrix( 1 , d, p)

  XX=X%*%t(X)
  Xy=X%*%t(Y)

  d1=rep(1,d)
  d2=rep(1,d)
  Wi=matrix(0,d,1)
  obj = c()
  Tol = 10e-4
  iter = 0
  W_change = 1


  while ((W_change > Tol)) {

    iter = iter + 1
    W_Wang_0 = W_Wang
    D1=diag(d1)
    D2=diag(d2)
    W_Wang=solve(XX + r1*D1 + r2*D2)%*%(Xy)

    for (k in 1:group_num){
      idx=which(group==group_set[k])
      idx
      W.k=W_Wang[idx,]
      di=sqrt(sum(W.k*W.k)+eps)
      Wi[idx]=di
      ob[k]=di
    }

    Wi=c(Wi)
    d1=0.5/(Wi)
    Wi2=sqrt(rowSums(W_Wang*W_Wang)+eps)
    d2=0.5/Wi2
    ob2=sum(Wi2)

    W_change = norm((W_Wang-W_Wang_0), type = "F")/max( norm(W_Wang_0, type = "F"), 1)
    obj[iter]=sum(diag((t(t(X)%*%W_Wang-t(Y))%*%(t(X)%*%W_Wang-t(Y)))))+r1*sum(ob)+r2*ob2
    if(iter > 100){break}
  }

  list_to_return = list("W_Wang" = W_Wang, "Wang_obj_func" = obj, 'num_iterations' = iter)

  return(list_to_return)

}



Wang_CV_tuning_values_FUN = function(X, Y, group){

  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))


  d = nrow(X)


  p = nrow(Y)

  gamma_grid=expand.grid( g_1 = c(10e-4, 10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2, 10e3),
                          g_2 = c(10e-4, 10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2, 10e3))

  CV_tuning_results = data.frame(gamma_grid)

  F = 5

  R = length(CV_tuning_results$g_1)

  CV_tuning_results$RMSE_fold_1 = rep(NA, R)
  CV_tuning_results$RMSE_fold_2 = rep(NA, R)
  CV_tuning_results$RMSE_fold_3 = rep(NA, R)
  CV_tuning_results$RMSE_fold_4 = rep(NA, R)
  CV_tuning_results$RMSE_fold_5 = rep(NA, R)
  CV_tuning_results$RMSE_mean = rep(NA, R)

  #data = data.frame(cbind(t(Y), t(X)),tuning_id = rep(NA,nrow(t(Y))))
  data = as.data.frame(cbind(t(Y), t(X)))
  data$tuning_id = rep(NA, nrow(data))

  data$tuning_id = sample(1:F, nrow(data), replace = TRUE)
  cvlist = 1:F


  for (r in 1:R){

    for (f in 1:F){

      training_set = data[data$tuning_id != f,]
      test_set = data[data$tuning_id ==f,]

      X_train = t(training_set[, (p+1): (p+d) ])
      X_test = t(test_set[, (p+1): (p+d)] )

      Y_train = t(training_set[, 1:p ])
      Y_test = t(test_set[, 1:p] )


      W_Wang_hat = Wang_algorithm_FUN( X = X_train, Y = Y_train, group = group,
                                       r1 = CV_tuning_results$g_1[r],
                                       r2 = CV_tuning_results$g_2[r] )$W_Wang

      Y_predict = t(W_Wang_hat)%*%X_test

      RMSE = sqrt(mean((Y_test-Y_predict)^2))


      CV_tuning_results[r, f+2] = RMSE

    }

  }

  CV_tuning_results$RMSE_mean = rowMeans(CV_tuning_results[, 3:(2+F)])

  gamma_values = CV_tuning_results[which(CV_tuning_results$RMSE_mean==min(CV_tuning_results$RMSE_mean)), c(1,2,8)]


  penalty_1 = as.numeric(gamma_values[1])
  penalty_2 = as.numeric(gamma_values[2])

  Wang_algo_results = Wang_algorithm_FUN( X, Y, group, penalty_1, penalty_2 )


  list_to_return = list('CV_tuning_results' = CV_tuning_results,
                        'selected_gamma_values' = gamma_values,
                        'W_Wang_from_tuning_CV' = Wang_algo_results$W_Wang)

  return(list_to_return)

}


bgsmtr.waic = function( X, Y, group, lam_1_fixed, lam_2_fixed, WAIC_opt = TRUE,
                        iter_num = 10000, burn_in = 5001){

  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))


  a_sig_prior=3
  b_sig_prior=1

  mean_prior_sig=b_sig_prior/(a_sig_prior-1)
  var_prior_sig=b_sig_prior^2/((a_sig_prior-1)^2*(a_sig_prior-2))

  Gibbs_setup_return = list( 'iter_num' = iter_num, 'burn_in' = burn_in,
                             'a_sig_prior' = a_sig_prior, 'b_sig_prior' = b_sig_prior,
                             'lam_1_fixed' = lam_1_fixed, 'lam_2_fixed' = lam_2_fixed)


  d = dim(X)[1]
  n = dim(X)[2]

  p = dim(Y)[1]


  group_set=unique(group)
  K = length(group_set)
  m=rep(NA, K)
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }

  idx=list()
  for (k in 1:K){
    idx[[k]]=which(group==group_set[k])
  }

  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p)
  }

  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}
  }


  Xnk=list()
  for (k in 1:K){
    Xnk[[k]]=X[-idx[[k]],]
  }


  Xk_Xk=list()
  for (k in 1:K){
    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }



  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE)
  }


  Xk_Y=list()
  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }



  p_waic = rep(0,n)

  log_p_waic = rep(0,n)

  log_p2_waic = rep(0,n)

  waic_iter = 0



  W_est=array(NA, c(iter_num,d,p))
  tau=array(NA, c(iter_num, K))
  omega=array(NA, c(iter_num, d))
  sig=c(rep(NA, iter_num))

  tau_init_value = 1
  omega_init_value = 1
  sig_init_value = 1


  tau[1,1:K] = rep(tau_init_value ,K)

  omega[1,1:d] = rep( omega_init_value,d)

  sig[1] = sig_init_value

  stm <- proc.time()
  cat('Computing Initial Values \n')

  W_Wang = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = 0.5, r2 = 0.5)$W_Wang
  end_time<-proc.time() - stm
  cat('time: ', end_time[3], 's \n')
  cat('Gibbs Sampler Initialized and Running \n')

  W_est[1,1:d,1:p]=W_Wang

  W_int=W_Wang



  Wk=list()
  for (k in 1:K){
    Wk[[k]]=W_int[idx[[k]],]
  }


  Wnk=list()
  for (k in 1:K){
    Wnk[[k]]=W_int[-idx[[k]],]
  }



  mkp=list()
  Ak=list()
  mu_k=list()
  vec.Wk_t=list()
  Wk_t=list()



  stm <- proc.time()
  for (iter in 1:(iter_num-1)){


    for (k in 1:K){

      Wnk[[k]]=W_int[-idx[[k]],]


      mkp[[k]]<-(1/tau[iter,k]) +(1/omega[iter,idx_mkp[[k]]])

      Ak[[k]]=Xk_Xk[[k]]+.symDiagonal(n=m[k]*p,x=mkp[[k]])

      CH<-Cholesky((1/sig[iter])*Ak[[k]])

      mu_k[[k]]<-solve(Ak[[k]],Xk_Y[[k]]-Xk_Xnk[[k]]%*%as.vector(t(Wnk[[k]])))


      vec.Wk_t[[k]]=rmvn.sparse(n=1, mu=mu_k[[k]], CH=CH,prec=TRUE)


      Wk_t[[k]]=matrix((vec.Wk_t[[k]]), p, m[k])
      Wk[[k]]=t(Wk_t[[k]])

      W_int[idx[[k]],]<-Wk[[k]]

    }


    W_est[(iter+1),1:d,1:p]=W_int


    for (k in 1:K){
      tau[(iter+1),k]=(rinvgauss(1, sqrt((lam_1_fixed*sig[iter])/(t(as.vector(Wk[[k]]))%*%as.vector(Wk[[k]]))), lam_1_fixed ))^-1
    }

    for (i in 1:d){
      omega[(iter+1),i]=(rinvgauss(1, sqrt((lam_2_fixed*sig[iter])/((t(W_est[(iter+1),i,1:p])%*%W_est[(iter+1),i,1:p]))), lam_2_fixed))^(-1)
    }


    Wij2_vk=0
    for (k in 1:K){
      for (i in 1:m[k]){
        Wij2_vk=(t(W_est[(iter+1),idx[[k]][i],1:p])%*%(W_est[(iter+1), idx[[k]][i], 1:p]))*
          ((1/tau[(iter+1),k])+ (1/omega[(iter+1),idx[[k]][i]])) +Wij2_vk
      }
    }

    a_sig=(p*n)/2 + (d*p)/2 + a_sig_prior
    b_sig=(norm((Y-t(W_est[(iter+1),1:d,1:p])%*%X), 'F'))^2/2 + Wij2_vk/2+ b_sig_prior

    sig[(iter+1)]=rinvgamma(1, a_sig, scale = b_sig )


    if (WAIC_opt){
      if(iter >= burn_in){
        waic_iter = waic_iter + 1

        lik.vec<-dmnorm(t(Y),crossprod(X,W_est[(iter+1),1:d,1:p]),varcov=(sig[(iter+1)]*diag(p)))
        p_waic<-p_waic + lik.vec
        log_p_waic<-log_p_waic+log(lik.vec)
        log_p2_waic<-log_p2_waic+(log(lik.vec)^2)
      }
    }

  }
  end_time<-proc.time() - stm
  cat('time: ', end_time[3], 's \n')


  if (WAIC_opt){
    approx_lpd = sum( log ( (1/waic_iter)*p_waic ) )

    approx_P_waic = sum(  (1/(waic_iter -1))*log_p2_waic  -
                            (1/(waic_iter*(waic_iter -1)))*(log_p_waic^2) )

    WAIC = -2*(approx_lpd - approx_P_waic)
  }


  mcmc_tau = mcmc(tau[burn_in:iter_num,1:K])
  mcmc_omega = mcmc(omega[burn_in:iter_num,1:d])
  mcmc_sig = mcmc(sig[burn_in:iter_num])
  mcmc_W_est = list()

  for(q in 1:d) {
    mcmc_W_est[[q]] = mcmc(as.matrix(W_est[burn_in:iter_num,q, 1:p]))
  }


  mcmc_W_est_summaries=list()
  for (q in 1:d){
    mcmc_W_est_summaries[[q]]=summary(mcmc_W_est[[q]])
  }

  mcmc_tau_summary = summary(mcmc_tau)
  mcmc_omega_summary = summary(mcmc_omega)
  mcmc_sig_summary = summary(mcmc_sig)


  W_post_mean=matrix(NA, d, p)
  W_post_sd=matrix(NA, d, p)

  for (q in 1:d){
    W_post_mean[q,]=mcmc_W_est_summaries[[q]][[1]][,1]
    W_post_sd[q,]=mcmc_W_est_summaries[[q]][[1]][,2]
  }


  W_2.5_quantile=matrix(NA, d, p)
  W_97.5_quantile=matrix(NA, d, p)

  for (q in 1:d){
    W_2.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,1]
    W_97.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,5]
  }

  stm <- proc.time()
  cat('Computing Approximate Posterior Mode \n')
  r1.final<-2*mean(sqrt(sig[burn_in:iter_num]))*sqrt(lam_1_fixed)
  r2.final<-2*mean(sqrt(sig[burn_in:iter_num]))*sqrt(lam_2_fixed)
  W_post_mode = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = r1.final, r2 = r2.final)$W_Wang
  end_time<-proc.time() - stm
  cat('time: ', end_time[3], 's \n')

  row.names(W_post_mode) = row.names(W_post_mean) = row.names(W_post_sd) = row.names(W_2.5_quantile) = row.names(W_97.5_quantile) = row.names(X)

  colnames(W_post_mode) = colnames(W_post_mean) = colnames(W_post_sd) = colnames(W_2.5_quantile) = colnames(W_97.5_quantile) = row.names(Y)



  Gibbs_W_summaries_return = list(  'W_post_mean' =  W_post_mean,
                                    'W_post_mode' = W_post_mode,
                                    'W_post_sd' = W_post_sd,
                                    'W_2.5_quantile' =  W_2.5_quantile,
                                    'W_97.5_quantile' =  W_97.5_quantile)


  if (WAIC_opt){ function_returns = list( 'WAIC' = WAIC,  'Gibbs_setup' =   Gibbs_setup_return,
                                          'Gibbs_W_summaries' = Gibbs_W_summaries_return)   } else { function_returns =
                                            list( 'Gibbs_setup' =   Gibbs_setup_return,
                                                  'Gibbs_W_summaries' = Gibbs_W_summaries_return) }



  return(function_returns)

}




bgsmtr.cv.mode = function( X, Y, group, iter_num = 10000, burn_in = 5001){


  WAIC_opt = FALSE

  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))


  a_sig_prior=3
  b_sig_prior=1

  mean_prior_sig=b_sig_prior/(a_sig_prior-1)
  var_prior_sig=b_sig_prior^2/((a_sig_prior-1)^2*(a_sig_prior-2))


  Gibbs_setup_return = list( 'iter_num' = iter_num, 'burn_in' = burn_in,
                             'a_sig_prior' = a_sig_prior, 'b_sig_prior' = b_sig_prior,
                             'lam_1_fixed' = NULL, 'lam_2_fixed' = NULL)



  d = dim(X)[1]
  n = dim(X)[2]

  p = dim(Y)[1]


  group_set=unique(group)
  K = length(group_set)
  m=rep(NA, K)
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }


  idx=list()
  for (k in 1:K){
    idx[[k]]=which(group==group_set[k])
  }

  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p)
  }

  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}
  }


  Xnk=list()
  for (k in 1:K){
    Xnk[[k]]=X[-idx[[k]],]
  }


  Xk_Xk=list()
  for (k in 1:K){

    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }



  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE)
  }


  Xk_Y=list()
  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }


  p_waic = rep(0,n)

  log_p_waic = rep(0,n)

  log_p2_waic = rep(0,n)

  waic_iter = 0


  W_est=array(NA, c(iter_num,d,p))
  tau=array(NA, c(iter_num, K))
  omega=array(NA, c(iter_num, d))
  sig=c(rep(NA, iter_num))

  tau_init_value = 1
  omega_init_value = 1
  sig_init_value = 1


  tau[1,1:K] = rep(tau_init_value ,K)

  omega[1,1:d] = rep( omega_init_value,d)

  sig[1] = sig_init_value

  stm <- proc.time()
  cat('Computing Initial Values and Estimating Tuning Parameters Using Five-Fold Cross-Validation \n')

  CV_results = Wang_CV_tuning_values_FUN(X = X , Y = Y, group = group)
  W_Wang = CV_results$W_Wang_from_tuning_CV
  W_post_mode = W_Wang
  r1.hat = CV_results$selected_gamma_values[1]
  r2.hat = CV_results$selected_gamma_values[2]
  end_time<-proc.time() - stm
  cat('time: ', end_time[3], 's \n')
  cat('Gibbs Sampler Initialized and Running \n')

  W_est[1,1:d,1:p]=W_Wang

  W_int=W_Wang


  lam_1_fixed = (r1.hat/(2*sqrt(sig[1])))^2
  lam_2_fixed = (r2.hat/(2*sqrt(sig[1])))^2



  Wk=list()
  for (k in 1:K){
    Wk[[k]]=W_int[idx[[k]],]
  }


  Wnk=list()
  for (k in 1:K){
    Wnk[[k]]=W_int[-idx[[k]],]
  }

  mkp=list()
  Ak=list()
  mu_k=list()
  vec.Wk_t=list()
  Wk_t=list()


  stm <- proc.time()
  for (iter in 1:(iter_num-1)){

    for (k in 1:K){

      Wnk[[k]]=W_int[-idx[[k]],]

      mkp[[k]]<-(1/tau[iter,k]) +(1/omega[iter,idx_mkp[[k]]])

      Ak[[k]]=Xk_Xk[[k]]+.symDiagonal(n=m[k]*p,x=mkp[[k]])

      CH<-Cholesky((1/sig[iter])*Ak[[k]])

      mu_k[[k]]<-solve(Ak[[k]],Xk_Y[[k]]-Xk_Xnk[[k]]%*%as.vector(t(Wnk[[k]])))


      vec.Wk_t[[k]]=rmvn.sparse(n=1, mu=mu_k[[k]], CH=CH,prec=TRUE)


      Wk_t[[k]]=matrix((vec.Wk_t[[k]]), p, m[k])
      Wk[[k]]=t(Wk_t[[k]])

      W_int[idx[[k]],]<-Wk[[k]]

    }


    W_est[(iter+1),1:d,1:p]=W_int

    for (k in 1:K){
      tau[(iter+1),k]=(rinvgauss(1, sqrt((lam_1_fixed*sig[iter])/(t(as.vector(Wk[[k]]))%*%as.vector(Wk[[k]]))), lam_1_fixed ))^-1
    }

    for (i in 1:d){
      omega[(iter+1),i]=(rinvgauss(1, sqrt((lam_2_fixed*sig[iter])/((t(W_est[(iter+1),i,1:p])%*%W_est[(iter+1),i,1:p]))), lam_2_fixed))^(-1)
    }

    Wij2_vk=0
    for (k in 1:K){
      for (i in 1:m[k]){
        Wij2_vk=(t(W_est[(iter+1),idx[[k]][i],1:p])%*%(W_est[(iter+1), idx[[k]][i], 1:p]))*
          ((1/tau[(iter+1),k])+ (1/omega[(iter+1),idx[[k]][i]])) +Wij2_vk
      }
    }

    a_sig=(p*n)/2 + (d*p)/2 + a_sig_prior
    b_sig=(norm((Y-t(W_est[(iter+1),1:d,1:p])%*%X), 'F'))^2/2 + Wij2_vk/2+ b_sig_prior

    sig[(iter+1)]=rinvgamma(1, a_sig, scale = b_sig )

    lam_1_fixed = (r1.hat/(2*sqrt(sig[(iter+1)])))^2
    lam_2_fixed = (r2.hat/(2*sqrt(sig[(iter+1)])))^2



    if (WAIC_opt){
      if(iter >= burn_in){  # start calculations after burn_in
        waic_iter = waic_iter + 1

        lik.vec<-dmnorm(t(Y),crossprod(X,W_est[(iter+1),1:d,1:p]),varcov=(sig[(iter+1)]*diag(p)))
        p_waic<-p_waic + lik.vec
        log_p_waic<-log_p_waic+log(lik.vec)
        log_p2_waic<-log_p2_waic+(log(lik.vec)^2)
      }
    }


  }
  end_time<-proc.time() - stm
  cat('time: ', end_time[3], 's \n')



  if (WAIC_opt){
    approx_lpd = sum( log ( (1/waic_iter)*p_waic ) )

    approx_P_waic = sum(  (1/(waic_iter -1))*log_p2_waic  -
                            (1/(waic_iter*(waic_iter -1)))*(log_p_waic^2) )

    WAIC = -2*(approx_lpd - approx_P_waic)
  }



  mcmc_tau = mcmc(tau[burn_in:iter_num,1:K])
  mcmc_omega = mcmc(omega[burn_in:iter_num,1:d])
  mcmc_sig = mcmc(sig[burn_in:iter_num])
  mcmc_W_est = list()

  for(q in 1:d) {
    mcmc_W_est[[q]] = mcmc(as.matrix(W_est[burn_in:iter_num,q, 1:p]))
  }


  mcmc_W_est_summaries=list()
  for (q in 1:d){
    mcmc_W_est_summaries[[q]]=summary(mcmc_W_est[[q]])
  }

  mcmc_tau_summary=summary(mcmc_tau)
  mcmc_omega_summary=summary(mcmc_omega)
  mcmc_sig_summary=summary(mcmc_sig)



  W_post_mean=matrix(NA, d, p)
  W_post_sd=matrix(NA, d, p)

  for (q in 1:d){
    W_post_mean[q,]=mcmc_W_est_summaries[[q]][[1]][,1]
    W_post_sd[q,]=mcmc_W_est_summaries[[q]][[1]][,2]
  }



  W_2.5_quantile=matrix(NA, d, p)
  W_97.5_quantile=matrix(NA, d, p)

  for (q in 1:d){
    W_2.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,1]
    W_97.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,5]
  }


  row.names(W_post_mode) = row.names(W_post_mean) = row.names(W_post_sd) = row.names(W_2.5_quantile) = row.names(W_97.5_quantile) = row.names(X)

  colnames(W_post_mode) = colnames(W_post_mean) = colnames(W_post_sd) = colnames(W_2.5_quantile) = colnames(W_97.5_quantile) = row.names(Y)


  Gibbs_W_summaries_return = list(  'W_post_mean' =  W_post_mean,
                                    'W_post_mode' = W_post_mode,
                                    'W_post_sd' = W_post_sd,
                                    'W_2.5_quantile' =  W_2.5_quantile,
                                    'W_97.5_quantile' =  W_97.5_quantile)


  if (WAIC_opt){ function_returns = list( 'WAIC' = WAIC,  'Gibbs_setup' =  Gibbs_setup_return,
                                          'Gibbs_W_summaries' = Gibbs_W_summaries_return)   }

  else { function_returns = list( 'Gibbs_setup' = Gibbs_setup_return, 'Gibbs_W_summaries' = Gibbs_W_summaries_return) }

  return(function_returns)

}

# Spatial Bayesian Group Sparse Multi-Task Regression model with mean field variational bayes method.
# Spatial Bayesian Group Sparse Multi-Task Regression model with mean field variational bayes method.
# Spatial Bayesian Group Sparse Multi-Task Regression model with mean field variational bayes method.
# Spatial Bayesian Group Sparse Multi-Task Regression model with mean field variational bayes method.
sp_bgsmtr_mfvb = function(X, Y, rho = NULL, lambdasq = NULL,
                          alpha = NULL, A = NULL, FDR_opt = TRUE,
                          iter_num = 10000)
{
  # Bayesian Group Sparse Multi-Task Regression Model with MCMC method
  #
  # Args:
  #   X: A d-by-n matrix; d is the number of SNPs and n is the number of subjects.
  #   Y: A c-by-n matrix; c is the number of phenotypes (brain imaging measures) and n is the            number of subjects.
  #   rho: spatial cohesion paramter for spatial correlation between different regions. Value is between 0 and 1.
  #   A: neighbourhood structure for brain region. if not provided, it will bee computed based on       correlation matrix of Y.
  #   lambdasq: tunning paramter. If not given, a default value 1000 is assigned.
  #   FDR_opt: logical operator for computing Bayesian False Dsicovery Rate(FDR).
  #   alpha: Bayesian FDR level. If not given, a default value 0.05 is assigned.
  #   iter_num: Number of iterations with 10000 default value.

  # Output:
  #   MFVB_summaries: mean field variational Bayes approximation results.
  #   FDR_summaries: Bayesian FDR summaries result, which includes significant SNPs, specficity and sensitivity rate for each region.

  num_sub <- dim(X)[2]

  if (num_sub <= 0 || num_sub != dim(Y)[2]) {
    stop("Arguments X and Y have different number of subjects : ",
         dim(X)[2], " and ", dim(Y)[2], ".")
  }

  if (TRUE %in% is.na(X) || TRUE %in% is.na(Y)) {
    stop(" Arguments X and Y must not have missing values.")
  }

  if(is.null(alpha)){
    alpha = 0.05
  }

  if(is.null(rho)){
    rho = 0.8
  }

  if( is.null(lambdasq)){
    lambdasq = 10000
  }



  # Data transformation
  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))

  # Number of Brain Measures
  p = dim(Y)[1]

  # Order the Y such that they are in left and right pairs

  Y = Y[order(rownames(Y)),]
  new_Y = matrix(0, nrow = dim(Y)[1],ncol = dim(Y)[2])

  left_index = seq(1, p, 2)
  right_index = seq(2, p, 2)

  left_names = rownames(Y[1:(p/2),])
  right_names = rownames(Y[(p/2 +1):p,])

  new_Y[left_index,] = Y[1:(p/2),]
  new_Y[right_index,] = Y[(p/2 +1):p,]
  rownames(new_Y) = c(rbind(left_names,right_names))
  colnames(new_Y) = colnames(Y)

  Y = new_Y

  # Neighborhood structure matrix
  if(is.null(A)) {

    A = cor(t(Y[,1:p/2]))
    diag(A) = 0
    A = abs(A)

  }else{
    if ( !is.matrix(A) ){
      stop("Neighborhood structure A has to be symmetric matrix format!")
    }
  }

  # Hyeper_parameters for Sigma
  S_Sig_prior = matrix(c(1,0,0,1),nrow = 2,byrow = TRUE)
  v_Sig_prior = 2


  # number SNPs
  d = dim(X)[1]
  # number subjects
  n = dim(X)[2]


  # unique gene names
  # We are grouping by each SNP, not by gene right now.
  group = 1:d
  group_set = unique(group)
  K = length(group_set)
  m = rep(NA, K) # vector that holds mk number of SNPs in kth gene
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }

  # Next a number of data objects are created
  # these are used repeatedly throughout the mcmc
  # NOTICE THE USE OF SPARSE MATRICES IN SOME PLACES
  # NOTICE KRONECKOR PRODUCT IS %X%

  idx=list()
  for (k in 1:K){
    idx[[k]] = which(group==group_set[k])
  }

  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p/2)
  }

  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}
  }


  Xnk=list()
  for (k in 1:K){
    Xnk[[k]]=X[-idx[[k]],]
  }


  Xk_Xk=list()
  for (k in 1:K){
    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }


  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE)
  }


  Xk_Y=list()
  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }

  Xk_Xk_new = list()
  for (k in 1:K){
    Xk_Xk_new[[k]] = matrix(0, nrow = m[k], ncol = m[k])
    for ( l in 1: n){
      Xk_Xk_new[[k]] = Xk_Xk_new[[k]] + tcrossprod(Xk[[k]][,l])
    }
  }


  Xk_Xnk_new  = list()
  for ( k in 1:K){
    Xk_Xnk_new[[k]] = matrix(0, nrow = m[k], ncol = (d - m[k]))
    for (l in 1:n){
      Xk_Xnk_new[[k]] = Xk_Xnk_new[[k]] + tcrossprod( Xk[[k]][,l], Xnk[[k]][,l] )
    }
  }

  Xk_Y_new  = list()
  for (k in 1:K){
    Xk_Y_new[[k]] =rep(0, m[k]*p)
    for ( l in 1:n){
      Xk_Y_new[[k]] =  Xk_Y_new[[k]] + Xk[[k]][,l] %x% Y[,l]
    }
  }


  D_A = diag(colSums(A))


  # Initial values for VB updating parameters.
  # For each W_k
  # qwk_mu = array(runif(d*p,0,1), c(d,p))
  #USE THE WANG ESTIMATOR TO GET THE
  stm = proc.time()
  cat('Computing Initial Values \n')

  W_Wang = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = 0.5, r2 = 0.5)$W_Wang
  end_time=proc.time() - stm
  cat('time: ', end_time[3], 's \n')
  cat('Mean Field Variational Bayes Initialized and Running \n')

  qwk_mu = W_Wang


  qwk_sigma = list()
  for ( i in 1:d){
    qwk_sigma[[i]] = diag(p)
  }

  # For Sigma
  S_Sig_prior = matrix(c(1,0,0,1),nrow = 2,byrow = TRUE)
  v_Sig_prior = 2

  qsigma_S = array(c(2,0,0,2), c(2,2))
  qsigma_v = 4

  # For omega_i, i=1, ..., d
  qeta_mu = rep(0.5,d)
  qeta_lambda = rep(2,d)
  qomega_mu = 1/qeta_mu + 1/qeta_lambda
  qomega_var = 1/(qeta_lambda*qeta_mu) + 2 / (qeta_lambda)^2

  # For lambda^2
  alpha_lambda_prior = 2
  beta_lambda_prior = 2
  qlambda_alpha = 3
  qlambda_beta = 3

  # For rho
  rho_range = seq(0.01, 0.99, 0.01)
  M = length(rho_range)
  qrho_prob = rep(1/M,M)
  qrho_mu = sum(qrho_prob*rho_range)

  qmu <- rep(0,1000)
  qmu[1] <-  qrho_mu

  log_delta <- 0.5
  iter <- 0
  log_prev <-0
  tol <- 10^(-4)
  mkp=list()
  slogp <- rep(0,1000)
  logd <- rep(0, 1000)
  s_alpha <- rep(0,1000)
  s_beta <- rep(0,1000)
  s_alpha[1] <- alpha_lambda_prior
  s_beta[1] <- beta_lambda_prior


  while( log_delta > tol )

  {
    stm_iter = proc.time()

    iter <- iter +1
    # Updating W_k for k = 1, ... d.
    L.Robert = (D_A - qrho_mu*A) %x% (qsigma_v*solve(qsigma_S))

    for (k in 1:d)
    {
      #  Hk = .symDiagonal(n= (p/2),x = 1/qeta_mu[k]) %x% (qsigma_v*solve(qsigma_S))
      Hk = .symDiagonal(n= (p/2),x = 1/qeta_mu[k]) %x% (qsigma_v*solve(qsigma_S))

      qwk_simga_inv = Matrix(Hk + Xk_Xk_new[[k]] %x% L.Robert, sparse = TRUE)

      qwk_sigma[[k]] = solve(qwk_simga_inv)

      Xk_Xnk_Wnk = Xk_Xnk_new[[k]]%x% L.Robert %*% as.vector(t(qwk_mu[-k,]))

      Xk_DY = diag(1, nrow = m[k]) %x% L.Robert %*% Xk_Y_new[[k]]

      qwk_mu[k,] = solve(qwk_simga_inv, Xk_DY - Xk_Xnk_Wnk, sparse=TRUE )
    }


    eps = 1

    qlambda_beta =  lambdasq / eps

    qlambda_alpha = (lambdasq)^2 / eps



    # Updating for Sigma

    B = D_A - qrho_mu*A

    rs_B = rowSums(B)

    S_term1 = matrix(0,2,2)

    c_idx=seq(1,p/2,1)

    for ( l in 1:n)
    {
      for ( i in c_idx)
      {
        S_term1 = rs_B[i]*Y[(2*i-1):(2*i),l]%*%t(Y[(2*i-1):(2*i),l]) + S_term1
      }
    }

    S_term2 = matrix(0,2,2)

    qwk_ij  = array(0, c(2,2,p/2,d))

    for ( i in 1:K )
    {
      for (j in c_idx)
      {
        qwk_ij[1,1,j,i] = (qwk_mu[i,2*j-1])^2 + qwk_sigma[[i]][(2*j-1), (2*j-1)]

        qwk_ij[2,2,j,i] = (qwk_mu[i,2*j])^2 + qwk_sigma[[i]][2*j,2*j]

        qwk_ij[1,2,j,i] = (qwk_mu[i,2*j])*(qwk_mu[i,(2*j-1)]) + qwk_sigma[[i]][(2*j-1),2*j]

        qwk_ij[2,1,j,i] = (qwk_mu[i,2*j])*(qwk_mu[i,(2*j-1)]) + qwk_sigma[[i]][(2*j-1),2*j]

      }
      S_term2 =  apply(qwk_ij[,,,i], c(1,2), sum)*qeta_mu[i] + S_term2
    }

    qsigma_S = S_term1 + S_term2 + S_Sig_prior

    qsigma_v = 2*n + p*d/2 + v_Sig_prior



    # Update the omega_i variables
    c_idx = seq(1,p/2,1)

    c_star = rep(0,d)

    W_Sig = matrix(0,2,2)

    # update the omega variables
    for (i in 1:d)
    {
      # Updating eta first

      c_star[i] = sum( diag(apply(qwk_ij[,,,i], c(1,2), sum) %*% (qsigma_v*solve(qsigma_S))))

      qeta_mu[i] = sqrt( (qlambda_alpha/qlambda_beta) / c_star[i])

      qeta_lambda[i] = qlambda_alpha / qlambda_beta

      # Updating omega_i,

      qomega_mu[i] = 1 /qeta_mu[i] + 1 / qeta_lambda[i]

      qomega_var[i] =  1 / (qeta_mu[i] * qeta_lambda[i]) + 2 / (qeta_lambda[i])^2
    }

    S_Y_W = list()
    for(l in 1:n)
    {
      S_Y_W[[l]] = tcrossprod(Y[,l] - t(qwk_mu) %*%X[,l])
    }
    S_Y_W = Reduce('+', S_Y_W)

    rho_idx = which(rho_range == rho)
    qrho_prob[rho_idx] = 0.999
    qrho_prob[-rho_idx] =  0.001 / (M-1)
    qrho_mu = sum(qrho_prob*rho_range)
    qmu[iter+1] = qrho_mu

    # Compute the Lower Bound
    LB_term1 = -n/2*log( det(solve(D_A - qrho_mu*A)%x% (1/(qsigma_v-2-1)*qsigma_S))) - 0.5* sum(  diag( S_Y_W %*% ( (D_A - qrho_mu*A) %x% (qsigma_v*solve(qsigma_S)))))

    LB_term2 = 0

    LB_term3 = 0

    for( i in 1:d)
    {

      LB_term2 = -1/2*log( det(qomega_mu[i]* (1/(qsigma_v-2-1)*qsigma_S))) - 1/2*(qeta_mu[i]*c_star[i]) + LB_term2

      LB_term3 = (p+1)/2 *(digamma(qlambda_alpha) - log(qlambda_beta)) + ((p+1)/2 -1)*(log(qomega_mu[i]) - 1/(2*qomega_mu[i]) * qomega_var[i]) - 1/2*qomega_mu[i]*(qlambda_alpha/qlambda_beta) + LB_term3

    }

    LB_term4 = -(v_Sig_prior+3)/2*log(det(1/(qsigma_v-3)*qsigma_S)) -1/2*sum(diag(S_Sig_prior%*%(qsigma_v*solve(qsigma_S))))


    LB_term5 = (alpha_lambda_prior -1)*(digamma(qlambda_alpha) - log(qlambda_beta)) - (beta_lambda_prior + 0.5*sum(qomega_mu))*(qlambda_alpha/qlambda_beta)

    LB_term6 = 1/M

    # Second part of lower bound
    LB_qlog_1 = 0
    LB_qlog_2 = 0

    for(i in 1:d)
    {

      LB_qlog_1 = -1/2* log(det(2*pi*qwk_sigma[[i]])) - 1/2*d + LB_qlog_1

      LB_qlog_2 = 1/2*( log(qeta_lambda[i]) - log(2*pi)) - log(qomega_mu[i]) - 1 / (2*qomega_mu[i]) *qomega_var[i] - qeta_lambda[i]*((qeta_mu[i] - 2) / (2*qeta_mu[i]^2) + qomega_mu[i] /2)
    }

    LB_qlog_3 = qsigma_v/2*log(det(qsigma_S)) - qsigma_v*log(2) - logmvgamma((qsigma_v/2), 2) - (qsigma_v+3)/2*log(det(qsigma_v*qsigma_S)) - 1/2*sum(diag( qsigma_S%*%(qsigma_v*solve(qsigma_S))))

    # LB_qlog_4 = qlambda_alpha*log(qlambda_beta) - log(gamma(qlambda_alpha)) - (qlambda_alpha +1)*(digamma(qlambda_alpha) - log(qlambda_beta))- qlambda_beta*(qlambda_beta/(qlambda_alpha -1) )

    LB_qlog_4 = qlambda_alpha*log(qlambda_beta) - (qlambda_alpha +1)*(digamma(qlambda_alpha) - log(qlambda_beta))- qlambda_beta*(qlambda_beta/(qlambda_alpha -1) )

    LB_qlog_5 = sum(qrho_prob*log(rho_range))

    elbo = (LB_term2+LB_term3+LB_term4+LB_term5+LB_term6)-(LB_qlog_1+LB_qlog_2+LB_qlog_3+LB_qlog_4+LB_qlog_5)
    slogp[iter] = elbo
    log_delta = abs(elbo - log_prev) / log_prev
    logd[iter] = log_delta

    print(c(iter,elbo,log_delta))

    #print(cor(c(qwk_mu),c(W_true)))

    log_prev = elbo

    end_time=proc.time() - stm

    cat('time: ', end_time[3], 's \n')
  }


  mfvb_return = list("Number of Iteration" = iter,
                     "W_post_mean" = qwk_mu,
                     "Sigma_post_mean" = qsigma_S / (qsigma_v - 3),
                     "omega_post_mean" = qomega_mu)


  if(FDR_opt){

    code <- '
    using namespace Rcpp;
    int n = as<int>(n_);
    arma::vec mu = as<arma::vec>(mu_);
    arma::mat sigma = as<arma::mat>(sigma_);
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
    '

    rmvnorm.rcpp <- cxxfunction(signature(n_="integer", mu_="numeric",
                        sigma_="matrix"), code, plugin="RcppArmadillo",
                        verbose=TRUE)

    W_est = array(NA, c(iter_num,p,d))
    n = iter_num
    min_std = rep(0,d)
    prob = matrix(0, nrow = d, ncol = p)

    for (j in 1:d)
    {
      temp = (rmvnorm.rcpp(n, as.vector(qwk_mu[j,]), as.matrix(qwk_sigma[[j]])))
      W_est[1:iter_num,1:p,j] = abs(temp)
      min_std[j] = min(apply(temp,MARGIN = 2, sd))
    }

    for(j in 1:d){

      W_est[1:iter_num,1:p,j] =  apply(W_est[1:iter_num,1:p,j], MARGIN = 2, function(x)  x > min(min_std))
      prob[j,] = apply( W_est[1:iter_num,1:p,j], MARGIN = 2, sum)/n
    }

    prob = replace(prob, prob == 1, 1 - (2*n)^(-1) )
    P.m = prob

    # 2. Get threshold of phi_alpha given alpha
    phi_v = NA
    phi_a = NA
    fdr_rate = rep(NA,p)      # estimated bayesian FDR
    sensitivity_rate = rep(NA,p)  # estimated senstivity
    specificity_rate = rep(NA,p)  # estimated specificity
    significant_snp_idx = matrix(0,d,p) # significant index : 0 for non-significant; 1 for significant
    for(i in 1:p){
      P.i = sort( P.m[,i], decreasing = TRUE)
      P.i2 = length( which( (cumsum(1-P.i)/c(1:d)) <= alpha))
      if(P.i2 > 0){
        phi_v = c(phi_v, max(which((cumsum(1-P.i) / c(1:d)) <= alpha)))
        phi_a = c(phi_a, P.i[max(which((cumsum(1-P.i) / c(1:d)) <= alpha))])
        significant_snp_idx[which(P.m[,i] >= phi_a[(i+1)]),i]=1
        phi.i = which(significant_snp_idx[,i]==1)      #Significant Region for ith ROI.
        fdr_rate[i] = round((cumsum(1-P.i)/c(1:d))[phi_v[(i+1)]],4)
        sensitivity_rate[i] = round(sum(P.m[phi.i,i]) / sum(sum(P.m[phi.i,i]), sum((1-P.m[,i])[-phi.i])), 4)
        specificity_rate[i] = round(sum(P.m[-phi.i,i]) / sum(sum(P.m[-phi.i,i]), sum((1-P.m[,i])[phi.i])), 4)
      }else{
        phi_v = c(phi_v,NA)
        phi_a = c(phi_a,NA)
        fdr_rate[i] = NA
        sensitivity_rate[i] = NA
        specificity_rate[i] = NA
      }
    }
    phi_v = phi_v[-1]
    phi_a = phi_a[-1]
    fdr_vb_result=list(
      fdr_rate = fdr_rate,
      sensitivity_rate = sensitivity_rate,
      specificity_rate = specificity_rate,
      significant_snp_idx = significant_snp_idx)


    function_returns = list("MFVB_summaries" = mfvb_return,
                            "FDR_summaries" = fdr_vb_result,
                            "lower_boud" = elbo)
  }

  else{
    function_returns = list("MFVB_summaries" = mfvb_return, "lower_boud" = elbo)
  }

  return(function_returns)

}
# Spatial Bayesian Group Sparse Multi-Task Regression model with Gibbs Sampling method
sp_bgsmtr_mcmc = function(X, Y, rho = NULL, lambdasq = NULL, alpha = NULL,
                          A = NULL, FDR_opt = TRUE,WAIC_opt = TRUE,
                          iter_num = 10000, burn_in=5001){

  # Bayesian Group Sparse Multi-Task Regression Model with MCMC method.
  #
  # Args:
  #   X: A d-by-n matrix; d is the number of SNPs and n is the number of subjects.
  #   Y: A c-by-n matrix; c is the number of phenotypes (brain imaging measures) and n is       the            number of subjects.
  #   rho: spatial cohesion paramter.
  #   A: neighbourhood structure for brain region. if not provided, it will bee computed        based on       correlation matrix of Y.
  #   lambdasq: tunning paramter.
  #   FDR_opt: logical operator for computing Bayesian False Dsicovery Rate(FDR).
  #   alpha: Bayesian FDR level. Defult is 0.05.
  #   WAIC_opt: logical operator for computing WAIC.
  #   iter_num: Number of iterations for Gibbs sampling.
  #   burn_in: index for burn in.
  #
  # Output:
  #   Gibbs_W_summaries: posterior distribution summaries for W.
  #   FDR_summaries: Bayesian FDR summaries result, which includes significant SNPs,            specficit       y and sensitivity rate for each region.
  #   WAIC: the computed value for waic.
  #   Gibbs_setp_return: set up values for Gibbs sampling.

  num_sub <- dim(X)[2]

  if (num_sub <= 0 || num_sub != dim(Y)[2]) {
    stop("Arguments X and Y have different number of subjects : ",
         dim(X)[2], " and ", dim(Y)[2], ".")
  }

  if (TRUE %in% is.na(X) || TRUE %in% is.na(Y)) {
    stop(" Arguments X and Y must not have missing values.")
  }

  if(is.null(alpha)){
    alpha = 0.05
  }

  if(is.null(rho)){
    rho = 0.8
  }

  if( is.null(lambdasq)){
    lambdasq = 1000
  }

  Gibbs_setup_return = list( 'iter_num' = iter_num, 'burn_in' = burn_in,
                             'rho' = rho, 'lambda_squre' = lambdasq)

  X = t(scale(t(X), center = TRUE, scale = FALSE))
  Y = t(data.frame(scale(t(Y), scale=TRUE, center=TRUE)))

  # Number of Brain Measures
  p = dim(Y)[1]

  # Order the Y such that they are in left and right pairs
  Y = Y[order(rownames(Y)),]
  new_Y = matrix(0, nrow = dim(Y)[1],ncol = dim(Y)[2])

  left_index = seq(1, p, 2)
  right_index = seq(2, p, 2)

  left_names = rownames(Y[1:(p/2),])
  right_names = rownames(Y[(p/2 +1):p,])

  new_Y[left_index,] = Y[1:(p/2),]
  new_Y[right_index,] = Y[(p/2 +1):p,]
  rownames(new_Y) = c(rbind(left_names,right_names))
  colnames(new_Y) = colnames(Y)

  Y = new_Y

  # Neighborhood structure matrix
  if(is.null(A)) {

    A = cor(t(Y[,1:p/2]))
    diag(A) = 0
    A = abs(A)

  }else{
    if (!is.matrix(A)){
      stop("Neighborhood structure A has to be symmetric matrix format!")
    }
  }


  # Hyeper_parameters for Sigma
  S_Sig_prior = matrix(c(1,0,0,1),nrow = 2,byrow = TRUE)
  v_Sig_prior = 2

  # number SNPs
  d = dim(X)[1]
  # number subjects
  n = dim(X)[2]

  group = 1:d
  group_set = unique(group)
  K = length(group_set)
  m = rep(NA, K)  # vector that holds mk number of SNPs in kth gene
  for (k in 1:K){
    m[k]=length(which(group==group_set[k]))
  }

  idx=list()
  for (k in 1:K){
    idx[[k]] = which(group==group_set[k])
  }

  idx_mkp=list()
  for (k in 1:K){
    idx_mkp[[k]]=rep(idx[[k]], each=p/2)
  }

  Xk=list()
  for (k in 1:K){
    Xk[[k]]=as.matrix(X[idx[[k]],])
    if (m[k]==1) {Xk[[k]]=t(Xk[[k]])}
  }


  Xnk=list()
  for (k in 1:K){
    Xnk[[k]] = X[-idx[[k]],]
  }

  Xk_Xk=list()
  for (k in 1:K){
    Xk_Xk[[k]] = Matrix(tcrossprod(Xk[[k]])%x%diag(p),sparse=TRUE)
  }

  Xk_Xnk=list()
  for (k in 1:K){
    Xk_Xnk[[k]] = Matrix(tcrossprod(x=Xk[[k]],y=Xnk[[k]])%x%diag(p), sparse=TRUE)
  }


  Xk_Y=list()

  for (k in 1:K){
    Xk_Y[[k]]=rep(0, m[k]*p)
    for (l in 1:n){
      Xk_Y[[k]]=((Xk[[k]][,l])%x%diag(p))%*%(Y[,l])+Xk_Y[[k]]
    }
  }

  Xk_Xk_new = list()
  for (k in 1:K){
    Xk_Xk_new[[k]] = matrix(0, nrow = m[k], ncol = m[k])
    for ( l in 1: n){
      Xk_Xk_new[[k]] = Xk_Xk_new[[k]] + tcrossprod(Xk[[k]][,l])
    }
  }


  Xk_Xnk_new = list()
  for ( k in 1:K){
    Xk_Xnk_new[[k]] = matrix(0, nrow = m[k], ncol = (d - m[k]))
    for (l in 1:n){
      Xk_Xnk_new[[k]] = Xk_Xnk_new[[k]] + tcrossprod( Xk[[k]][,l], Xnk[[k]][,l] )
    }
  }

  Xk_Y_new = list()
  for (k in 1:K){
    Xk_Y_new[[k]] =rep(0, m[k]*p)
    for ( l in 1:n){
      Xk_Y_new[[k]] =  Xk_Y_new[[k]] + Xk[[k]][,l] %x% Y[,l]
    }
  }

  p_waic_est = array(NA,c(iter_num,n))
  waic_iter = 0

  # Arrays to hold the mcmc samples
  W_est = array(NA, c(iter_num,d,p))
  omega = array(NA, c(iter_num, d))
  Sig  = array(NA,c(iter_num,2,2))

  # Initial values
  omega_init_value = 1
  omega[1,1:d] = rep(omega_init_value,d)

  # Hyeper_parameters for Sigma
  S_Sig_prior = matrix(c(1,0,0,1),nrow = 2,byrow = TRUE)
  v_Sig_prior = 2
  Sig_init_value = rwishart(v_Sig_prior,S_Sig_prior)
  Sig[1 , , ] = Sig_init_value
  inv_Sig = solve( Sig[1 , , ] )

  # USE THE WANG ESTIMATOR TO GET THE
  stm = proc.time()
  cat('Computing Initial Values \n')

  W_Wang = Wang_algorithm_FUN( X = X, Y = Y, group = group, r1 = 0.5, r2 = 0.5)$W_Wang
  end_time=proc.time() - stm
  cat('time: ', end_time[3], 's \n')
  cat('Gibbs Sampler Initialized and Running \n')

  W_est[1,1:d,1:p] = W_Wang

  W_int = W_Wang

  Wk=list()
  for (k in 1:K){
    Wk[[k]]=W_int[idx[[k]],]
  }


  Wnk=list()
  for (k in 1:K){
    Wnk[[k]]=W_int[-idx[[k]],]
  }


  mkp=list()
  Ak=list()
  mu_k=list()
  vec.Wk_t=list()
  Wk_t=list()

  D_A = diag(colSums(A))

  # start mcmc sampler
  stm = proc.time()

  for (iter in 1:(iter_num-1)){

    stm_iter = proc.time()

    # Update each W^k from MVN full conditional

    L.Robert = (D_A - rho*A) %x% inv_Sig

    for (k in 1:K){

      Wnk[[k]] = W_int[-idx[[k]],]

      mkp[[k]] = (1/omega[iter,idx_mkp[[k]]])

      Hk = .symDiagonal(n= (m[k]*p/2),x=mkp[[k]]) %x% inv_Sig

      Wk_prec = Matrix(Hk + Xk_Xk_new[[k]] %x% L.Robert, sparse = TRUE)

      Xk_Xnk_Wnk = Xk_Xnk_new[[k]]%x% L.Robert %*% as.vector(t(Wnk[[k]]))

      Xk_DY = diag(1, nrow = m[k]) %x% L.Robert %*% Xk_Y_new[[k]]

      CH = Cholesky(Wk_prec)

      mu_k[[k]] = solve(Wk_prec, Xk_DY - Xk_Xnk_Wnk, sparse=TRUE )

      # This function generates multivariate normal quickly if precision matrix is sparse
      vec.Wk_t[[k]]=rmvn.sparse(n=1, mu=mu_k[[k]], CH=CH,prec=TRUE)

      # store the sampled values where required
      Wk_t[[k]]=matrix((vec.Wk_t[[k]]), p, m[k])
      Wk[[k]]=t(Wk_t[[k]])
      W_int[idx[[k]],]= Wk[[k]]
    }

    # store the values for this iteration
    W_est[(iter+1),1:d,1:p] = W_int


    # BEGIN UPDATING OMEGA.
    c_idx = seq(1,p/2,1)
    c_star = rep(0,d)

    # update the omega variables
    for (i in 1:d)
    {
      W_Sig = matrix(0,2,2)
      for ( j in c_idx)
      {
        W_Sig = W_est[iter + 1,i,(2*j-1):(2*j)] %*% t( W_est[iter + 1,i,(2*j-1):(2*j)]) %*% inv_Sig + W_Sig
      }
      c_star[i] = sum(diag(W_Sig))
      omega[(iter+1),i] = (rinvgauss(1, sqrt(lambdasq /c_star[i]), lambdasq))^(-1)
    }

    # END UPDATTING FOR OMEGA

    # BEGIN UPDATING FOR SIGMA.

    B = D_A - rho*A
    rs_B = rowSums(B)

    S_term1 = matrix(0,2,2)

    c_idx=seq(1,p/2,1)

    for ( l in 1:n)
    {
      for ( i in c_idx)
      {
        S_term1 = rs_B[i]*Y[(2*i-1):(2*i),l]%*%t(Y[(2*i-1):(2*i),l]) + S_term1
      }
    }

    S_term2 = matrix(0,2,2)

    for (k in 1:K){
      for ( i in idx[[k]]){
        for (j in c_idx){
          S_term2 =  W_est[iter + 1,i,(2*j-1):(2*j)] %*% t(W_est[iter + 1,i,(2*j-1):(2*j)]) * (1/omega[iter + 1, i]) + S_term2
        }
      }
    }

    S_Sig = S_Sig_prior + S_term1 + S_term2
    v_Sig = 2*n + p*d/2 + v_Sig_prior
    Sig[iter+ 1, , ] = rinvwishart(v_Sig, S_Sig)
    inv_Sig = solve(Sig[iter+1 ,,])


    # END UPDATING SIGMA

    S_Y_W = list()
    for(l in 1:n)
    {
      S_Y_W[[l]] = tcrossprod(Y[,l] - t(W_est[(iter+1),1:d,1:p]) %*%X[,l])
    }
    S_Y_W = Reduce('+', S_Y_W)


    # Save lik.vec for computation of WAIC
    if (WAIC_opt){
      if(iter >= burn_in){
        waic_iter = waic_iter + 1
        var_cov=(solve(D_A - rho*A)%x%Sig[iter+1,,])
        s = var_cov
        s.diag = diag(s)
        s[lower.tri(s,diag=T)] = 0
        s = s + t(s) + diag(s.diag)
        lik.vec=dmnorm(t(Y),crossprod(X,W_est[(iter+1),1:d,1:p]),varcov=s)
        p_waic_est[iter+1,]=lik.vec
      }
    }

    end_time=proc.time() - stm
    cat('time: ', end_time[3], 's \n')

  }


  if(length(which((is.na(p_waic_est[-c(1:burn_in),1])*1)==1))==0){
    end_ind=iter_num
  }else{
    end_ind=min(which((is.na(p_waic_est[-c(1:burn_in),1])*1)==1))+burn_in-1
  }

  pwaic = sum(apply(log(p_waic_est[c((burn_in+1):end_ind),]),MARGIN = 2,FUN = var))
  lppd = sum(log(apply(p_waic_est[c((burn_in+1):end_ind),],MARGIN = 2,FUN = mean)))
  waic = 2*pwaic-2*lppd

  mcmc_W_est = list()
  for(q in 1:d) {
    mcmc_W_est[[q]] = mcmc(as.matrix(W_est[(burn_in+1):end_ind,q,1:p]))
  }

  mcmc_W_est_summaries=list()
  for (q in 1:d){
    mcmc_W_est_summaries[[q]]=summary(mcmc_W_est[[q]])
  }


  W_post_mean=matrix(NA, d, p)
  W_post_sd=matrix(NA, d, p)

  for (q in 1:d){
    W_post_mean[q,]=mcmc_W_est_summaries[[q]][[1]][,1]
    W_post_sd[q,]=mcmc_W_est_summaries[[q]][[1]][,2]
  }

  W_2.5_quantile=matrix(NA, d, p)
  W_97.5_quantile=matrix(NA, d, p)

  for (q in 1:d){
    W_2.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,1]
    W_97.5_quantile[q,]=mcmc_W_est_summaries[[q]][[2]][,5]
  }

  row.names(W_post_mean) = row.names(W_post_sd) = row.names(W_2.5_quantile) = row.names(W_97.5_quantile) = row.names(X)

  colnames(W_post_mean) = colnames(W_post_sd) = colnames(W_2.5_quantile) = colnames(W_97.5_quantile) = row.names(Y)

  Gibbs_W_summaries_return = list(
    'W_post_mean' =  W_post_mean,
    'W_post_sd' = W_post_sd,
    'W_2.5_quantile' =  W_2.5_quantile,
    'W_97.5_quantile' =  W_97.5_quantile)


  if (FDR_opt){
    # By default, c.star=min(c(W_post_sd))
    c.star = min(c(W_post_sd))
    # 1. Compute p_{ij}
    P.m = matrix(NA,d,p)
    for(i in 1:d){
      for(j in 1:p){
        P.m[i,j] = sum(abs(W_est[(burn_in+1):end_ind,i,j])>c.star)/(end_ind-burn_in)
        if(P.m[i,j] == 1){
          P.m[i,j]=1- 1/2 /(end_ind-burn_in)
        }
      }
    }

    # 2. Get threshold of phi_alpha given alpha
    phi_v = NA
    phi_a = NA
    fdr_rate = rep(NA,p)      # estimated bayesian FDR
    sensitivity_rate = rep(NA,p)  # estimated senstivity
    specificity_rate = rep(NA,p)  # estimated specificity
    significant_snp_idx = matrix(0,d,p) # significant index : 0 for non-significant; 1 for significant
    for(i in 1:p){
      P.i = sort( P.m[,i], decreasing = TRUE)
      P.i2 = length( which( (cumsum(1-P.i)/c(1:d)) <= alpha))
      if(P.i2 > 0){
        phi_v = c(phi_v, max(which((cumsum(1-P.i) / c(1:d)) <= alpha)))
        phi_a = c(phi_a, P.i[max(which((cumsum(1-P.i) / c(1:d)) <= alpha))])
        significant_snp_idx[which(P.m[,i] >= phi_a[(i+1)]),i] = 1
        phi.i = which(significant_snp_idx[,i] == 1)      #Significant Region for ith ROI.
        fdr_rate[i] = round((cumsum(1-P.i)/c(1:d))[phi_v[(i+1)]],4)
        sensitivity_rate[i] = round(sum(P.m[phi.i,i]) / sum(sum(P.m[phi.i,i]), sum((1-P.m[,i])[-phi.i])), 4)
        specificity_rate[i] = round(sum(P.m[-phi.i,i]) / sum(sum(P.m[-phi.i,i]), sum((1-P.m[,i])[phi.i])), 4)
      }else{
        phi_v = c(phi_v,NA)
        phi_a = c(phi_a,NA)
        fdr_rate[i] = NA
        sensitivity_rate[i] = NA
        specificity_rate[i] = NA
      }
    }
    phi_v = phi_v[-1]
    phi_a = phi_a[-1]
    fdr_mcmc_result=list(
                         fdr_rate = fdr_rate,
                         sensitivity_rate = sensitivity_rate,
                         specificity_rate = specificity_rate,
                         significant_snp_idx = significant_snp_idx)
  }

  if(WAIC_opt){
    if(FDR_opt){
      function_returns = list( 'WAIC' = waic,
                               'Gibbs_setup' = Gibbs_setup_return,
                               'Gibbs_W_summaries' = Gibbs_W_summaries_return,
                               'FDR_summaries' =fdr_mcmc_result)
    }else{
      function_returns = list('WAIC' = waic,
                              'Gibbs_setup' = Gibbs_setup_return,
                              'Gibbs_W_summaries' = Gibbs_W_summaries_return)
    }
  }else{
    if(FDR_opt){
      function_returns = list( 'Gibbs_setup' = Gibbs_setup_return,
                               'Gibbs_W_summaries' = Gibbs_W_summaries_return,
                               'FDR_summaries' =fdr_mcmc_result)
    }else{
      function_returns = list('Gibbs_setup' = Gibbs_setup_return,
                              'Gibbs_W_summaries' = Gibbs_W_summaries_return)
    }
  }

  return(function_returns)
}


#' Bayesian Group Sparse Multi-Task Regression for Imaging Genetics
#'
#' Runs the the Gibbs sampling algorithm to fit a Bayesian group sparse multi-task regression model.
#' Tuning parameters can be chosen using either the MCMC samples and the WAIC (multiple runs) or using an approximation to
#' the posterior mode and five-fold cross-validation (single run).
#'
#'
#' @param X A d-by-n matrix; d is the number of SNPs and n is the number of subjects. Each row of X should correspond to a particular SNP
#' and each column should correspond to a particular subject. Each element of X should give the number of minor alleles for the corresponding
#' SNP and subject. The function will center each row of X to have mean zero prior to running the Gibbs sampling algorithm.
#' @param Y A c-by-n matrix; c is the number of phenotypes (brain imaging measures) and n is the number of subjects. Each row of
#' Y should correspond to a particular phenotype and each column should correspond to a particular subject. Each element of Y should give
#' the measured value for the corresponding phentoype and subject. The function will center and scale each row of Y to have mean zero and unit
#' variance prior to running the Gibbs sampling algorithm.
#' @param group A vector of length d; d is the number of SNPs. Each element of this vector is a string representing a gene or group
#' label associated with each SNP. The SNPs represented by this vector should be ordered according to the rows of X.
#' @param tuning A string, either 'WAIC' or 'CV.mode'. If 'WAIC', the Gibbs sampler is run with fixed values of the tuning
#' parameters specified by the arguments \emph{lam_1_fixed} and  \emph{lam_2_fixed} and the WAIC is computed based on the sampling output. This
#' can then be used to choose optimal values for \emph{lam_1_fixed} and \emph{lam_2_fixed} based on multiple runs with each run using different
#' values of \emph{lam_1_fixed} and \emph{lam_2_fixed}. This option is best suited for either comparing a small set of tuning parameter values or
#' for computation on a high performance computing cluster where different nodes can be used to run the function with different
#' values of \emph{lam_1_fixed} and \emph{lam_2_fixed}. Posterior inference is then based on the run that produces the lowest value for the WAIC.
#' The option 'CV.mode', which is the default, is best suited for computation using just a single processor. In this case the
#' tuning parameters are chosen based on five-fold cross-validation over a grid of possible values with out-of-sample prediction based on an
#' approximate posterior mode. The Gibbs sampler is then run using the chosen values of the tuning parameters. When tuning = 'CV.mode' the values
#' for the arguments \emph{lam_1_fixed} and \emph{lam_2_fixed} are not required.
#' @param lam_1_fixed Only required if tuning = 'WAIC'. A positive number giving the value for the gene-specific tuning parameter. Larger values lead to a larger
#' degree of shrinkage to zero of estimated regression coefficients at the gene level (across all SNPs and phenotypes).
#' @param lam_2_fixed Only required if tuning = 'WAIC'. A positive number giving the value for the SNP-specific tuning parameter. Larger values lead to a larger
#' degree of shrinkage to zero of estimated regression coefficients at the SNP level (across all phenotypes).
#' @param iter_num Positive integer representing the total number of iterations to run the Gibbs sampler. Defaults to 10,000.
#' @param burn_in Nonnegative integer representing the number of MCMC samples to discard as burn-in. Defaults to 5001.
#'
#' @return A list with the elements
#' \item{WAIC}{If tuning = 'WAIC' this is the value of the WAIC computed from the MCMC output. If tuning = 'CV.mode' this component is excluded.}
#' \item{Gibbs_setup}{A list providing values for the input parameters of the function.}
#' \item{Gibbs_W_summaries}{A list with five components, each component being a d-by-c matrix giving some posterior summary of the regression parameter
#' matrix W, where the ij-th element of W represents the association between the i-th SNP and j-th phenotype.
#'
#' -Gibbs_W_summaries$W_post_mean is a d-by-c matrix giving the posterior mean of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_post_mode is a d-by-c matrix giving the posterior mode of W.
#'
#'
#' -Gibbs_W_summaries$W_post_sd is a d-by-c matrix giving the posterior standard deviation for each element of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_2.5_quantile is a d-by-c matrix giving the posterior 2.5 percent quantile for each element of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_97.5_quantile is a d-by-c matrix giving the posterior 97.5 percent quantile for each element of W.'}
#'
#'
#' @author Farouk S. Nathoo, \email{nathoo@uvic.ca}
#' @author Keelin Greenlaw  \email{keelingreenlaw@gmail.com}
#' @author Mary Lesperance  \email{mlespera@uvic.ca}
#'
#' @examples
#' data(bgsmtr_example_data)
#' names(bgsmtr_example_data)
#' \dontshow{
#' ## Toy example with a small subset of the data for routine CRAN testing
#' ## reduce the number of SNPs to 5, subjects to 5, phenotypes to 5
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data[1:5,1:5], Y = bgsmtr_example_data$BrainMeasures[1:5,1:5],
#' group = bgsmtr_example_data$SNP_groups[1:5], tuning = 'WAIC', lam_1_fixed = 2, lam_2_fixed = 2,
#' iter_num = 5, burn_in = 1)
#' }
#'
#' \dontrun{
#' ## test run the sampler for 100 iterations with fixed tunning parameters and compute WAIC
#' ## we recomend at least 5,000 iterations for actual use
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures,
#' group = bgsmtr_example_data$SNP_groups, tuning = 'WAIC', lam_1_fixed = 2, lam_2_fixed = 2,
#' iter_num = 100, burn_in = 50)
#' ## posterior mean for regression parameter relating 100th SNP to 14th phenotype
#' fit$Gibbs_W_summaries$W_post_mean[100,14]
#' ## posterior mode for regression parameter relating 100th SNP to 14th phenotype
#' fit$Gibbs_W_summaries$W_post_mode[100,14]
#' ## posterior standard deviation for regression parameter relating 100th SNP to 14th phenotype
#' fit$Gibbs_W_summaries$W_post_sd[100,14]
#' ## 95% equal-tail credible interval for regression parameter relating 100th SNP to 14th phenotype
#' c(fit$Gibbs_W_summaries$W_2.5_quantile[100,14],fit$Gibbs_W_summaries$W_97.5_quantile[100,14])
#'}
#'
#'\dontrun{
#' ## run the sampler for 10,000 iterations with tuning parameters set using cross-validation
#' ## On a standard computer with a small numer of cores this is the recomended option
#' fit = bgsmtr(X = bgsmtr_example_data$SNP_data, Y = bgsmtr_example_data$BrainMeasures,
#' group = bgsmtr_example_data$SNP_groups, tuning = 'CV.mode',iter_num = 10000, burn_in = 5000)
#'}
#'
#' @references Greenlaw, Keelin, Elena Szefer, Jinko Graham, Mary Lesperance, and Farouk S. Nathoo. "A Bayesian Group Sparse Multi-Task Regression Model for Imaging Genetics." arXiv preprint arXiv:1605.02234 (2016).
#' @references Nathoo, Farouk S., Keelin Greenlaw, and Mary Lesperance. "Regularization Parameter Selection for a Bayesian Multi-Level Group Lasso Regression Model with Application to Imaging Genomics." arXiv preprint arXiv:1603.08163 (2016).
#'
#' @import Matrix mvtnorm
#' @importFrom sparseMVN rmvn.sparse
#' @importFrom statmod rinvgauss
#' @importFrom EDISON rinvgamma
#' @importFrom coda mcmc
#' @importFrom mnormt dmnorm
#'
#'@export
bgsmtr = function(X, Y, group, tuning = 'CV.mode', lam_1_fixed = NULL, lam_2_fixed = NULL, iter_num = 10000, burn_in = 5001)
{
  if (tuning == 'WAIC')
  {
    result = bgsmtr.waic( X=X, Y=Y, group=group, lam_1_fixed=lam_1_fixed, lam_2_fixed=lam_2_fixed, WAIC_opt = TRUE, iter_num = iter_num, burn_in = burn_in)
  }
  else
  {
    result = bgsmtr.cv.mode( X=X, Y=Y, group=group, iter_num = iter_num, burn_in = burn_in)
  }
  return(result)
}

#' @title Spatial Bayesian Group Sparse Multi-Task Regression for Imaging Genetics
#' @description Bayesian Group Sparse Multi-Task Regression that allows for two types of correlation typically seen in structural brain imaging data. First, the spatial correlation in the imaging phenotypes obtained from neighbouring regions of the brain. Second, the correlation between corresponding measures on opposite hemispheres.
#' @param X A d-by-n matrix; d is the number of SNPs and n is the number of subjects. Each row of X should correspond to a particular SNP
#' and each column should correspond to a particular subject. Each element of X should give the number of minor alleles for the corresponding
#' SNP and subject. The function will center each row of X to have mean zero prior to running the Gibbs sampling algorithm.
#' @param Y A c-by-n matrix; c is the number of phenotypes (brain imaging measures) and n is the number of subjects. Each row of
#' Y should correspond to a particular phenotype and each column should correspond to a particular subject. Each element of Y should give
#' the measured value for the corresponding phentoype and subject. The function will center and scale each row of Y to have mean zero and unit
#' variance prior to running the Gibbs sampling algorithm.
#' @param method A string, either 'MCMC' or 'MFVB'. If 'MCMC', the Gibbs sampling method will be used. If 'MFVB', mean field variational bayes method will be used.
#' @param lambdasq A tuning paratmeter. If no value has beeen assigned to it, it takes 1000 by default.
#' @param rho spatial cohesion paramter. If no value has been assigned to it, it takes 0.8 by default.
#' @param alpha Bayesian False Discovery Rate (FDR) level. Default level is 0.05.
#' @param A A c/2 by c/2 neighborhood structure matrix for different brain regions.
#' @param FDR_opt A logical operator for computing Bayesian FDR. By default, it's TRUE.
#' @param WAIC_opt A logical operator for computing WAIC from MCMC method. By default, it's TRUE.
#' @param iter_num Positive integer representing the total number of iterations to run the Gibbs sampler. Defaults to 10,000.
#' @param burn_in Nonnegative integer representing the number of MCMC samples to discard as burn-in. Defaults to 5001.
#' @return A list with the elements
#' \item{WAIC}{WAIC is computed from the MCMC output if "MCMC" is chosen for method.}
#' \item{lower_boud}{Lower bound from MFVB output if "MFVB is choosen for method.}
#' \item{Gibbs_setup}{A list providing values for the input parameters of the function.}
#' \item{Gibbs_W_summaries}{A list with five components, each component being a d-by-c matrix giving some posterior summary of the regression parameter
#' matrix W, where the ij-th element of W represents the association between the i-th SNP and j-th phenotype.
#'
#' -Gibbs_W_summaries$W_post_mean is a d-by-c matrix giving the posterior mean of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_post_mode is a d-by-c matrix giving the posterior mode of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_post_sd is a d-by-c matrix giving the posterior standard deviation for each element of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_2.5_quantile is a d-by-c matrix giving the posterior 2.5 percent quantile for each element of W.
#'
#'
#'
#' -Gibbs_W_summaries$W_97.5_quantile is a d-by-c matrix giving the posterior 97.5 percent quantile for each element of W.'}
#'
#' \item{FDR_summaries}{A list with four components providing the summaries for estimated Bayesian FDR results for both MCMC and MFVB methods. Details for Bayesian FDR computation could be found at Morris et al.(2008).
#'
#' -fdr_rate is the estimated bayesian FDR rate for each region.
#'
#' -sensitivity_rate is the estimated sensitivity rate for each region.
#'
#' -specificity_rate is the estimated specificity rate for each region.
#'
#' -significant_snp_idx is the index of estimated significant/important SNPs for each region.
#' }
#'
#' \item{MFVB_summaries}{A list with four components, each component is the mean field variational bayes approximation summary of model paramters.
#'
#' -Number of Iteration is how many iterations it takes for convergence.
#'
#' -W_post_mean is MFVB approximation of W.
#'
#' -Sigma_post_mean is MFVB approximation of Sigma.
#'
#' -omega_post_mean is MFVB approximation of Omega.
#'
#' }
#'
#'
#' @author Yin Song, \email{yinsong@uvic.ca}
#' @author Shufei Ge \email{shufeig@sfu.ca}
#' @author Farouk S. Nathoo, \email{nathoo@uvic.ca}
#' @author Liangliang Wang \email{lwa68@sfu.ca}
#'
#' @examples
#' data(sp_bgsmtr_example_data)
#' names(sp_bgsmtr_example_data)
#'
#'
#'\dontrun{
#'
#' # Run the example data with Gibbs sampling and compute Bayesian FDR as follow:
#'
#' fit_mcmc = sp_bgsmtr(X = sp_bgsmtr_example_data$SNP_data,
#' Y = sp_bgsmtr_example_data$BrainMeasures, method = "MCMC",
#' A = bgsmtr_example_data$neighborhood_structure, rho = 0.8,
#' FDR_opt = TRUE, WAIC_opt = TRUE,lambdasq = 1000, iter_num = 10000.)
#'
#' # MCMC estimation results for regression parameter W and estimated Bayesian FDR summaries
#'
#' fit_mcmc$Gibbs_W_summaries
#' fit_mcmc$FDR_summaries
#'
#' # The WAIC could be also obtained as:
#'
#' fit_mcmc$WAIC
#'
#' # Run the example data with mean field variational Bayes and compute Bayesian FDR as follow:
#'
#' fit_mfvb = sp_bgsmtr(X = sp_bgsmtr_example_data$SNP_data,
#' Y = sp_bgsmtr_example_data$BrainMeasures, method = "MFVB",
#' A = bgsmtr_example_data$neighborhood_structure, rho = 0.8,FDR_opt = TRUE,
#' lambdasq = 1000, iter_num = 10000.)
#'
#' # MFVB estimated results for regression parameter W and estimated Bayesian FDR summaries
#' fit_mfvb$MFVB_summaries
#' fit_mfvb$FDR_summaries
#'
#' # The corresponding lower bound of MFVB method after convergence is obtained as:
#' fit_mfvb$lower_boud
#'
#'}
#'
#'
#' @import Matrix Rcpp mvtnorm TargetScore miscTools matrixcalc
#' @importFrom LaplacesDemon rinvwishart rwishart
#' @importFrom sparseMVN rmvn.sparse
#' @importFrom inline cxxfunction
#' @importFrom statmod rinvgauss
#' @importFrom EDISON rinvgamma
#' @importFrom coda mcmc
#' @importFrom mnormt dmnorm
#' @importFrom methods  signature
#' @importFrom stats cor sd var
#' @export
sp_bgsmtr = function(X, Y, method = "MCMC", rho = NULL, lambdasq = NULL, alpha = NULL, A = NULL, FDR_opt = TRUE, WAIC_opt = TRUE, iter_num = 10000, burn_in=5001)
{
    if(method == 'MCMC')
    {
       result = sp_bgsmtr_mcmc(X, Y, rho = rho, lambdasq = lambdasq,
                               alpha = alpha,A = A, FDR_opt = TRUE,
                               WAIC_opt = TRUE,iter_num = iter_num, burn_in = burn_in)
    }
    else{
     result = sp_bgsmtr_mfvb(X, Y, rho = rho, lambdasq = lambdasq,
                             alpha = alpha, A = A, FDR_opt = TRUE,
                             iter_num = iter_num)

   }
  return(result)
}
