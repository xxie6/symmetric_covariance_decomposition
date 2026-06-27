refit_lambda <- function(S, sym_ebcovmf_obj, maxiter = 100, tol = 10^(-6), remove_null = TRUE){
  K <- length(sym_ebcovmf_obj$lambda)
  if (K <= 1){
    print('Cannot refit lambda')
  } else {
    sym_ebcovmf_obj.old <- sym_ebcovmf_obj
    iter <- 1
    obj_diff <- Inf
    curr_elbo <- -Inf
    while((iter <= maxiter) && (obj_diff > tol)){
      # print(iter)
      # Update lambdas
      lambda.old <- sym_ebcovmf_obj$lambda
      for (k in 1:K){
        Rk <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
        sym_ebcovmf_obj$lambda[k] <- max(t(sym_ebcovmf_obj$L_pm[,k, drop = FALSE]) %*% Rk %*% sym_ebcovmf_obj$L_pm[,k, drop = FALSE], 0)
        # print(sym_ebcovmf_obj$lambda[k])
      }
      # print(sym_ebcovmf_obj$lambda)
      
      # Update resid_s2
      resid_s2.old <- sym_ebcovmf_obj$resid_s2
      sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                    L = sym_ebcovmf_obj$L_pm,
                                                    lambda = sym_ebcovmf_obj$lambda,
                                                    n = sym_ebcovmf_obj$n,
                                                    K = K)
      
      # Update elbo
      curr_elbo.old <- curr_elbo
      # print(curr_elbo.old)
      curr_elbo <- compute_elbo(S = S,
                                L = sym_ebcovmf_obj$L_pm,
                                lambda = sym_ebcovmf_obj$lambda,
                                resid_s2 = sym_ebcovmf_obj$resid_s2,
                                n = sym_ebcovmf_obj$n,
                                K = K,
                                KL = sym_ebcovmf_obj$KL)
      # print(curr_elbo)
      
      # Check convergence
      if (iter > 1){
        obj_diff <- curr_elbo - curr_elbo.old
      }
      if (obj_diff < 0){
        sym_ebcovmf_obj$lambda <- lambda.old
        sym_ebcovmf_obj$resid_s2 <- resid_s2.old
        curr_elbo <- curr_elbo.old
        print(paste('elbo decreased by', abs(obj_diff)))
        break
      }
      iter <- iter + 1
    }
    # nullcheck
    if (any(sym_ebcovmf_obj$lambda == 0) & remove_null == TRUE){
      idx <- which(sym_ebcovmf_obj$lambda != 0)
      # print(idx)
      K <- length(idx) # I just added this
      sym_ebcovmf_obj$lambda <- sym_ebcovmf_obj$lambda[idx]
      sym_ebcovmf_obj$L_pm <- sym_ebcovmf_obj$L_pm[,idx]
      sym_ebcovmf_obj$KL <- sym_ebcovmf_obj$KL[idx]
      sym_ebcovmf_obj$fitted_gs <- sym_ebcovmf_obj$fitted_gs[idx]
      sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                    L = sym_ebcovmf_obj$L_pm,
                                                    lambda = sym_ebcovmf_obj$lambda,
                                                    n = sym_ebcovmf_obj$n,
                                                    K = K)
      curr_elbo <- compute_elbo(S = S,
                                L = sym_ebcovmf_obj$L_pm,
                                lambda = sym_ebcovmf_obj$lambda,
                                resid_s2 = sym_ebcovmf_obj$resid_s2,
                                n = sym_ebcovmf_obj$n,
                                K = K,
                                KL = sym_ebcovmf_obj$KL)
    }
    # check objective function
    # print(curr_elbo)
    # print(sym_ebcovmf_obj.old$elbo)
    if ((curr_elbo - sym_ebcovmf_obj.old$elbo) < 0){
      sym_ebcovmf_obj <- sym_ebcovmf_obj.old
    } else {
      sym_ebcovmf_obj$elbo <- curr_elbo
      sym_ebcovmf_obj$vec_elbo_K[K] <- curr_elbo
    }
  }
  return(sym_ebcovmf_obj)
}

nullcheck_factors <- function(S, sym_ebcovmf_obj, L2_tol = 10^(-8)){
  null_lambda_idx <- which(sym_ebcovmf_obj$lambda == 0)
  factor_L2_norms <- apply(sym_ebcovmf_obj$L_pm, 2, function(v){sqrt(sum(v^2))})
  null_factor_idx <- which(factor_L2_norms < L2_tol)
  null_idx <- unique(c(null_lambda_idx, null_factor_idx))
  
  keep_idx <- setdiff(c(1:length(sym_ebcovmf_obj$lambda)), null_idx)
  
  if (length(keep_idx) < length(sym_ebcovmf_obj$lambda)){
    # remove factors
    sym_ebcovmf_obj$L_pm <- sym_ebcovmf_obj$L_pm[,keep_idx]
    sym_ebcovmf_obj$lambda <- sym_ebcovmf_obj$lambda[keep_idx]
    sym_ebcovmf_obj$KL <- sym_ebcovmf_obj$KL[keep_idx]
    sym_ebcovmf_obj$fitted_gs <- sym_ebcovmf_obj$fitted_gs[keep_idx]
  }
  
  K <- length(keep_idx)
  # recompute residual variance estimate
  sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                L = sym_ebcovmf_obj$L_pm,
                                                lambda = sym_ebcovmf_obj$lambda,
                                                n = sym_ebcovmf_obj$n,
                                                K = K)
  
  # recompute elbo
  curr_elbo <- compute_elbo(S = S,
                            L = sym_ebcovmf_obj$L_pm,
                            lambda = sym_ebcovmf_obj$lambda,
                            resid_s2 = sym_ebcovmf_obj$resid_s2,
                            n = sym_ebcovmf_obj$n,
                            K = K,
                            KL = sym_ebcovmf_obj$KL)
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_backfit_alt <- function(S, sym_ebcovmf_obj, ebnm_fn, backfit_maxiter = 100, backfit_tol = 10^(-8), optim_maxiter= 500, optim_tol = 10^(-8)){
  K <- length(sym_ebcovmf_obj$lambda)
  kset <- c(1:K)
  iter <- 1
  obj_diff <- Inf
  sym_ebcovmf_obj$backfit_vec_elbo_full <- NULL
  sym_ebcovmf_obj$backfit_iter_elbo_vec <- NULL
  
  # refit lambda
  # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 25)
  
  while((iter <= backfit_maxiter) && (obj_diff > backfit_tol)){
    # print(iter)
    obj_old <- sym_ebcovmf_obj$elbo
    # loop through each factor
    for (k in kset){
      # print(k)
      print(paste('Updating factor', k))
      # compute residual matrix
      R <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k, drop = FALSE] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
      R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm[,-k, drop = FALSE], sym_ebcovmf_obj$lambda[-k], (K-1)) #this is right but I have one instance where the values don't match what I expect
      
      # optimize factor
      factor_proposed <- optimize_factor(R, ebnm_fn, optim_maxiter, optim_tol, sym_ebcovmf_obj$L_pm[,k], sym_ebcovmf_obj$lambda[k], sym_ebcovmf_obj$fitted_gs[[k]], R2k, sym_ebcovmf_obj$n, sym_ebcovmf_obj$KL[-k])
      
      # update object
      # check if update leads to increase in objective function
      if ((factor_proposed$curr_elbo > sym_ebcovmf_obj$elbo) | (iter == 1)){
        sym_ebcovmf_obj$L_pm[,k] <- factor_proposed$v
        sym_ebcovmf_obj$KL[k] <- factor_proposed$rank_one_KL
        sym_ebcovmf_obj$lambda[k] <- factor_proposed$lambda_k
        sym_ebcovmf_obj$resid_s2 <- factor_proposed$resid_s2
        sym_ebcovmf_obj$fitted_gs[[k]] <- factor_proposed$fitted_g_k
        sym_ebcovmf_obj$elbo <- factor_proposed$curr_elbo
        sym_ebcovmf_obj$backfit_vec_elbo_full <- c(sym_ebcovmf_obj$backfit_vec_elbo_full, factor_proposed$vec_elbo_full)
      } else {
        obj_diff <- sym_ebcovmf_obj$elbo - factor_proposed$curr_elbo
        print(paste('update to factor', k, 'decreased the elbo by', abs(obj_diff)))
      }
      
      #print(sym_ebcovmf_obj$elbo)
      sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 1, remove_null = FALSE) # add refitting step?
      # print(sym_ebcovmf_obj$lambda)
      #print(sym_ebcovmf_obj$elbo)
    }
    kset <- which(sym_ebcovmf_obj$lambda != 0)
    sym_ebcovmf_obj$backfit_iter_elbo_vec <- c(sym_ebcovmf_obj$backfit_iter_elbo_vec, sym_ebcovmf_obj$elbo)
    
    iter <- iter + 1
    obj_diff <- abs(sym_ebcovmf_obj$elbo - obj_old)
  }
  # nullcheck
  sym_ebcovmf_obj <- nullcheck_factors(S, sym_ebcovmf_obj)
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_backfit_alt_v2 <- function(S, sym_ebcovmf_obj, ebnm_fn, backfit_maxiter = 100, backfit_tol = 10^(-8), optim_maxiter= 500, optim_tol = 10^(-8)){
  K <- length(sym_ebcovmf_obj$lambda)
  kset <- c(1:K)
  iter <- 1
  obj_diff <- Inf
  sym_ebcovmf_obj$backfit_vec_elbo_full <- NULL
  sym_ebcovmf_obj$backfit_iter_elbo_vec <- NULL
  
  # refit lambda
  # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 25)
  
  while((iter <= backfit_maxiter) && (obj_diff > backfit_tol)){
    # print(iter)
    obj_old <- sym_ebcovmf_obj$elbo
    # loop through each factor
    for (k in kset){
      if (sym_ebcovmf_obj$lambda[k] == 0){
        print(paste('Skipping factor', k, 'because lambda is 0'))
        next
      }
      # print(k)
      print(paste('Updating factor', k))
      # compute residual matrix
      R <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k, drop = FALSE] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
      R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm[,-k, drop = FALSE], sym_ebcovmf_obj$lambda[-k], (K-1)) #this is right but I have one instance where the values don't match what I expect
      
      # optimize factor
      factor_proposed <- optimize_factor(R, ebnm_fn, optim_maxiter, optim_tol, sym_ebcovmf_obj$L_pm[,k], sym_ebcovmf_obj$lambda[k], sym_ebcovmf_obj$fitted_gs[[k]], R2k, sym_ebcovmf_obj$n, sym_ebcovmf_obj$KL[-k])
      
      # update object
      # check if update leads to increase in objective function
      if ((factor_proposed$curr_elbo > sym_ebcovmf_obj$elbo) | (iter == 1)){
        sym_ebcovmf_obj$L_pm[,k] <- factor_proposed$v
        sym_ebcovmf_obj$KL[k] <- factor_proposed$rank_one_KL
        sym_ebcovmf_obj$lambda[k] <- factor_proposed$lambda_k
        sym_ebcovmf_obj$resid_s2 <- factor_proposed$resid_s2
        sym_ebcovmf_obj$fitted_gs[[k]] <- factor_proposed$fitted_g_k
        sym_ebcovmf_obj$elbo <- factor_proposed$curr_elbo
        sym_ebcovmf_obj$backfit_vec_elbo_full <- c(sym_ebcovmf_obj$backfit_vec_elbo_full, factor_proposed$vec_elbo_full)
      } else {
        obj_diff <- sym_ebcovmf_obj$elbo - factor_proposed$curr_elbo
        print(paste('update to factor', k, 'decreased the elbo by', abs(obj_diff)))
      }
      
      #print(sym_ebcovmf_obj$elbo)
      sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 1, remove_null = FALSE) # add refitting step?
      # print(sym_ebcovmf_obj$lambda)
      #print(sym_ebcovmf_obj$elbo)
    }
    # kset <- intersect(kset, which(sym_ebcovmf_obj$lambda != 0))
    kset <- which(sym_ebcovmf_obj$lambda != 0)
    sym_ebcovmf_obj$backfit_iter_elbo_vec <- c(sym_ebcovmf_obj$backfit_iter_elbo_vec, sym_ebcovmf_obj$elbo)
    
    iter <- iter + 1
    obj_diff <- abs(sym_ebcovmf_obj$elbo - obj_old)
  }
  # nullcheck
  sym_ebcovmf_obj <- nullcheck_factors(S, sym_ebcovmf_obj)
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_backfit_alt_v3 <- function(S, sym_ebcovmf_obj, ebnm_fn, backfit_maxiter = 100, backfit_tol = 10^(-8), optim_maxiter= 500, optim_tol = 10^(-8)){
  K <- length(sym_ebcovmf_obj$lambda)
  kset <- c(1:K)
  iter <- 1
  obj_diff <- Inf
  sym_ebcovmf_obj$backfit_vec_elbo_full <- NULL
  sym_ebcovmf_obj$backfit_iter_elbo_vec <- NULL
  
  # refit lambda
  # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 25)
  
  while((iter <= backfit_maxiter) && (obj_diff > backfit_tol)){
    # print(iter)
    obj_old <- sym_ebcovmf_obj$elbo
    # loop through each factor
    for (k in kset){
      R <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k, drop = FALSE] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
      
      if (sym_ebcovmf_obj$lambda[k] == 0 | as.numeric(t(sym_ebcovmf_obj$L_pm[,k, drop = FALSE]) %*% R %*% sym_ebcovmf_obj$L_pm[,k, drop = FALSE]) <= 0){
        print(paste('Skipping factor', k, 'because lambda is 0'))
        next
      }
      # print(k)
      print(paste('Updating factor', k))
      # compute residual matrix
      R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm[,-k, drop = FALSE], sym_ebcovmf_obj$lambda[-k], (K-1)) #this is right but I have one instance where the values don't match what I expect
      
      # optimize factor
      factor_proposed <- optimize_factor(R, ebnm_fn, optim_maxiter, optim_tol, sym_ebcovmf_obj$L_pm[,k], sym_ebcovmf_obj$lambda[k], sym_ebcovmf_obj$fitted_gs[[k]], R2k, sym_ebcovmf_obj$n, sym_ebcovmf_obj$KL[-k])
      
      # update object
      # check if update leads to increase in objective function
      if ((factor_proposed$curr_elbo > sym_ebcovmf_obj$elbo) | (iter == 1)){
        sym_ebcovmf_obj$L_pm[,k] <- factor_proposed$v
        sym_ebcovmf_obj$KL[k] <- factor_proposed$rank_one_KL
        sym_ebcovmf_obj$lambda[k] <- factor_proposed$lambda_k
        sym_ebcovmf_obj$resid_s2 <- factor_proposed$resid_s2
        sym_ebcovmf_obj$fitted_gs[[k]] <- factor_proposed$fitted_g_k
        sym_ebcovmf_obj$elbo <- factor_proposed$curr_elbo
        sym_ebcovmf_obj$backfit_vec_elbo_full <- c(sym_ebcovmf_obj$backfit_vec_elbo_full, factor_proposed$vec_elbo_full)
      } else {
        obj_diff <- sym_ebcovmf_obj$elbo - factor_proposed$curr_elbo
        print(paste('update to factor', k, 'decreased the elbo by', abs(obj_diff)))
      }
      
      #print(sym_ebcovmf_obj$elbo)
      sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 1, remove_null = FALSE) # add refitting step?
      # print(sym_ebcovmf_obj$lambda)
      #print(sym_ebcovmf_obj$elbo)
    }
    # kset <- intersect(kset, which(sym_ebcovmf_obj$lambda != 0))
    kset <- which(sym_ebcovmf_obj$lambda != 0)
    sym_ebcovmf_obj$backfit_iter_elbo_vec <- c(sym_ebcovmf_obj$backfit_iter_elbo_vec, sym_ebcovmf_obj$elbo)
    
    iter <- iter + 1
    obj_diff <- abs(sym_ebcovmf_obj$elbo - obj_old)
  }
  # nullcheck
  sym_ebcovmf_obj <- nullcheck_factors(S, sym_ebcovmf_obj)
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_backfit_alt_v4 <- function(S, sym_ebcovmf_obj, ebnm_fn, backfit_maxiter = 100, backfit_tol = 10^(-8), optim_maxiter= 500, optim_tol = 10^(-8)){
  K <- length(sym_ebcovmf_obj$lambda)
  kset <- c(1:K)
  iter <- 1
  obj_diff <- Inf
  sym_ebcovmf_obj$backfit_vec_elbo_full <- NULL
  sym_ebcovmf_obj$backfit_iter_elbo_vec <- NULL
  
  # refit lambda
  # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 25)
  
  while((iter <= backfit_maxiter) && (obj_diff > backfit_tol)){
    # print(iter)
    obj_old <- sym_ebcovmf_obj$elbo
    # loop through each factor
    for (k in kset){
      if (sym_ebcovmf_obj$lambda[k] == 0){
        print(paste('Skipping factor', k, 'because lambda is 0'))
        next
      }
      # print(k)
      print(paste('Updating factor', k))
      # compute residual matrix
      R <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k, drop = FALSE] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
      R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm[,-k, drop = FALSE], sym_ebcovmf_obj$lambda[-k], (K-1)) #this is right but I have one instance where the values don't match what I expect
      
      # optimize factor
      factor_proposed <- optimize_factor(R, ebnm_fn, optim_maxiter, optim_tol, sym_ebcovmf_obj$L_pm[,k], sym_ebcovmf_obj$lambda[k], sym_ebcovmf_obj$fitted_gs[[k]], R2k, sym_ebcovmf_obj$n, sym_ebcovmf_obj$KL[-k])
      
      # update object
      # check if update leads to increase in objective function
      if ((factor_proposed$curr_elbo > sym_ebcovmf_obj$elbo) | (iter == 1)){
        sym_ebcovmf_obj$L_pm[,k] <- factor_proposed$v
        sym_ebcovmf_obj$KL[k] <- factor_proposed$rank_one_KL
        sym_ebcovmf_obj$lambda[k] <- factor_proposed$lambda_k
        sym_ebcovmf_obj$resid_s2 <- factor_proposed$resid_s2
        sym_ebcovmf_obj$fitted_gs[[k]] <- factor_proposed$fitted_g_k
        sym_ebcovmf_obj$elbo <- factor_proposed$curr_elbo
        sym_ebcovmf_obj$backfit_vec_elbo_full <- c(sym_ebcovmf_obj$backfit_vec_elbo_full, factor_proposed$vec_elbo_full)
      } else {
        obj_diff <- sym_ebcovmf_obj$elbo - factor_proposed$curr_elbo
        print(paste('update to factor', k, 'decreased the elbo by', abs(obj_diff)))
      }
      
      #print(sym_ebcovmf_obj$elbo)
      sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 100, remove_null = FALSE) # add refitting step?
      # print(sym_ebcovmf_obj$lambda)
      #print(sym_ebcovmf_obj$elbo)
    }
    # kset <- intersect(kset, which(sym_ebcovmf_obj$lambda != 0))
    kset <- which(sym_ebcovmf_obj$lambda != 0)
    sym_ebcovmf_obj$backfit_iter_elbo_vec <- c(sym_ebcovmf_obj$backfit_iter_elbo_vec, sym_ebcovmf_obj$elbo)
    
    iter <- iter + 1
    obj_diff <- abs(sym_ebcovmf_obj$elbo - obj_old)
  }
  # nullcheck
  sym_ebcovmf_obj <- nullcheck_factors(S, sym_ebcovmf_obj)
  return(sym_ebcovmf_obj)
}
