compute_elbo <- function(S = NULL, L = NULL, lambda = NULL, resid_s2, n, K = NULL, KL, R2 = NULL){
  if (is.null(R2) == TRUE){
    EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
    elbo <- -((n*(n-1)/4)*log(2*pi*resid_s2)) - 
      (n/2)*log(2*pi*2*resid_s2) -
      (1/(4*resid_s2))*(sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
                          sum(lambda^2) - sum(EvTEv^2 * lambda^2)) + sum(KL)
  } else {
    elbo <- -((n*(n-1)/4)*log(2*pi*resid_s2)) - (n/2)*log(2*pi*2*resid_s2) -
      (1/(4*resid_s2))*R2 + sum(KL)
  }
  return(elbo)
}

compute_R2 <- function(S, L, lambda, K){
  EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
  R2 <- sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
    sum(lambda^2) - sum(EvTEv^2 * lambda^2)
  return(R2)
}

estimate_resid_s2 <- function(S = NULL, L = NULL, lambda = NULL, n, K = NULL, R2 = NULL){
  if (is.null(R2) == TRUE){
    EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
    obj_func_fit <- sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
      sum(lambda^2) - sum(EvTEv^2 * lambda^2)
    resid_s2 <- obj_func_fit/(n*(n+1))
  } else {
    resid_s2 <- R2/(n*(n+1))
  }
  return(resid_s2)
}

normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}

refit_lambda <- function(S, sym_ebcovmf_obj, maxiter = 100, tol = 10^(-6)){
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
      }
      
      # Update resid_s2
      resid_s2.old <- sym_ebcovmf_obj$resid_s2
      sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                    L = sym_ebcovmf_obj$L_pm,
                                                    lambda = sym_ebcovmf_obj$lambda,
                                                    n = sym_ebcovmf_obj$n,
                                                    K = K)
      
      # Update elbo
      curr_elbo.old <- curr_elbo
      curr_elbo <- compute_elbo(S = S,
                                L = sym_ebcovmf_obj$L_pm,
                                lambda = sym_ebcovmf_obj$lambda,
                                resid_s2 = sym_ebcovmf_obj$resid_s2,
                                n = sym_ebcovmf_obj$n,
                                K = K,
                                KL = sym_ebcovmf_obj$KL)
      
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
    if (any(sym_ebcovmf_obj$lambda == 0)){
      idx <- which(sym_ebcovmf_obj$lambda != 0)
      sym_ebcovmf_obj$lambda <- sym_ebcovmf_obj$lambda[idx]
      sym_ebcovmf_obj$L_pm <- sym_ebcovmf_obj$L_pm[,idx]
      sym_ebcovmf_obj$KL <- sym_ebcovmf_obj$KL[idx]
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
    if ((curr_elbo - sym_ebcovmf_obj.old$elbo) < 0){
      sym_ebcovmf_obj <- sym_ebcovmf_obj.old
    } else {
      sym_ebcovmf_obj$elbo <- curr_elbo
      sym_ebcovmf_obj$vec_elbo_K[K] <- curr_elbo
    }
  }
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_fit <- function(S, ebnm_fn, Kmax, maxiter, rank_one_tol, tol, refit_lam = FALSE){
  #initialize object
  sym_ebcovmf_obj <- sym_ebcovmf_init(S)
  
  curr_rank <- 0
  obj_diff <- Inf
  while ((curr_rank < Kmax) & (obj_diff > tol)){
    # add factor
    sym_ebcovmf_obj <- sym_ebcovmf_r1_fit(S, sym_ebcovmf_obj, ebnm_fn, maxiter, rank_one_tol)
    
    # check if new factor was added
    if (length(sym_ebcovmf_obj$vec_elbo_K) == curr_rank){
      print(paste('Adding factor', (curr_rank + 1), 'does not improve the objective function'))
      break
    } else {
      if (curr_rank > 0){
        if (refit_lam == TRUE){
          sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj)
        }
        obj_diff <- sym_ebcovmf_obj$vec_elbo_K[curr_rank + 1] - sym_ebcovmf_obj$vec_elbo_K[curr_rank]
      }
    }
    curr_rank <- curr_rank + 1
  }
  
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_init <- function(S){
  sym_ebcovmf_obj <- list(n = ncol(S), 
                          L_pm = NULL, 
                          resid_s2 = 0, 
                          KL = c(), 
                          lambda = c(), 
                          elbo = 0, 
                          vec_elbo_K = c(), 
                          vec_elbo_full = c(), 
                          fitted_gs = list())
  return(sym_ebcovmf_obj)
}


sym_ebcovmf_r1_fit <- function(S, sym_ebcovmf_obj, ebnm_fn, maxiter, tol, v_init = NULL){
  if (is.null(sym_ebcovmf_obj$L_pm) == FALSE){
    K <- length(sym_ebcovmf_obj$lambda) + 1
    R <- S - tcrossprod(sym_ebcovmf_obj$L_pm %*% diag(sqrt(sym_ebcovmf_obj$lambda), ncol = (K-1)))
    R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm, sym_ebcovmf_obj$lambda, (K-1))
  } else {
    K <- 1
    R <- S
    R2k <- sum(S^2)
  }
  sym_ebcovmf_obj.old <- sym_ebcovmf_obj
  
  # initialize estimate for l
  if (is.null(v_init) == TRUE){
    # need to add something that checks non-negativity of prior
    sym_ebcovmf_v_init <- sym_ebcovmf_r1_init(R)
    v <- sym_ebcovmf_v_init$v
    lambda_k <- sym_ebcovmf_v_init$lambda_k
  } else {
    v <- v_init
    v <- v/sqrt(sum(v^2))
    lambda_k <- drop(t(v) %*% R %*% v)
  }
  
  # initialize other values
  R2 <- R2k - lambda_k^2
  resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)
  rank_one_KL <- 0
  curr_elbo <- -Inf
  obj_diff <- Inf
  fitted_g_k <- NULL
  iter <- 1
  
  sym_ebcovmf_obj$vec_elbo_full <- c(sym_ebcovmf_obj$vec_elbo_full, K)
  while((iter <= maxiter) && (obj_diff > tol)){
    # update l; power iteration step
    v.old <- v
    x <- R %*% v
    e <- ebnm_fn(x = x, s = sqrt(resid_s2), g_init = fitted_g_k)
    scaling_factor <- sqrt(sum(e$posterior$mean^2) + sum(e$posterior$sd^2))
    if (scaling_factor == 0){ # check if scaling factor is zero
      scaling_factor <- Inf
      v <- e$posterior$mean/scaling_factor
      print('Warning: scaling factor is zero')
      break
    }
    v <- e$posterior$mean/scaling_factor
    
    # update lambda and R2
    lambda_k.old <- lambda_k
    lambda_k <- max(as.numeric(t(v) %*% R %*% v), 0)
    R2 <- R2k - lambda_k^2
    
    #store estimate for g
    fitted_g_k.old <- fitted_g_k
    fitted_g_k <- e$fitted_g
    
    # store KL
    rank_one_KL.old <- rank_one_KL
    rank_one_KL <- as.numeric(e$log_likelihood) +
      - normal_means_loglik(x, sqrt(resid_s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)
    
    # update resid_s2
    resid_s2.old <- resid_s2
    resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)
    
    # check convergence
    curr_elbo.old <- curr_elbo
    curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                              n = sym_ebcovmf_obj$n,
                              KL = c(sym_ebcovmf_obj$KL, rank_one_KL),
                              R2 = R2)
    if (iter > 1){
      obj_diff <- curr_elbo - curr_elbo.old
    }
    if (obj_diff < 0){ # check if convergence_val < 0
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('elbo decreased by', abs(obj_diff)))
      break
    }
    sym_ebcovmf_obj$vec_elbo_full <- c(sym_ebcovmf_obj$vec_elbo_full, curr_elbo)
    iter <- iter + 1
  }
  
  # nullcheck
  if((lambda_k == 0) | (sqrt(sum(v^2)) < 10^(-8))){
    #print('additional factor does not improve fit')
    sym_ebcovmf_obj <- sym_ebcovmf_obj.old
  } else {
    sym_ebcovmf_obj$L_pm <- cbind(sym_ebcovmf_obj$L_pm, v)
    sym_ebcovmf_obj$KL[K] <- rank_one_KL
    sym_ebcovmf_obj$lambda[K] <- lambda_k
    sym_ebcovmf_obj$resid_s2 <- resid_s2
    sym_ebcovmf_obj$fitted_gs[[K]] <- fitted_g_k
    sym_ebcovmf_obj$elbo <- curr_elbo
    sym_ebcovmf_obj$vec_elbo_K <- c(sym_ebcovmf_obj$vec_elbo_K, curr_elbo)
  }
  return(sym_ebcovmf_obj)
}

sym_ebcovmf_r1_init <- function(R, nonnegative = TRUE){
  svd1 <- RSpectra::eigs_sym(R, k = 1)
  v <- svd1$vectors # scaled such that v'v = 1
  # if (abs(min(v)) > abs(max(v))){
  #   v <- -1*v
  # }
  lambda_k <- svd1$values
  if(nonnegative == TRUE){
    svd_v <- v
    v <- pmax(svd_v, 0)
    if (sqrt(sum(v^2)) > 0){
      v <- v/sqrt(sum(v^2))
    }
    
    minus_v <- pmax(-svd_v, 0)
    if(sqrt(sum(minus_v^2)) > 0){
      minus_v <- minus_v/sqrt(sum(minus_v^2))
    }
  }
  lambda_options <- c(drop(t(v) %*% R %*% v), drop(t(minus_v) %*% R %*% minus_v))
  if(which.max(lambda_options) == 1){
    v <- v
    lambda_k <- lambda_options[1]
  } else {
    v <- minus_v
    lambda_k <- lambda_options[2]
  }
  return(list(v = v, lambda_k = lambda_k))
}