renew_DR <- function(B, tempdatadir, family, p1, p2, delta0, test_freq, cutoffs){
    

    tol <- 1e-6
    maxit <- 100
    time_load <- 0
    decision <- 0
    
    decision_vec <- c()
    delta_vec <- c()
    sd_vec <- c()

    d <- p1 + p2 + 1

    theta_new <- rep(0, d)

    S_accum <- matrix(rep(0, d * d), nrow = d)
    V_accum <- matrix(rep(0, d * d), nrow = d)
    
    ptm <- proc.time()

    for(b in 1:B){

    	load <- proc.time()
    	load(paste(tempdatadir, "/", b, ".RData", sep = ""))
    	time_load <- time_load + (proc.time() - load)[3]

    	theta_old <- theta_new
      
      X <- as.matrix(X)
      W <- as.matrix(W)
      Z <- as.matrix(Z)
      Y <- as.matrix(Y)
    
      estimate <- increDR(X, W, Z, Y, family, theta_old, S_accum, V_accum, p1, p2, maxit, tol)

      theta_new <- estimate$theta
      S_accum <- estimate$S_accum
      V_accum <- estimate$V_accum

      J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
      sd_theta <- sqrt(diag(J_accum_inv))
      
      # sequential testing 
      if(b %% test_freq == 0){
        z_stat <- abs(theta_new[d] - delta0) / sd_theta[d]

          if(decision == 0){ # continue testing only if not rejected at previous step
            if(z_stat >= cutoffs[b / test_freq]){ # reject, if rejected, should stop
                decision <- 1
            } else {
                decision <- 0
            }
          }
      }
      decision_vec <- c(decision_vec, decision)
      delta_vec <- c(delta_vec, theta_new[d])
      sd_vec <- c(sd_vec, sd_theta[d])      
    }
    zvec <- delta_vec / sd_vec
    pvec <- pnorm(-abs(zvec))
    out_delta <- cbind(delta_vec, sd_vec, zvec = zvec, pvec = pvec)

    zscore = theta_new / sd_theta
    pvalue <- pnorm(-abs(zscore))

    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = zscore, pvalue = pvalue)

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()

    out$delta <- out_delta
    out$theta <- out_theta
    out$decision <- decision_vec
    out$time <- time_total
    out$run <- time_run
   
    return(out)
}

renew_RM <- function(B, tempdatadir, family, p1, delta0, test_freq, cutoffs){

    tol <- 1e-6
    maxit <- 100
    time_load <- 0
    decision <- 0

    decision_vec <- c()
    delta_vec <- c()
    sd_vec <- c()

    d <- p1 + 1

    theta_new <- rep(0, d)

    S_accum <- matrix(rep(0, d * d), nrow = d)
    V_accum <- matrix(rep(0, d * d), nrow = d)
    
    ptm <- proc.time()

    for(b in 1:B){

        load <- proc.time()
        load(paste(tempdatadir, "/", b, ".RData", sep = ""))
        time_load <- time_load + (proc.time() - load)[3]

        theta_old <- theta_new

        X <- as.matrix(X)
        Z <- as.matrix(Z)
        Y <- as.matrix(Y)
   
      estimate <- increRM(X, Z, Y, family, theta_old, S_accum, V_accum, p1, maxit, tol)

      theta_new <- estimate$theta
      S_accum <- estimate$S_accum
      V_accum <- estimate$V_accum

      J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
      sd_theta <- sqrt(diag(J_accum_inv))

      # sequential testing 
      if(b %% test_freq == 0){
        z_stat <- abs(theta_new[d] - delta0) / sd_theta[d]
  
          if(decision == 0){ # continue testing only if not rejected at previous step
            if(z_stat >= cutoffs[b / test_freq]){ # reject, if rejected, should stop
                decision <- 1
            } else {
                decision <- 0
            }
          }
      }
      decision_vec <- c(decision_vec, decision)
      delta_vec <- c(delta_vec, theta_new[d])
      sd_vec <- c(sd_vec, sd_theta[d])      
    }

    zvec <- delta_vec / sd_vec
    pvec <- pnorm(-abs(zvec))
    out_delta <- cbind(delta_vec, sd_vec, zvec = zvec, pvec = pvec)

    zscore = theta_new / sd_theta
    pvalue <- pnorm(-abs(zscore))
    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = zscore, pvalue = pvalue)
    

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()
    
    out$delta <- out_delta
    out$theta <- out_theta
    out$decision <- decision_vec
    out$time <- time_total
    out$run <- time_run
   
    return(out)
}

renew_IW <- function(B, tempdatadir, family, p2, delta0, test_freq, cutoffs){
    

    tol <- 1e-6
    maxit <- 100
    time_load <- 0
    decision <- 0

    decision_vec <- c()
    delta_vec <- c()
    sd_vec <- c()
    
    d <- p2 + 1

    theta_new <- rep(0, d)

    S_accum <- matrix(rep(0, d * d), nrow = d)
    V_accum <- matrix(rep(0, d * d), nrow = d)
    
    ptm <- proc.time()

    for(b in 1:B){

        load <- proc.time()
        load(paste(tempdatadir, "/", b, ".RData", sep = ""))
        time_load <- time_load + (proc.time() - load)[3]

        theta_old <- theta_new
      
      X <- as.matrix(X)
      W <- as.matrix(W)
      Z <- as.matrix(Z)
      Y <- as.matrix(Y)
    
      estimate <- increIW(W, Z, Y, family, theta_old, S_accum, V_accum, p2, maxit, tol)

      theta_new <- estimate$theta
      S_accum <- estimate$S_accum
      V_accum <- estimate$V_accum

      J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
      sd_theta <- sqrt(diag(J_accum_inv))
      
      # sequential testing 
      if(b %% test_freq == 0){
        z_stat <- abs(theta_new[d] - delta0) / sd_theta[d]
  
          if(decision == 0){ # continue testing only if not rejected at previous step
            if(z_stat >= cutoffs[b / test_freq]){ # reject, if rejected, should stop
                decision <- 1
            } else {
                decision <- 0
            }
          }
      }
      decision_vec <- c(decision_vec, decision)
      delta_vec <- c(delta_vec, theta_new[d])
      sd_vec <- c(sd_vec, sd_theta[d])      
    }
    
    zvec <- delta_vec / sd_vec
    pvec <- pnorm(-abs(zvec))
    out_delta <- cbind(delta_vec, sd_vec, zvec = zvec, pvec = pvec)

    zscore = theta_new / sd_theta
    pvalue <- pnorm(-abs(zscore))

    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = zscore, pvalue = pvalue)

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()
    
    out$delta <- out_delta
    out$theta <- out_theta
    out$decision <- decision_vec
    out$time <- time_total
    out$run <- time_run
   
    return(out)
}

eval.func <- function(theta, thetahat, sd){

    bias <- abs(thetahat - theta)
    pvalue <- 2 * pnorm(- bias / sd)
    
    covprob <- ifelse(pvalue >= 0.05, 1, 0)

    c(bias = mean(bias), sd = mean(sd), covprob = mean(covprob))
}


