renew_DR <- function(B, tempdatadir, family, p1, p2){
    

    tol <- 1e-6
    maxit <- 200
    time_load <- 0
    
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
    }

  #  J_accum <- t(S_accum) %*% solve(V_accum) %*% S_accum
  #  sd_theta <- sqrt(diag(solve(J_accum)))
    J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
    sd_theta <- sqrt(diag(J_accum_inv))

    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = theta_new / sd_theta)

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()

    out$theta <- out_theta
    out$time <- time_total
    out$run <- time_run
   
    return(out)
}

renew_RM <- function(B, tempdatadir, family, p1){

    tol <- 1e-6
    maxit <- 200
    time_load <- 0
    
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
    }

    J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
    sd_theta <- sqrt(diag(J_accum_inv))

    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = theta_new / sd_theta)

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()

    out$theta <- out_theta
    out$time <- time_total
    out$run <- time_run
   
    return(out)
}

renew_IW <- function(B, tempdatadir, family, p2){
    

    tol <- 1e-6
    maxit <- 200
    time_load <- 0
    
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
      
      W <- as.matrix(W)
      Z <- as.matrix(Z)
      Y <- as.matrix(Y)
    
      estimate <- increIW(W, Z, Y, family, theta_old, S_accum, V_accum, p2, maxit, tol)

      theta_new <- estimate$theta
      S_accum <- estimate$S_accum
      V_accum <- estimate$V_accum
    }

 #   J_accum <- t(S_accum) %*% solve(V_accum) %*% S_accum
 #   sd_theta <- sqrt(diag(solve(J_accum)))

    J_accum_inv <- solve(S_accum) %*% V_accum %*% solve(t(S_accum))
    sd_theta <- sqrt(diag(J_accum_inv))

    out_theta <- cbind(Estimate = drop(theta_new), StdErr = sd_theta, Zscore = theta_new / sd_theta)

    time_total <- (proc.time() - ptm)[3]
    time_run <- time_total - time_load
   
    out <- list()

    out$theta <- out_theta
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


