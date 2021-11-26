datagenerator_fix <- function(alpha, beta, n, p1, p2, B, family, construct, rho, 
  tempdatadir, interaction = FALSE, seed = NA){

  alpha       # true parameter vector 1
  beta        # true parameter vector 2
  n           # sample size in each data batch
  p1          # dimension of alpha
  p2          # dimension of beta
  B           # number of data batches
  family      # c("gaussian", "binomial")
  construct   # c("ind", "cs", "ar1") correlation structure of X
  rho         # rho is associated with corstruct
  tempdatadir # name of the director to store simulated data
  seed        # set random seed
  

  if(length(n)!=1 & length(n)!= B){ stop("n must be a number of a vector of length B.") }
  if(!is.na(seed)){ set.seed(seed) }
#  if(length(n)==1){ n <- rep(n, B) } #same sample size
  seed.list <- sample(1:1e8, 1, replace = FALSE) # seed is for sampling the seed list for a sequence of batches
  N <- n * B

  dir.create(tempdatadir)
    #for a fixed total sample size, the dataset remains the same
    set.seed(seed.list)
    if(construct =='ind'){Sigma <- diag(rep(1, p2))}
    if(construct =='cs'){Sigma <- matrix(rho, p2, p2); diag(Sigma) <- 1}
    if(construct =='ar1'){Sigma <- rho^(abs(outer(1 : p2, 1 : p2, "-")))}

  Xp <- mvrnorm(n = N, mu = rep(0, p2), Sigma = Sigma)  # N * p2, covariate matrix for propensity score
  Xp[, 2] <- rbinom(n = N, size = 1, p = 0.5) # interaction term is binary
  Xp[, 1] <- 1 

  # for generating group assignment based on propensity score model
  zc <- rbinom(n = N, size = 1, p = exp(Xp %*% beta) / (1 + exp(Xp %*% beta)))
  
  Xp <- Xp[, -1] # remove intercept for the stored matrix


  if (interaction == TRUE){
    Xc <- cbind(Xp, zc * Xp[, 2])  # intercept not included
    Xz1 <- cbind(1, 1, Xp, 1 * Xp[, 2]) # design matrix for generating outcome in treatment group
    Xz0 <- cbind(1, 0, Xp, 0)
  } else{
    Xc <- Xp  # intercept not included
    Xz1 <- cbind(1, 1, Xp) # design matrix for generating outcome in treatment group
    Xz0 <- cbind(1, 0, Xp)
  }

  set.seed(seed.list + 1234567)
  if(family == "gaussian"){
    yc1 <- rnorm(n = N, mean = Xz1 %*% alpha) # point-wise multiplication
    yc0 <- rnorm(n = N, mean = Xz0 %*% alpha) # point-wise multiplication
  } else if (family == "binomial"){
    yc1 <- rbinom(n = N, size = 1, prob = exp(Xz1 %*% alpha) / (1 + exp(Xz1 %*% alpha)))
    yc0 <- rbinom(n = N, size = 1, prob = exp(Xz0 %*% alpha) / (1 + exp(Xz0 %*% alpha)))
  } else {stop("Unknown family distribution.")}
  
  yc <- yc1 * zc + yc0 * (1 - zc) # assign outcome based on treatment group assignment vec zc

  delta <- mean(yc1 - yc0) #~alpha[2] in linear regression model
 
  yc <- drop(yc)
  
 
    # store N samples in a total of B data batches
    for(b in 1 : B){
        Y <- yc[((b -1) * n +1) : (b * n)]
        Z <- zc[((b -1) * n +1) : (b * n)]
        X <- Xc[((b -1) * n +1) : (b * n), ] # include main effect and interaction (if exists)
        W <- Xp[((b -1) * n +1) : (b * n), ] # include only main effects in X
        save(Y, Z, X, W, file = paste(tempdatadir, "/", b, ".RData", sep = ""))
    }
  save(delta, file = "delta.RData")
}

fulldata <- function(B, tempdatadir){
  X.full <- c() 
  Y.full <- c()
  Z.full <- c() 
  W.full <- c()
  for(b in 1 : B){
    load(paste(tempdatadir, "/", b, ".RData", sep=""))
    X.full <- rbind(X.full, X)
    Y.full <- c(Y.full, Y)
    Z.full <- c(Z.full, Z) 
    W.full <- rbind(W.full, W)
  }
  X <- X.full
  W <- W.full
  Y <- Y.full
  Z <- Z.full
  save(Y, Z, X, W, file = paste(tempdatadir, "/", "full.RData", sep = ""))
}


