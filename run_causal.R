setwd("/Users/lanluo/Dropbox/Online Causal Inference2021/code")

#setwd("/Users/lluo10/OneDrive - University of Iowa/Online causal inference/Simulation")

library(MASS)
library(Rcpp)

source("datagenerator.R")
source("offline_ATE.R")
source("renew_ATE.R")
sourceCpp("increATE.cpp")


################

nsim <- 500

N <- 10000
B <- 50
n <- N / B   # sample size in every single data batch in block k

interaction <- TRUE


p2 <- 4
if(interaction == TRUE){
    p1 <- p2 + 2
 
    alpha <- c(-0.63, 0.18, -0.84, 0.59, 0.32, -0.82)
    }else{
    p1 <- p2 + 1
    alpha <- c(-0.63, 0.18, -0.84, 0.59, 0.32)
    }
d <- p1 + p2 + 1

family <- "binomial"  # family type for the regression model  "binomial"

beta <- c(-0.82, 0.49, -0.73, 0.58)

rho_x <- 0.5
outputfilename <- paste(family, "_", "N", N, "_", "B", B, "_", "d", d, sep = "")

# approximate ATE
if(family == "gaussian" && interaction == FALSE){
    delta_ref <- alpha[2]} else {

delta_vec <- c()
for(s in 1 : 2000){
    print(s)
    tempdatadir <- paste("Temp_", outputfilename, "_", s, sep = "")

    datagenerator_fix(alpha = alpha, beta = beta, n = n, p1 = p1, p2 = p2, B = B, family = family,
        construct = "cs", rho = rho_x, tempdatadir = tempdatadir, interaction, seed = s)
    load("delta.RData")
    delta_vec <- c(delta_vec, delta)
    unlink(tempdatadir, recursive = TRUE)
}
delta_ref <- mean(delta_vec)

}

#delta_ref <- 0.03782009
#delta_ref <- -0.05033158



for(s in 1:nsim){
    
    print(s)

    tempdatadir <- paste("Temp_", outputfilename, "_", s, sep = "")
    
    # Generate B data batches
    datagenerator_fix(alpha = alpha, beta = beta, n = n, p1 = p1, p2 = p2, B = B, family = family,
        construct = "cs", rho = rho_x, tempdatadir = tempdatadir, interaction, seed = s)

    fulldata(B = B, tempdatadir = tempdatadir)

    result.A.RM <- offline_RM(tempdatadir, family, p1 = p1)
    print("offline causal RM: done!")
    
    result.A.IW <- offline_IW(tempdatadir, family, p2 = p2)
    print("offline causal IW: done!")

    result.A.DR <- offline_DR(tempdatadir, family, p1 = p1, p2 = p2)
    print("offline causal DR: done!")

    result.B.RM <- renew_RM(B, tempdatadir, family, p1 = p1)
    print("online causal RM: done!")

    result.B.IW <- renew_IW(B, tempdatadir, family, p2 = p2)
    print("online causal IW: done!")

    result.B.DR <- renew_DR(B, tempdatadir, family, p1 = p1, p2 = p2)
    print("online causal DR: done!")
    
    # offline estimation results
    A.theta.est.RM <- result.A.RM$theta[, 1]
    A.theta.sd.RM <- result.A.RM$theta[, 2]

    A.theta.est.IW <- result.A.IW$theta[, 1]
    A.theta.sd.IW <- result.A.IW$theta[, 2]

    A.theta.est.DR <- result.A.DR$theta[, 1]
    A.theta.sd.DR <- result.A.DR$theta[, 2]

    # online estimation results
    B.theta.est.RM <- result.B.RM$theta[, 1]
    B.theta.sd.RM <- result.B.RM$theta[, 2]

    B.theta.est.IW <- result.B.IW$theta[, 1]
    B.theta.sd.IW <- result.B.IW$theta[, 2]

    B.theta.est.DR <- result.B.DR$theta[, 1]
    B.theta.sd.DR <- result.B.DR$theta[, 2]

    unlink(tempdatadir, recursive = TRUE)
    

 #   out <- c(eval.func(c(alpha, theta[d]), theta_est_RM, theta_sd_RM), result_RM$time, result_RM$run,
 #       eval.func(c(beta, theta[d]), theta_est_IW, theta_sd_IW), result_IW$time, result_IW$run,
 #       eval.func(c(alpha, beta, theta[d]), theta_est_DR, theta_sd_DR), result_DR$time, result_DR$run)
    
 # other parameters exclude ATE
 #    out <- c(eval.func(alpha, theta_est_RM[1:p1], theta_sd_RM[1:p1]), result_RM$time, result_RM$run,
 #       eval.func(beta, theta_est_IW[1:p2], theta_sd_IW[1:p2]), result_IW$time, result_IW$run,
 #       eval.func(theta[1:d-1], theta_est_DR[1:d-1], theta_sd_DR[1:d-1]), result_DR$time, result_DR$run)
 # ATE 

     out_offline <- c(
        eval.func(delta_ref, A.theta.est.RM[p1+1], A.theta.sd.RM[p1+1]), result.A.RM$time, result.A.RM$run,
        eval.func(delta_ref, A.theta.est.IW[p2+1], A.theta.sd.IW[p2+1]), result.A.IW$time, result.A.IW$run,
        eval.func(delta_ref, A.theta.est.DR[d], A.theta.sd.DR[d]), result.A.DR$time, result.A.DR$run
        )

     out_online <- c(
        eval.func(delta_ref, B.theta.est.RM[p1+1], B.theta.sd.RM[p1+1]), result.B.RM$time, result.B.RM$run,
        eval.func(delta_ref, B.theta.est.IW[p2+1], B.theta.sd.IW[p2+1]), result.B.IW$time, result.B.IW$run,
        eval.func(delta_ref, B.theta.est.DR[d], B.theta.sd.DR[d]), result.B.DR$time, result.B.DR$run
        )
    
 #   out <- c(eval.func(theta[1], theta_est_RM[1], theta_sd_RM[1]), result_RM$time, result_RM$run,
 #       eval.func(theta[1], theta_est_IW[1], theta_sd_IW[1]), result_IW$time, result_IW$run,
 #       eval.func(theta[1], theta_est_DR[1], theta_sd_DR[1]), result_DR$time, result_DR$run)
    

    out_ese_offline <- c(A.theta.est.RM, A.theta.est.IW, A.theta.est.DR)
    out_ese_online <- c(B.theta.est.RM, B.theta.est.IW, B.theta.est.DR)


    write.table(as.matrix(t(out_offline)), file=paste(outputfilename, "offline.csv", sep=""), sep=",",
        col.names=FALSE, row.names=s, append=TRUE)
    write.table(as.matrix(t(out_ese_offline)), file=paste(outputfilename, "offline.ese.csv", sep=""), sep=",",
      col.names=FALSE, row.names=s, append=TRUE)

    write.table(as.matrix(t(out_online)), file=paste(outputfilename, "online.csv", sep=""), sep=",",
        col.names=FALSE, row.names=s, append=TRUE)
    write.table(as.matrix(t(out_ese_online)), file=paste(outputfilename, "online.ese.csv", sep=""), sep=",",
      col.names=FALSE, row.names=s, append=TRUE)
}


result <- read.table(paste(outputfilename, "offline.csv", sep = ""), sep = ",")[, -1]
colnames(result) <- c(
    "A.theta_RM_bias", "A.theta_RM_sd", "A.theta_RM_cp", "A.t_RM.time", "A.r_RM.time",
    "A.theta_IW_bias", "A.theta_IW_sd", "A.theta_IW_cp", "A.t_IW.time", "A.r_IW.time",
	"A.theta_DR_bias", "A.theta_DR_sd", "A.theta_DR_cp", "A.t_DR.time", "A.r_DR.time")
result_summary <- apply(result, 2, mean)
print(result_summary)
    

ese <- read.table(paste(outputfilename, "offline.ese.csv", sep = ""), sep = ",")[, -1]
bias_offline <- apply(ese, 2, mean)[c(p1 + 1, p1 + 1 + p2 +1, p1 + p2 + 2 + d)] - delta_ref
bias_offline
bias_offline / delta_ref

ese_col <- apply(ese, 2, sd)

#mean(ese_col[1 : p1])
#mean(ese_col[(p1 + 2): (p1 + p2 + 1)])
#mean(ese_col[(p1 + p2 + 2):(p1 + p2 +2 +d)])

ese_col[p1 + 1]

ese_col[p1 + 1 + p2 +1]

ese_col[p1 + 1 +p2 + 1 + d]

######## online 
result <- read.table(paste(outputfilename, "online.csv", sep = ""), sep = ",")[, -1]
colnames(result) <- c(
    "B.theta_RM_bias", "B.theta_RM_sd", "B.theta_RM_cp", "B.t_RM.time", "B.r_RM.time",
    "B.theta_IW_bias", "B.theta_IW_sd", "B.theta_IW_cp", "B.t_IW.time", "B.r_IW.time",
    "B.theta_DR_bias", "B.theta_DR_sd", "B.theta_DR_cp", "B.t_DR.time", "B.r_DR.time")
result_summary <- apply(result, 2, mean)
print(result_summary)
    

ese <- read.table(paste(outputfilename, "online.ese.csv", sep = ""), sep = ",")[, -1]

bias_online <- apply(ese, 2, mean)[c(p1 + 1, p1 + 1 + p2 +1, p1 + p2 + 2 + d)] - delta_ref
bias_online
bias_online / delta_ref


ese_col <- apply(ese, 2, sd)

ese_col[p1 + 1]

ese_col[p1 + 1 + p2 +1]

ese_col[p1 + 1 + p2 + 1 + d]





