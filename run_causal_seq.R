setwd("/Users/lanluo/Dropbox/Online Causal Inference2021/code")

#setwd("/Users/lluo10/OneDrive - University of Iowa/Online causal inference/Simulation")

library(MASS)
library(Rcpp)
library(ldbounds)

source("datagenerator.R")
source("offline_ATE.R")
source("renew_ATE_seq.R")
sourceCpp("increATE_no_interaction.cpp")


#################

nsim <- 500

N <- 10000
B <- 10
n <- N / B   # sample size in every single data batch in block k

interaction <- FALSE


p2 <- 4
if(interaction == TRUE){
    p1 <- p2 + 2
    alpha <- c(-0.63, 0.18, -0.84, 0.59, 0.32, -0.82)
    }else{
    p1 <- p2 + 1
    alpha <- c(-0.63, 0.18, -0.84, 0.59, 0.32)
    }
d <- p1 + p2 + 1

family <- "gaussian"  # family type for the regression model  "binomial"


alpha[2] <- 0.12 # true delta under H0 or H1
beta <- c(-0.82, 0.49, -0.73, 0.58)


rho_x <- 0.5

############ for hypothesis testing ############
delta0 <- 0
test_freq <- 2
cutoffs <- bounds(seq(1, B, test_freq) / B, alpha = 0.05)$upper.bounds
###########

outputfilename <- paste(family, "_", "N", N, "_", "B", B, "_", "d", d, "_", "delta", alpha[2], sep = "")

# approximate ATE
if(family == "gaussian" && interaction == FALSE){
    delta_ref <- alpha[2]} else {

delta_vec <- c()
for(s in 1 : 1000){
    print(s)
    tempdatadir <- paste("Temp_", outputfilename, "_", s, sep = "")

    datagenerator_fix(alpha = alpha, beta = beta, n = n, p1 = p1, p2 = p2, B = B, family = family, construct = "cs", 
        rho = rho_x, tempdatadir = tempdatadir, interaction, seed = s)
    load("delta.RData")
    delta_vec <- c(delta_vec, delta)
    unlink(tempdatadir, recursive = TRUE)
}
delta_ref <- mean(delta_vec)

}


for(s in 1 : nsim){
    
    print(s)

    tempdatadir <- paste("Temp_", outputfilename, "_", s, sep = "")
    
    # Generate B data batches
    datagenerator_fix(alpha = alpha, beta = beta, n = n, p1 = p1, p2 = p2, B = B, family = family, construct = "cs", 
        rho = rho_x, tempdatadir = tempdatadir, interaction, seed = s)

    result.B.RM <- renew_RM(B, tempdatadir, family, p1 = p1, delta0, test_freq, cutoffs)
    print("online causal RM: done!")

    result.B.IW <- renew_IW(B, tempdatadir, family, p2 = p2, delta0, test_freq, cutoffs)
    print("online causal IW: done!")

    result.B.DR <- renew_DR(B, tempdatadir, family, p1 = p1, p2 = p2, delta0,
        test_freq, cutoffs)
    print("online causal DR: done!")
    
    # online estimation results
    B.theta.est.RM <- result.B.RM$theta[, 1]
    B.theta.sd.RM <- result.B.RM$theta[, 2]

    B.theta.est.IW <- result.B.IW$theta[, 1]
    B.theta.sd.IW <- result.B.IW$theta[, 2]

    B.theta.est.DR <- result.B.DR$theta[, 1]
    B.theta.sd.DR <- result.B.DR$theta[, 2]
    
    RM.seq.test <- result.B.RM$decision
    IW.seq.test <- result.B.IW$decision
    DR.seq.test <- result.B.DR$decision

    unlink(tempdatadir, recursive = TRUE)


     out_online <- c(
        eval.func(delta_ref, B.theta.est.RM[p1+1], B.theta.sd.RM[p1+1]), result.B.RM$time, result.B.RM$run,
        eval.func(delta_ref, B.theta.est.IW[p2+1], B.theta.sd.IW[p2+1]), result.B.IW$time, result.B.IW$run,
        eval.func(delta_ref, B.theta.est.DR[d], B.theta.sd.DR[d]), result.B.DR$time, result.B.DR$run,
        RM.seq.test[B], IW.seq.test[B], DR.seq.test[B]
        )
    
    out_ese_online <- c(B.theta.est.RM, B.theta.est.IW, B.theta.est.DR)


    write.table(as.matrix(t(out_online)), file=paste(outputfilename, "online.csv", sep = ""), sep=",",
        col.names=FALSE, row.names=s, append=TRUE)
    write.table(as.matrix(t(out_ese_online)), file=paste(outputfilename, "online.ese.csv", sep=""), sep=",",
      col.names=FALSE, row.names=s, append=TRUE)
}

######## online 
result <- read.table(paste(outputfilename, "online.csv", sep = ""), sep = ",")[, -1]
colnames(result) <- c(
    "B.theta_RM_bias", "B.theta_RM_sd", "B.theta_RM_cp", "B.t_RM.time", "B.r_RM.time",
    "B.theta_IW_bias", "B.theta_IW_sd", "B.theta_IW_cp", "B.t_IW.time", "B.r_IW.time",
    "B.theta_DR_bias", "B.theta_DR_sd", "B.theta_DR_cp", "B.t_DR.time", "B.r_DR.time",
    "RM.test", "IPTW.test", "DR.test")
result_summary <- apply(result, 2, mean)
print(result_summary)
    

ese <- read.table(paste(outputfilename, "online.ese.csv", sep = ""), sep = ",")[, -1]
ese_col <- apply(ese, 2, sd)

ese_col[p1 + 1]

ese_col[p1 + 1 + p2 +1]

ese_col[p1 + 1 + p2 + 1 + d]





