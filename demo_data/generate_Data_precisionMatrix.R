##############################################################################################
# Two tasks:
# Task 1. Generate a 50 X 50 precision matrix "sigma_inv" whose off-diagonal elements are nonzero with probability 5%;
# Task 2. Simulate 100 samples of data from a 50-dimensional multivariate normal distribution with mean 0 and precision matrix  
# "sigma_inv" generated from Task 1.

##############################################################################################

# Task 1. Generate precision matrix p=50; (5% nonzero) (random structure)

##############################################################################################
rm(list=ls())
if (!require(mvtnorm)) install.packages('mvtnorm')
if (!require(moments)) install.packages('moments')
if (!require(Matrix)) install.packages('Matrix')
library(mvtnorm)
library(moments)
library(Matrix)
set.seed(2021)
path_0 <- './demo_data/'
p <- 50
n <- 100 
# Off-diagonal elements are nonzero with probability 5%
prob_nozero <- 0.05
# Generate off-diagonal elements in "sigma_inv"
a = rep(NA,p*(p-1)/2)
for (i in 1:(p*(p-1)/2)) {
  if (runif(1)>prob_nozero) {a[i] = 0}    
  else{
    a[i] = runif(1, 3, 6)
  }
}
A <- matrix(NA,nrow=p,ncol=p)
A[upper.tri(A)] <- a
diag(A) <- rgamma(p,shape=2,rate=0.05) # for p =50
sigma_inv <- forceSymmetric(A)
sigma_inv <- as.matrix(sigma_inv)
# if "sigma_inv" is not positive-definite, generate again till it is
while (min(eigen(sigma_inv)$values) < 0.01) {
  for (i in 1:(p*(p-1)/2)) {
    if (runif(1)>prob_nozero) {a[i] = 0}
    else{
      a[i] = runif(1, 3, 6)
    }
  }
  A <- matrix(NA,nrow=p,ncol=p)
  A[upper.tri(A)] <- a
  diag(A) <- rgamma(p,shape=2,rate=0.05) # for p =50
  sigma_inv <- forceSymmetric(A)
  sigma_inv <- as.matrix(sigma_inv)
}
# eigen(sigma_inv)$values
sigma <- solve(sigma_inv)
sigma <- as.matrix(sigma)
write.table(sigma_inv,file=paste0(path_0,"GNBP_sim_p50_sigmainv5(3-6).csv"),sep=",",row.names=FALSE,col.names=FALSE)

##############################################################################################

# Task 2. Simulate m data sets with each data set containing 100 samples

##############################################################################################

m <- 1 # generate one data set for demo
omega_elements <- t(sigma_inv)[lower.tri(sigma_inv,diag=FALSE)]
sum(omega_elements!=0)
xx_list <- list(NA)
for (i in 1:m) {
  set.seed(2050+i)
  xx_list[[i]] <- rmvnorm(n=n,mean=rep(0,p),sigma=sigma)
  write.table(xx_list[[i]],file=paste0(path_0,"GNBP_sim_p50n100_5nonzero_data(3-6)", i, ".csv"),sep=",",row.names=FALSE,col.names=FALSE)
}

