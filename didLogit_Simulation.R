##################################################################################
# SIMULATION STUDY for difference-in-differences estimation with binary response data
# Three Algorithms for computing the standard error: Delta method; Bootstrap; Bayesian 
# Author: Taeyong Park
##################################################################################
library(foreign); library(car); library(MASS); library(nnet); library(lme4); library(arm); library(coda); library(mice); library(xtable); library(texreg); library(apsrtable); library(MCMCpack)


####################
# TRUE VALUES 
####################

## The number of simulation
nSim <- 1000
set.seed(0520)

## TRUE BETA
## A single stimulation study should run only one latentZ
## among the LARGE, MEDIUM, and SMALL cases.
## Once one study is done, move to the next case.

# LARGE latentZ
beta0 <- sample(seq(0, 2, length.out=11), nSim, replace=TRUE)
beta1 <- sample(seq(0, 2, length.out=11), nSim, replace=TRUE)
beta2 <- sample(seq(1, 2, length.out=11), nSim, replace=TRUE)
beta3 <- sample(seq(1, 2, length.out=11), nSim, replace=TRUE)
beta4 <- sample(seq(0, 2, length.out=11), nSim, replace=TRUE)
beta5 <- sample(seq(0, 2, length.out=11), nSim, replace=TRUE)

# MEDIUM latentZ
beta0 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)
beta1 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)
beta2 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)
beta3 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)
beta4 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)
beta5 <- sample(seq(-1, 1, length.out=11), nSim, replace=TRUE)

# SMALL latentZ
beta0 <- sample(seq(-2, 0, length.out=11), nSim, replace=TRUE)
beta1 <- sample(seq(-2, 0, length.out=11), nSim, replace=TRUE)
beta2 <- sample(seq(-2, -1, length.out=11), nSim, replace=TRUE)
beta3 <- sample(seq(-2, -1, length.out=11), nSim, replace=TRUE)
beta4 <- sample(seq(-2, 0, length.out=11), nSim, replace=TRUE)
beta5 <- sample(seq(-2, 0, length.out=11), nSim, replace=TRUE)



# Fixed Xs 
n <- 500 # This is the number of rows for the fake data that we deal with
V1 <- seq(0, 1, length.out=n) + rnorm(n, 0, 0.25)
V2 <- seq(0, 1, length.out=n) + rnorm(n, 0, 0.25)
X <- c(mean(V1),mean(V2))

## TRUE TAU
#X <- matrix(NA, nSim, 2)
tau <- numeric(0)
for (i in 1:nSim) {
  #X[i,] <- sample(seq(0, 10, length.out=11), 2, replace=TRUE) 
  tau[i] <- invlogit(beta0[i]
                     + beta1[i] 
                     + beta2[i]
                     + beta3[i] 
                     + beta4[i]*X[1] 
                     + beta5[i]*X[2]) - invlogit(beta0[i] 
                                                 + beta1[i] 
                                                 + beta2[i] 
                                                 + beta4[i]*X[1] 
                                                 + beta5[i]*X[2])
}

## DATA GENERATING PROCESS ##
set.seed(0520)

# Time
t0 <- rep(0, n/2)
t1 <- rep(1, n/2)
Time <- sample(c(t0,t1), size=n, replace=FALSE)
# Group
control <- rep(0, n/2)
treat <- rep(1, n/2)
Group <- sample(c(control,treat), size=n, replace=FALSE)
# TG
TG <- Time*Group
# V1 and V2
V1 <- V1
V2 <- V2
# Error
e <- rlogis(n, 0, 1)

# Model
latentZ <- probY1 <- Y <- matrix(NA, n, nSim)
for(i in 1:nSim){
  latentZ[,i] <- beta0[i] + beta1[i]*Time + beta2[i]*Group + beta3[i]*TG + beta4[i]*V1 + beta5[i]*V2 + e
  probY1[,i] <- invlogit(latentZ[,i])
  Y[,i] <- rbinom(n, 1, probY1[,i])
}

round(quantile(latentZ, c(0.25, 0.5, 0.75)), 3)



####################
# ESTIMATION 
####################

## 1. APPROXIMATION-BASED (DELTA METHOD - BASED) ANALYTICAL APPROACH ##

# Using latentZ as the response variable
betaHat<-matrix(NA, nSim, 6)
vcovDelta<-array(NA, dim=c(6,6,nSim))
for(i in 1:nSim){
  betaHat[i,]<-coef(lm(latentZ[,i] ~ Time*Group
                       + V1
                       + V2))
  vcovDelta[,,i]<-vcov(lm(latentZ[,i] ~ Time*Group
                          + V1
                          + V2)) # NOTE THAT vcovDelta[6,6,i] is for TIME*GROUP
}

# Using Y as the response variable
#betaHat<-matrix(NA, nSim, 6)
#vcovDelta<-array(NA, dim=c(6,6,nSim))
#for(i in 1:nSim){
#  betaHat[i,]<-coef(glm(Y[,i] ~ Time*Group
#                        + V1
#                        + V2),
#                    family=binomial(link = "logit"))
#  vcovDelta[,,i]<-vcov(glm(Y[,i] ~ Time*Group
#                           + V1
#                           + V2),
#                       family=binomial(link = "logit")) # NOTE THAT vcovDelta[6,6,i] is for TIME*GROUP
#}

beta0Hat<-betaHat[,1]
beta1Hat<-betaHat[,2]
beta2Hat<-betaHat[,3]
beta3Hat<-betaHat[,6]
beta4Hat<-betaHat[,4]
beta5Hat<-betaHat[,5]

## Estimate of DID
tauHatDelta <- numeric(0)
for(i in 1:nSim){
  tauHatDelta[i] <- invlogit(beta0Hat[i]
                             +beta1Hat[i]
                             +beta2Hat[i]
                             +beta3Hat[i]
                             +beta4Hat[i]*X[1]
                             +beta5Hat[i]*X[2])-invlogit(beta0Hat[i]
                                                         +beta1Hat[i]
                                                         +beta2Hat[i]
                                                         +beta4Hat[i]*X[1]
                                                         +beta5Hat[i]*X[2])
}

## Calculate the standard error of tauHat

expon1<-expon2<-partialWithRespectToBeta3<-partialWithRespectToBeta<-numeric(0)
partialWithRespectToTheta1<-partialWithRespectToTheta2<-numeric(0)
vectorLeft<-matrix(NA, nSim, ncol(betaHat))
vectorRight<-matrix(NA, ncol(betaHat), nSim)
tauHatDeltaSE<-confIntDeltaLower<-confIntDeltaUpper<-numeric(0)
coverageDelta<-numeric(0)
for(i in 1:nSim){
  expon1[i]<-exp(beta0Hat[i]
                 +beta1Hat[i]
                 +beta2Hat[i]
                 +beta3Hat[i]
                 +beta4Hat[i]
                 +beta5Hat[i])
  expon2[i]<-exp(beta0Hat[i]
                 +beta1Hat[i]
                 +beta2Hat[i]
                 +beta4Hat[i]
                 +beta5Hat[i])
  partialWithRespectToBeta3[i] <- expon1[i]*((expon1[i]+1)^(-1)+expon1[i]*(-1)*(expon1[i]+1)^(-2))
  partialWithRespectToBeta[i] <- expon1[i]*((expon1[i]+1)^(-1)+expon1[i]*(-1)*(expon1[i]+1)^(-2))-expon2[i]*((expon2[i]+1)^(-1)+expon2[i]*(-1)*(expon2[i]+1)^(-2))
  partialWithRespectToTheta1[i] <- X[1]*expon1[i]*((expon1[i]+1)^(-1)+expon1[i]*(-1)*(expon1[i]+1)^(-2))-X[1]*expon2[i]*((expon2[i]+1)^(-1)+expon2[i]*(-1)*(expon2[i]+1)^(-2))
  partialWithRespectToTheta2[i] <- X[2]*expon1[i]*((expon1[i]+1)^(-1)+expon1[i]*(-1)*(expon1[i]+1)^(-2))-X[2]*expon2[i]*((expon2[i]+1)^(-1)+expon2[i]*(-1)*(expon2[i]+1)^(-2))
  vectorLeft[i,]<-matrix(c(rep(partialWithRespectToBeta[i],3),
                       partialWithRespectToTheta1[i],
                       partialWithRespectToTheta2[i],
                       partialWithRespectToBeta3[i]), 1, ncol(betaHat))
  vectorRight[,i]<-matrix(c(rep(partialWithRespectToBeta[i],3),
                          partialWithRespectToTheta1[i],
                          partialWithRespectToTheta2[i],
                          partialWithRespectToBeta3[i]), ncol(betaHat), 1)
  tauHatDeltaSE[i] <- sqrt(vectorLeft[i,]%*%vcovDelta[,,i]%*%vectorRight[,i])
  confIntDeltaLower[i] <- tauHatDelta[i]-1.96*tauHatDeltaSE[i]
  confIntDeltaUpper[i] <- tauHatDelta[i]+1.96*tauHatDeltaSE[i]
  coverageDelta[i] <- ifelse(tau[i]>confIntDeltaLower[i] & tau[i]<confIntDeltaUpper[i], 1, 0)
}

errorDelta<-sqrt((tauHatDelta-tau)^2)
confIntDelta<-confIntDeltaUpper-confIntDeltaLower
outputDelta <- data.frame(tau, 
                          tauHatDelta,
                          errorDelta,
                          tauHatDeltaSE,
                          confIntDeltaLower,
                          confIntDeltaUpper,
                          confIntDelta,
                          coverageDelta)
round(quantile(outputDelta$tauHatDeltaSE, c(0.25, 0.5, 0.75)), 3)

## Save the output
# outputDeltaLarge<-outputDelta
# outputDeltaSmall<-outputDelta
# outputDeltaMedium<-outputDelta



## 2. BOOTSTRAPPING ##

## Consider one set of data that produces one set of betaHat
## This one set of data consists of n=500 rows and 5 columns (latentZ, Time, Group, V1, V2) 
## To compute the bootstrap estimates and standard errors, we want to simulate this one set of data B times with allowing replacement.
## Then, we will end up with B sets of betaHat instead of one set of betaHat --> The mean of B sets of betaHat is the bootstrap mean, 
## and the s.d. of B sets of betaHat gives you the bootstrap SE

## NOTE that the above process is essentially the same as the bootstrpping done by didLogit.R
## Here, we should have one more step because we are now dealing with nSim sets of parameters. 

## We apply this process to nSim sets. 
## Each of nSim sets produces B sets of betaHat --> one bootstrap mean and one bootstrap SE  
## As a result, nSim sets produce nSim sets of bootstrap mean and nSim sets of bootstrap SE.
## This final output will be corresponding to the final output by the Delta method, outputDelta


B <- 1000 # Number of bootstrap simulations
eachNsimData <- array(NA, dim=c(n, 5, nSim))
bootData <- array(NA, dim=c(n, 5, B)) # n=500; 5 is the number of the columns
coefBoot<-list()
for(j in 1:nSim){
  eachNsimData[,,j] <- cbind(latentZ[,j], Time, Group, V1, V2)
  for (i in 1:B){
    sampleObs<-sample(1:n, n, replace=TRUE)
    bootData[,,i]<-eachNsimData[,,j][sampleObs, ]
  }
  colnames(bootData) <- c("latentZ", "Time", "Group", "V1", "V2")
  
  linearBootstrap <- function (i) {
    bootModel<-lm(latentZ~Time*Group
                  +V1
                  +V2,
                  data=data.frame(bootData[,,i]))
    return(list(coef(bootModel)))
  }
  coefBoot[[j]] <- sapply(1:B, FUN=linearBootstrap)
}

tauHatBoot<-tauHatBootMean<-tauHatBootSE<-numeric(0)
confIntBootLower<-confIntBootUpper<-coverageBoot<-numeric(0)
for(j in 1:nSim){
  for(i in 1:B){
    tauHatBoot[i]<-invlogit(coefBoot[[j]][[i]][1]
                            +coefBoot[[j]][[i]][2]
                            +coefBoot[[j]][[i]][3]
                            +coefBoot[[j]][[i]][6]
                            +coefBoot[[j]][[i]][4]*X[1]
                            +coefBoot[[j]][[i]][5]*X[2])-invlogit(coefBoot[[j]][[i]][1]
                                                                  +coefBoot[[j]][[i]][2]
                                                                  +coefBoot[[j]][[i]][3]
                                                                  +coefBoot[[j]][[i]][4]*X[1]
                                                                  +coefBoot[[j]][[i]][5]*X[2])
  }
  tauHatBootMean[j]<-mean(tauHatBoot)
  tauHatBootSE[j]<-sd(tauHatBoot)
  confIntBootLower[j] <- tauHatBootMean[j]-1.96*tauHatBootSE[j]
  confIntBootUpper[j] <- tauHatBootMean[j]+1.96*tauHatBootSE[j]
  coverageBoot[j] <- ifelse(tau[j]>confIntBootLower[j] & tau[j]<confIntBootUpper[j], 1, 0)
}

errorBoot<-sqrt((tauHatBootMean-tau)^2)
confIntBoot<-confIntBootUpper-confIntBootLower
outputBoot <- data.frame(tau, 
                         tauHatBootMean,
                         errorBoot,
                         tauHatBootSE,
                         confIntBootLower,
                         confIntBootUpper,
                         confIntBoot,
                         coverageBoot)

round(quantile(outputBoot$tauHatBootSE, c(0.25, 0.5, 0.75)), 3) 

## Save the output
# outputBootLarge <- outputBoot
#outputBootSmall <- outputBoot
# outputBootMedium <- outputBoot

## 3. BAYESIAN ## 
nIter<-3000
posterior <- array(NA, dim=c(nIter, ncol(betaHat)+1, nSim))
beta0HatPosterior<-matrix(NA, nSim, nIter)
beta1HatPosterior<-matrix(NA, nSim, nIter)
beta2HatPosterior<-matrix(NA, nSim, nIter)
beta3HatPosterior<-matrix(NA, nSim, nIter)
beta4HatPosterior<-matrix(NA, nSim, nIter)
beta5HatPosterior<-matrix(NA, nSim, nIter)
tauHatPosterior<-matrix(NA, nSim, nIter)
tauHatPosteriorMean<-numeric(0)
tauHatPosteriorSE<-numeric(0)
lowerPosterior90<-lowerPosterior95<-lowerPosterior99<-numeric(0)
upperPosterior90<-upperPosterior95<-upperPosterior99<-numeric(0)
coveragePosterior<-numeric(0)
for(i in 1:nSim){
  posterior[,,i] <- mcmc(MCMCregress(latentZ[,i]~Time*Group
                                   + V1
                                   + V2, b0=0, B0=0, mcmc=nIter))
  beta0HatPosterior[i,] <- posterior[,,i][,1]
  beta1HatPosterior[i,] <- posterior[,,i][,2]
  beta2HatPosterior[i,] <- posterior[,,i][,3]
  beta3HatPosterior[i,] <- posterior[,,i][,6]
  beta4HatPosterior[i,] <- posterior[,,i][,4]
  beta5HatPosterior[i,] <- posterior[,,i][,5]
  tauHatPosterior[i,]<- invlogit(beta0HatPosterior[i,]
                                 +beta1HatPosterior[i,]
                                 +beta2HatPosterior[i,]
                                 +beta3HatPosterior[i,]
                                 +beta4HatPosterior[i,]*X[1]
                                 +beta5HatPosterior[i,]*X[2])-invlogit(beta0HatPosterior[i,]
                                                                       +beta1HatPosterior[i,]
                                                                       +beta2HatPosterior[i,]
                                                                       +beta4HatPosterior[i,]*X[1]
                                                                       +beta5HatPosterior[i,]*X[2])
  tauHatPosteriorMean[i]<-mean(tauHatPosterior[i,])
  tauHatPosteriorSE[i]<-sd(tauHatPosterior[i,])
  lowerPosterior90[i]<-quantile(tauHatPosterior[i,], 0.05)
  upperPosterior90[i]<-quantile(tauHatPosterior[i,], 0.95)
  lowerPosterior95[i]<-quantile(tauHatPosterior[i,], 0.025)
  upperPosterior95[i]<-quantile(tauHatPosterior[i,], 0.975)
  lowerPosterior99[i]<-quantile(tauHatPosterior[i,], 0.005)
  upperPosterior99[i]<-quantile(tauHatPosterior[i,], 0.995)
}

errorPosterior<-sqrt((tauHatPosteriorMean-tau)^2)
coveragePosterior90 <- ifelse(tau>lowerPosterior90 & tau<upperPosterior90, 1, 0)
coveragePosterior95 <- ifelse(tau>lowerPosterior95 & tau<upperPosterior95, 1, 0)
coveragePosterior99 <- ifelse(tau>lowerPosterior99 & tau<upperPosterior99, 1, 0)
credInt90<-upperPosterior90-lowerPosterior90
credInt95<-upperPosterior95-lowerPosterior95
credInt99<-upperPosterior99-lowerPosterior99
outputPosterior <- data.frame(tau, 
                              tauHatPosteriorMean,
                              errorPosterior,
                              tauHatPosteriorSE,
                              lowerPosterior90,
                              upperPosterior90,
                              credInt90,
                              lowerPosterior95,
                              upperPosterior95,
                              credInt95,
                              lowerPosterior99,
                              upperPosterior99,
                              credInt99,
                              coveragePosterior90,
                              coveragePosterior95,
                              coveragePosterior99)
round(quantile(outputPosterior$tauHatPosteriorSE, c(0.25, 0.5, 0.75)), 3)

## Save the output
# outputPosteriorLarge<-outputPosterior
# outputPosteriorSmall<-outputPosterior
# outputPosteriorMedium<-outputPosterior

## COMPARISONS OF THE TRHEE APPROACHES ##

comparison<-cbind(outputDelta$coverageDelta, 
                  outputBoot$coverageBoot,
                  outputPosterior$coveragePosterior95)
compareSE<-cbind(outputDelta$tauHatDeltaSE, 
                 outputBoot$tauHatBootSE,
                 outputPosterior$tauHatPosteriorSE)
deltaIsMostConservative<-numeric(0)
bootIsMoreConservative<-numeric(0)
for(i in 1:nSim){
  deltaIsMostConservative[i]<-ifelse(outputDelta$tauHatDeltaSE[i]<outputBoot$tauHatBootSE[i] & outputDelta$tauHatDeltaSE[i]<outputPosterior$tauHatPosteriorSE[i],
                                     1, 0)
  bootIsMoreConservative[i]<-ifelse(outputBoot$tauHatBootSE[i]<outputPosterior$tauHatPosteriorSE[i],
                                    1, 0)
}









