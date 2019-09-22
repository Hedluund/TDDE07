# Lab 3
################
# Assignment 1 #
################

# Clear plots and data and import data
library(LaplacesDemon)
library(coda)
rm(list = ls())
data <- read.table("rainfall.dat")
data$Rainfall <- data$V1
data$V1 <- c()
n <- nrow(data)
mean_original_data <- mean(data$Rainfall)
sd_original_data <- sd(data$Rainfall)
set.seed(12345)


initValuesGibbs <- c(1,1)

hist(data$Rainfall, breaks=50, freq = FALSE)
dataSamples <- data.frame(matrix(ncol = 3, nrow = 0))

# Set initial values for variables and prior hyperparameters

# Arbitrary start values for the variables
mu <- 1
sigmasq <- 1

# Values that represent our prior belief (how certain we are about the different values)
v0 <- 10
sigmasq0 <- 10
tausq0 <- 1

for (i in 1:1000) {
  
  
  # Calculate parameter values tau n square
  tausqn <- 1/((n/sigmasq) + (1/tausq0))
  
  # Calculate parameter value mu n
  w <- (n/sigmasq) /((n/sigmasq) + (1/tausq0))
  mun <- w*mean_original_data + (1-w)*mu
  
  # Calculate parameter value v n
  vn <- v0 + n
  
  if(i%%2 == 0){
    # Sample mu from the updated parameters
    mu = rnorm(1, mun, sqrt(tausqn))
  }else {
    # Sample sigma from the updated parameters
    sumDiffMuData <- sum((data$Rainfall - mu)^2)
    varianceInvChiSq <- ((v0*sigmasq0 + sumDiffMuData)/ (n + v0))
    sigmasq = rinvchisq(1, df = vn , scale = varianceInvChiSq)
  }
  sample = rnorm(1,mu, sqrt(sigmasq))
  samplerow = c(sample, mu, sigmasq)
  dataSamples = rbind(dataSamples, samplerow)
  
}

datanames <- c("Sampled Rainfall", "mu", "sigmasq")
colnames(dataSamples) <- datanames


# Plotting to see how the values converge and leaving out the first 50 data points.
# This to get a better graphic representation.
par(mfrow = c(1,2))
plot(dataSamples[1:1000,]$mu, dataSamples[1:1000,]$sigmasq, type="o", col="darkgreen",
     main="Mu and sigma without taking burn-in values into account", xlab="Mu", ylab="Sigma squared")
plot(dataSamples[50:1000,]$mu, dataSamples[50:1000,]$sigmasq, type="o", col="darkgreen",
     main="Mu and sigma taking burn-in values into account", xlab="Mu", ylab="Sigma squared")
plot(dataSamples[1:1000,]$`Sampled Rainfall`, type="l", col="darkblue")

############################################## part b ####################################################
# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com

##########    BEGIN USER INPUT #################
# Data options

rawData <- read.table("rainfall.dat")
rawData$Rainfall <- rawData$V1
rawData$V1 = c()

x <- as.matrix(rawData['Rainfall'])

# Model options
# Choosing 2 components according to assignement
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
# Choose 10 as our prior mean
muPrior <- rep(10,nComp) # Prior mean of mu
# Choosing 20 as our prior std
tau2Prior <- rep(20,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 500 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.001 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
#ylim <- c(0,2*max(hist(x)$density))

piVector <- c()
sigma2Vector <- c()
muVector <- c()
mixedDensVector <- c()
for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    # hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      # lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    # lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    # legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
    #       col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  sigma2Vector <- rbind(sigma2Vector, sigma2)
  muVector <- rbind(muVector, mu)
  piVector <- rbind(piVector, pi[1,])
  mixedDensVector <- rbind(mixedDensVector, mixDensMean)
}

par(mfrow = c(1,1))
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

#########################    Helper functions    ##############################################

par(mfrow = c(1,2))
plot(muVector[,1],sqrt(sigma2Vector[,1]), type="o", col="darkblue", main="First component with burnin",
     xlab="Mu", ylab="sigma squared")
plot(muVector[100:500,1],sqrt(sigma2Vector[100:500,1]), type="o", col="darkblue", 
     main="First component without burnin", xlab="Mu", ylab="Sigma squared")

plot(muVector[,2],sqrt(sigma2Vector[,2]), type="o", col="darkred", main="Second component with burnin",
     xlab="Mu", ylab="Sigma squared")
plot(muVector[100:500,2],sqrt(sigma2Vector[100:500,2]), type="o", col="darkred", 
     main="Second component without burnin", xlab="Mu", ylab="Sigma squared")

########################### part c ############################################################
# Calculating posterior mean and std for normal model in assignement a
muhata = mean(dataSamples$mu)
sigma2hata = mean(dataSamples$sigmasq)

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = muhata, sd = sqrt(sigma2hata)) ,type = "l", lwd = 2, col = "blue")
legend("topright", c("Original Data", "Normal model", "Mixed Normal Model"), col= c("black", "darkred", "darkblue"), lty=1, lwd=2)



################
# Assignment 3 #
################


rm(list = ls())
data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
data2 <- data[,-2]

### Part a)
maxLikeHoodModel <- glm(nBids ~. ,data = data2, family = "poisson")

logLik(maxLikeHoodModel)
# The significant features are (> +-0.1): MinBidShare, Sealed, VerifyId, MajBlem and LogBook



### Part b)

postPoisson <- function(betaVect,y,X,mu,Sigma){
  
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  # Calculating log likelihood for
  # logLik <- sum(y*ppois(linPred, log.p = TRUE) + (1-y)*ppois(linPred, log.p = TRUE, lower.tail = FALSE))
  logLik <- sum(y*linPred - exp(linPred))
  
  
  # evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
  
}

# Defining values as input for optim
chooseCov <- c(2:10)

X <- as.matrix(data[,1:10])
nPara <- dim(X)[2];
y <- as.vector(data[,1])
covNames <- names(data)[1:length(names(data))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
mu <- as.vector(rep(0,nPara))
Sigma <- 100*solve(as.matrix(t(X))%*%as.matrix(X))
initVal <- as.vector(rep(0,dim(X)[2])); 

OptimResults<-optim(initVal,postPoisson,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
postMode <- OptimResults$par
postCov <- -solve(OptimResults$hessian) # Posterior covariance matrix is -inv(Hessian)
names(postMode) <- covNames # Naming the coefficient by covariates
# approxPostStd <- sqrt(diag(postCov)) # Computing approximate standard deviations.
# names(approxPostStd) <- covNames # Naming the coefficient by covariates
print('The posterior mode is:')
print(postMode)
# print('The approximate posterior standard deviation is:')
# print(approxPostStd)
print('The posteror covariance matrix (-inv(Hessian)) is:')
print(postCov)

## Part c)
logPostPoisson <- function(theta,y,X,mu,Sigma){
  
  nPara <- length(theta);
  linPred <- X%*%theta;
  
  # Calculating log likelihood for
  # logLik <- sum(y*ppois(linPred, log.p = TRUE) + (1-y)*ppois(linPred, log.p = TRUE, lower.tail = FALSE))
  logLik <- sum(y*linPred - exp(linPred))
  
  
  # evaluating the prior
  logPrior <- dmvnorm(theta, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
  
}

mhRandomWalk <- function(postCov, c, logPostFunc, ...){
  theta <-  as.vector(rep(0,dim(X)[2]))
  thetaVector <- c()
  for(i in 1: 1000){
    newTheta <- as.vector(rmvnorm(1, theta, c*postCov))
    
    logPost <- logPostFunc(theta, ...)
    logPostNew <- logPostFunc(newTheta, ...)
    
    # Calculating the acceptance probability
    alpha <- min(1, exp(logPost - logPostNew))
    alphaCompare <- runif(1, 0, 1)
    if(alphaCompare > alpha){
      theta <- newTheta
    }
    thetaVector <- rbind(thetaVector, theta)
  }
  
  return(thetaVector)
}

c <- 0.1
res <- mhRandomWalk(postCov, c, logPostPoisson, y, X, mu, Sigma)

# Plot how all coefficients converge (or walk)
par(mfrow = c(1,3))
plot(res[,1], type="o", col="blue", cex=0.2, main="How Constant converge", 
     xlab="iteration", ylab="constant value")    
plot(res[,2], type="o", col="blue", cex=0.2, main="How PowerSeller converge", 
     xlab="iteration", ylab="PowerSeller value")  
plot(res[,3], type="o", col="blue", cex=0.2, main="How VerifyID converge", 
     xlab="iteration", ylab="VerifyID value") 
plot(res[,4], type="o", col="blue", cex=0.2, main="How Sealed converge", 
     xlab="iteration", ylab="Sealed value") 
plot(res[,5], type="o", col="blue", cex=0.2, main="How Minblem converge", 
     xlab="iteration", ylab="Minblem value") 
plot(res[,6], type="o", col="blue", cex=0.2, main="How MajBlem converge", 
     xlab="iteration", ylab="MajBlem value") 
plot(res[,7], type="o", col="blue", cex=0.2, main="How LargNeg converge", 
     xlab="iteration", ylab="LargeNeg value") 
plot(res[,8], type="o", col="blue", cex=0.2, main="How LogBook converge", 
     xlab="iteration", ylab="LogBook value") 
plot(res[,9], type="o", col="blue", cex=0.2, main="How MinBidShare converge", 
     xlab="iteration", ylab="MinBidShare value") 