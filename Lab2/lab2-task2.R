rm(list = ls())
library(mvtnorm)
library(LaplacesDemon)


data = read.table("WomenWork.dat", header = TRUE)

# The "0" in the model is so that the model will not add an intercepts since the data already has one
glmModel <- glm(Work~0 + . , data = data , family = binomial)

# 1b)
# Setting y and X to vector and matrix to be able to input into optim function
y <- as.vector(data[,1])
X <- as.matrix(data[,2:9])

# Setting the names of the covariates to be able to analyse the output
covNames <- names(data)[2:9]
nParam <- dim(X)[2]
mu <- rep(0,nParam)
tao <- 10
sigma <- tao^2*diag(nParam)

# the log posterior function
logPostLogistic <- function(betaVector,y,X,mu,sigma){
  
  # Computes the number of parameters
  nPara <- length(betaVector);
  
  # Predicts the y value by multiplying X and beta
  linPrediction <- X%*%betaVector;
  
  # Calculating the log-likelihood
  # Ask about this formula!
  logLikelihood <- sum( linPrediction*y -log(1 + exp(linPrediction)));
  
  # Likelihood is not finite, stear the optimizer away from here!
  if (abs(logLikelihood) == Inf) logLikelihood = -20000; 
  
  # Draw the prior beta (log)
  logPrior <- dmvnorm(betaVector, matrix(0,nPara,1), sigma, log=TRUE);
  
  # Add the log prior and log-likelihood together to get log posterior
  return(logLikelihood + logPrior)
}

# Set randomized initial values for beta
initBetas <- as.vector(rnorm(dim(X)[2]))

# Get the posterior mode and the Hessian by using the optim function
optimRes <- optim(initBetas, logPostLogistic, gr=NULL, y, X, mu, sigma, method=c("BFGS"), control=list(fnscale=-1), hessian=TRUE)

# Get the posterior mode and the posterior covariance matrix (-inverse Hessian) from the optimal results
postMode <- optimRes$par
postCovariance <- -solve(optimRes$hessian) 
names(postMode) <- covNames
approxPostStd <- sqrt(diag(postCovariance))
names(approxPostStd) <- covNames

# Draw simulations of NSmallChild
nDraws <- 1000

NSmallChildDraws <- rnorm(nDraws, postMode["NSmallChild"], approxPostStd["NSmallChild"])

# Compute the approximate credible interval of NSmallChild
lowerLimit <- quantile(NSmallChildDraws, 0.025)
upperLimit <- quantile(NSmallChildDraws, 0.975)

# Plotting the distribution
hist(NSmallChildDraws, probability = TRUE)
lines(density(NSmallChildDraws))
abline(v = c(lowerLimit, upperLimit), col = "red")

#Since the intervall does not include zero and has an expected value of -1.359 we can say that it is important


# 1c)

predDist <- function(xValues,postMod, postCov) {
  simulatedBetas <- as.vector(rmvnorm(1, postMod,postCov))
  
  #Calculate approximate probability
  expValue <-  exp(xValues %*%  simulatedBetas)
  prob <- expValue/(1 + expValue) 
  # Draws from bernoulli from our apporximated probability 
  bernDraw <- rbern(1,prob)
  return(bernDraw)
}

predDisitributionVector <- c()
womanFeatures <- c(1, 10, 8, 10, 1, 40, 1, 1)

for(i in 1:nDraws){
  pred <- predDist(womanFeatures, postMode, postCovariance)
  predDisitributionVector <- c(predDisitributionVector, pred)
}

hist(predDisitributionVector, probability = TRUE, breaks = 2)
legend("hej hej")
