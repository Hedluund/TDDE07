rm(list = ls())
library(mvtnorm)
library(matlib)

data <- read.table("TempLinkoping.txt", header = TRUE)
set.seed(12345)
omega0 <- 0.01*diag(3)
v0 <- 4
sigmaSquare0 <- 1
mu0 <- matrix(c(-10,100,-100))
timeVector = c(1:366)

# Function that draws sigramsquare from inverse chi 2 distribution
scaledInvChiSqGenerator<-function(NDraws, sigmaSquare, df){
  
  PostDraws=matrix(0,NDraws,1)
  PostDraws[,1]<-((df)*sigmaSquare)/rchisq(NDraws, df)
  
  return(PostDraws)
}
# Simulate 10 sigma draws over a sequens of 366 days

ndraws = 10

calcTemp = function(beta, sigmaSquare){
  temp = c()
  for (t in data$time) {
    error <- rnorm(1 , 0, sqrt(sigmaSquare))
    newtemp <-beta[1]+beta[2]*t+beta[3]*t^2 + error
    
    temp=c(temp,newtemp)
  }
  return(temp)
}


# Plotting original values
par(mfrow = c(1,2))
tempVector = c()
tempinfoVector = c()
plot(tempVector, ylim=c(-20,40), xlim = c(0,370), main = "Prior from original parameter values")

for (i in 1:ndraws) {
  sigmaDraw <- scaledInvChiSqGenerator(1, sigmaSquare0, v0 )
  betadraw <- rmvnorm(1, mean = mu0, sigma = sigmaDraw[1]*solve(omega0))
  temps <- calcTemp(betadraw[1,], sigmaSquare0)
  tempVector = rbind(tempVector,temps)
  lines(temps)
}

# Plotting modified values
# We see that the vaiance is large, therefore we want to increase omega0. Also we change the mu0 to get a little curvier curve. ie some more range
# in lowest and highest temperatures.
omega0 <- 1*diag(3)
mu0 <- matrix(c(-10,130,-130))
tempVector = c()
tempinfoVector = c()
plot(tempVector, ylim=c(-20,40), xlim = c(0,370), main = "Prior from modified parameter values")

  
for (i in 1:ndraws) {
  sigmaDraw <- scaledInvChiSqGenerator(1, sigmaSquare0, v0 )
  betadraw <- rmvnorm(1, mean = mu0, sigma = sigmaDraw[1]*solve(omega0))
  temps <- calcTemp(betadraw[1,], sigmaSquare0)
  tempVector = rbind(tempVector,temps)
  lines(temps)
}


#seems kinda resonable
# 1b) joint posterior distribution

# Function to caculate X
# Solve() is the inverse of a matrix
X = c()
for (time in data$time) {
  xNew = c(1, time, time^2)
  X = rbind(X, xNew)
}

# Computing all parameters
X_X = t(X) %*% X
betaHat <- solve(X_X) %*% t(X) %*% data$temp
mun <- solve(X_X + omega0) %*% (X_X %*% betaHat + omega0 %*% mu0)
omegan <- X_X + omega0
vn <- v0 + length(data[,2])
sigmaSquaren <- (v0*sigmaSquare0 + (t(data$temp) %*% data$temp + t(mu0) %*% omega0 %*% mu0 - t(mun) %*% omegan %*% mun))/vn

# Draw beta and sigmasquare from joint posterior distribution and plot
tempVector = c()
sigmaVector = c()
betaVector = c()
par(mfrow = c(1,1))
plot(c(), ylim=c(-20,40), col = "blue", xlim = c(0,370), main="Draws from the posterior")
legend("topleft", c("temp values", "posterior draws"), col = c("blue", "black"), pch = 16)

ndraws <- 100
for (i in 1:ndraws) {
  sigmaDraw <- scaledInvChiSqGenerator(1, sigmaSquaren, vn)
  sigmaVector = c(sigmaVector, sigmaDraw)
  
  betadraw <- rmvnorm(1, mean = mun, sigma = sigmaDraw[1]*solve(omegan))
  betaVector = rbind(betaVector, betadraw)
  
  temps <- calcTemp(betadraw[1,], sigmaDraw)
  tempVector <- rbind(tempVector,temps)
  #lines(temps)
}

# Calculating the median betas to be able to calculate the median posterior
medianBeta0 <- median(betaVector[,1])
medianBeta1 <- median(betaVector[,2])
medianBeta2 <- median(betaVector[,3])

# Computing the median posterior of the regression function
medianPosteriorTemp <- medianBeta0 + medianBeta1*data$time + medianBeta2*data$time^2
lines(medianPosteriorTemp)

# Plotting the original values
points(data$temp, col ="blue")

# Compute posterior credible interval
lowerLimitVector <- c()
upperLimitVector <- c()
for (i in 1:ncol(tempVector)) {
  lowerLimit <- quantile(tempVector[,i], probs = c(0.025))
  lowerLimitVector <- c(lowerLimitVector, lowerLimit[[1]])
  
  upperLimit <- quantile(tempVector[,i], probs = c(0.975))
  upperLimitVector <- c(upperLimitVector, upperLimit[[1]])
}

lines(lowerLimitVector, col = "red")
lines(upperLimitVector, col = "red")

# Plot histogram of sigma and beta
hist(sigmaVector, main = "Marginal posterior for sigma", breaks = 20)
hist(betaVector[,1], main = "Marginal posterior for Beta 0", breaks = 20)
hist(betaVector[,2], main = "Marginal posterior for Beta 1", breaks = 20)
hist(betaVector[,3], main = "Marginal posterior for Beta 2", breaks = 10)


# 1c) Posterior distribution for x~ using median beta
xtildeVector = c()
for (row in 1:nrow(betaVector)) {
  # using derivate of bet*x = 0
  betas <- betaVector[row,]
  xtilde <- (-betas[2]/(2*betas[3])*366)
  xtildeVector <- c(xtildeVector,xtilde)
}

hist(xtildeVector,breaks = 10, freq = FALSE)

# we can see that the xtilde is normally distributed with an expected value of 197

# 1d)
# For omeg07 we set the values for the high beta indexes lower and lower because of the uncertainty.
# For mu07 we do the same and set for the high beta indexes lower and lower because of the uncertainty
# By lower and lower we say that they should converge to zero
