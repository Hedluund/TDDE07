rm(list = ls())
set.seed(12345)
data <- c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
n <- length(data)

mu <- 3.5

taosq <- sum((log(data) - mu)^2)/n 

# in ownage of Mattias Villani
NormalNonInfoPrior<-function(NDraws,Data, tao){

  n<-length(Data)
  PostDraws=matrix(0,NDraws,1)
  PostDraws[,1]<-((n-1)*tao)/rchisq(NDraws,n-1)
 
  return(PostDraws)
}

logNormalDraw <- NormalNonInfoPrior(10000, log(data), taosq)


theoreticalScaledInv <- (((taosq * n/2) ^(n/2))/gamma(n/2)) * exp(-n*taosq/(seq(0, 2, 0.05)*2))/ seq(0, 2, 0.05) ^ (1+n/2) 
par(mfrow = c(1,1))

hist(logNormalDraw, probability = TRUE, breaks = 100 )

lines(seq(0, 2, 0.05) ,theoreticalScaledInv )

## 2 b Culumative normal distribution

cdf = function( x) {
  mu <- 0
  sigma <- 1
  cdfValues = 0.5*(1 + erf((x - mu)/(sigma*sqrt(2))))
  return(cdfValues)
}


# Plotting the CDF values
cdfValues <- 2 * cdf(sqrt(logNormalDraw/2)) -1

# The graph shows that the Gini index is around 0.2 - 0.4 which states relatively equal. The Large coverage in the graph 
# is due top a large variance which is due to a small data set

#Plotting the density of the cdf values
hist(cdfValues, probability = TRUE, breaks = 100, xlim = c(0.1, 0.6))
lines(density(cdfValues))

# Computing the equal tail interval
credibleInterval = p.interval(cdfValues, HPD = FALSE)

# adding the credible interval into the graph
abline(v = credibleInterval, col = "red", lwd = 2, lty = 2)

# Computing the highest posterior density
densityLimit95 <- p.interval(cdfValues)

# When plotting the HDP it is larger (the coverage ) compared to the credible interval
abline(v = densityLimit95, col = "blue", lwd = 2, lty = 2)

#Alternative way of computing the highest posterior density
densityCdfValues = density(cdfValues)
summary(densityCdfValues)

sortedDensityCdf = sort(densityCdfValues$y, decreasing = TRUE)

#Find the limit for 95 %
densityTmp = 0
count = 0
while (densityTmp <= 0.95) {
  count = count + 1
  densityTmp = sum(sortedDensityCdf[1:count])/sum(sortedDensityCdf)
}

densityLimit = sortedDensityCdf[count]

# Find the limit in the original density set
densityLimit95 <- which(densityCdfValues$y >= densityLimit)
# Generate a vector of lower and upper limit

limits <- c(densityCdfValues$x[densityLimit95[1]], densityCdfValues$x[densityLimit95[length(densityLimit95)]])
