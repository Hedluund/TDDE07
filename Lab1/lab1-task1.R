rm(list = ls())
set.seed(12345)
alpha0 <- 2
beta0 <- 2 

n1 <- 20
s <- 14
f <- n1-s

randomDrawsBeta <- rbeta(20, alpha0 + s, beta0 + f)

alphaPost <- alpha0 + s
betaPost <- beta0 + f

drawMean <- mean(randomDrawsBeta)

drawSd <- sd(randomDrawsBeta)

expectedMean <- alphaPost /(alphaPost + betaPost)

expectedSd <- ((alphaPost * betaPost) / ((alphaPost + betaPost)^2 * (alphaPost + betaPost + 1))) ^0.5

drawMeanVector <- c()
drawSdVector <- c()

nDraws <- seq(1, 10000, 20)

for(i in nDraws){
  
  randomDrawsBeta <- rbeta(i, alpha0 + s, beta0 + f)
  newMean <- mean(randomDrawsBeta)
  newSd <- sd(randomDrawsBeta)
  drawMeanVector <- c(drawMeanVector, newMean)
  drawSdVector <- c(drawSdVector, newSd)
}
par(mfrow = c(1,2))
plot(nDraws, drawMeanVector, ylim = c(0.64, 0.7), pch = 16, col = "darkblue", main="Mean values", 
     xlab = "Number of draws", ylab = "Mean")
abline(h = expectedMean, col = "red")
plot(nDraws, drawSdVector , pch = 16, col = "darkgreen", main="Standard deviation values", 
     xlab = "Number of draws", ylab = "Standard deviation")
abline(h = expectedSd, col = "red")

# Emma can se that the points very converge

# --------------------------------------- 
# 1 b)


exactPostProbs0.4 <- pbeta(c(0.4), alphaPost, betaPost)

randomDrawsBeta10000 <- rbeta(10000, alpha0 + s, beta0 + f)

amountBelow0.4 <- length(which(randomDrawsBeta10000 < 0.4))

percentage0.4 <- amountBelow0.4/10000

# 1 c)

logodds1000 <- log(randomDrawsBeta10000/(1 - randomDrawsBeta10000))

hist(logodds1000 )
logoddsdens <- density(logodds1000)


