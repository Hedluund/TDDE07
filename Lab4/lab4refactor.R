rm(list = ls())
library(rstan)
library(coda)

data = read.table("campy.dat")
set.seed(12345)


######################################### Task 1a ########################################################
T <- 200
sigmasq <- 2
mu <- 10 

x_pre <- mu

phiVector <- seq(-1, 1, 0.25)
phiLength <- length(phiVector)
x_matrix <- c()
par(mfrow = c(3,3))
for (k in 1:phiLength) {
  x_vector <- c()
  
  phi <- phiVector[k]
  for (i in 1:T) {
    error <- rnorm(1, 0, sqrt(sigmasq))
    
    x <- mu + phi*(x_pre - mu) + error
    
    x_pre <- x
    x_vector <- rbind(x_vector, x)
  }
  
  plot(x_vector, type = "l", col = "darkgreen", main = paste0("phi : ", + phi))
  
  x_matrix <- cbind(x_matrix, x_vector)
}

# phi affects of dependent the current simulation data is on the last simulated data

######################################## Task 1b #########################################################

# Sampling data for Stan
par(mfrow = c(1,2))
phiVector <- c(0.3, 0.95)
phiLength <- length(phiVector)
mu <- 10
sigmasqb <- 2
x_pre <- mu
x_init <- mu
x_matrix <- c()

for (k in 1:phiLength) {
  x_vector <- c()
  phi <- phiVector[k]
  for (i in 1:T) {
    error <- rnorm(1, 0, sqrt(sigmasqb))
    
    x <- mu + phi*(x_pre - mu) + error
    x_pre <- x
    x_vector <- rbind(x_vector, x)
  }
  
  plot(x_vector, type = "l", col = "darkblue", main = paste0("phi : ", + phi))
  
  x_matrix <- cbind(x_matrix, x_vector)
}

par(mfrow = c(1,1))

# Stan model data
x <- as.vector(x_matrix[,1])
y <- as.vector(x_matrix[,2])
N <- length(x)

# Stan model for phi = 0.3
StanModelx = '
data {
  int<lower=0> N; //Number of sampled observations
  real<lower=0> x[N]; // Value of the observations
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real phi;
}
model {

  for(i in 2:N)
  x[i] ~ normal(mu + phi*(x[i-1] - mu), sqrt(sigma2));
}
'


dataPhi30 = list(N=N, x=x)
burnin = 300
niter = 1000
fit30 = stan(model_code=StanModelx,data=dataPhi30,
             warmup=burnin,iter=niter,chains=4)
# Print the fitted model
print(fit30,digits_summary=2)
# Extract posterior samples
postDraws30 <- extract(fit30)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws30$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains


# Stan model for phi = 0.95

dataPhi95 = list(N=N, x=y)

fit95 = stan(model_code=StanModelx,data=dataPhi95,
             warmup=burnin,iter=niter,chains=4, control = list(max_treedepth=18))
# Print the fitted model
print(fit95,digits_summary=2)
# Extract posterior samples
postDraws95 <- extract(fit95)
# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws95$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains



# Calculating posterior mean credible intervals and effective posterior samples for
# mu sigma and phi

postMeanMu30 <- mean(postDraws30$mu)
postMeanMu95 <- mean(postDraws95$mu)
postSigma30 <- mean(postDraws30$sigma2)
postSigma95 <- mean(postDraws95$sigma2)
postPhi30 <- mean(postDraws30$phi)
postPhi95 <- mean(postDraws95$phi)


credIntervalMu30 <- quantile(postDraws30$mu, c(0.025, 0.975))
credIntervalMu95 <- quantile(postDraws95$mu, c(0.025, 0.975))
credIntervalSigma30 <- quantile(postDraws30$sigma2, c(0.025, 0.975))
credIntervalSigma95 <- quantile(postDraws95$sigma2, c(0.025, 0.975))
credIntervalPhi30 <- quantile(postDraws30$phi, c(0.025, 0.975))
credIntervalPhi95 <- quantile(postDraws95$phi, c(0.025, 0.975))

par(mfrow = c(1,2))
################# plotting for mu
plot(density(postDraws30$mu), main="mean mu for phi =0.3")
abline(v=c(postMeanMu30), col = "darkblue")
abline(v= c(credIntervalMu30), col = "red")

plot(density(postDraws95$mu), main="mean mu for phi =0.95")
abline(v=c(postMeanMu95), col = "darkblue")
abline(v= c(credIntervalMu95), col = "red")

################## plotting for phi
plot(density(postDraws30$phi), main="mean phi for phi =0.3")
abline(v=c(postPhi30), col = "darkblue")
abline(v= c(credIntervalPhi30), col = "red")

plot(density(postDraws95$phi), main="mean phi for phi =0.95")
abline(v=c(postPhi95), col = "darkblue")
abline(v= c(credIntervalPhi95), col = "red")

################## plotting for sigma
plot(density(postDraws30$sigma2), main="mean sigma for phi =0.3")
abline(v=c(postSigma30), col = "darkblue")
abline(v= c(credIntervalSigma30), col = "red")

plot(density(postDraws95$sigma2), main="mean sigma for phi =0.95")
abline(v=c(postSigma95), col = "darkblue")
abline(v= c(credIntervalSigma95), col = "red")



#ess <- ((niter - burnin)/( 1+ 2*)
# ess taken from the print 
# essmu30 = 2094 esssigma30= 3033 essphi30 = 2086
# essmu95 = 1473 esssigma95= 2114 essphi95 = 1493


################################# FRAGA: Hur data vi fran joint posterior och ar det resultatet vi far rimligt?
###### Answer: Plot phi vs mu in same plot
#plotting posterior

plot(postDraws30$mu, postDraws30$phi)
plot(postDraws95$mu, postDraws95$phi)

# comment: phi seems to be centerd around 0.4 abd mu around 10 in the phi30 plot,
# comment: phi seems to be centerd around  0.9 and mu around 9 in the phi 90 plot