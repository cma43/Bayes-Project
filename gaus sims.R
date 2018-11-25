mu <- rep(1,5); mu[3] <- mu[3] + 0.1
cov <- diag(rep(1,5))
mu2 <- rep(1,5); mu[3] <- mu[3] + 0.01

reg.I <- lapply( rep(0.5,1000),gaussian_bandit,t=1000,mu=mu,cov=cov)

library(Matrix)
cov2 = matrix(rep(1,5^2),nrow=5)
epsilon = .90
for(i in 1:5){
  for(j in 1:5){
    cov2[i,j] <- epsilon^(j-i+1)  
  }
  
}
cov2 <- forceSymmetric(cov2, "U")

reg.quad <- lapply( rep(0.5,1000),gaussian_bandit,t=1000,mu=mu,cov=cov2)

plot(Reduce('+', reg.I)/1000, type="l")
lines(Reduce('+', reg.quad)/1000)

