mix_weight = t(array(c(0.5,0.5,0.3,0.7),dim=c(2,2)))
mu = t(array(c(1,2,0,3),dim=c(2,2)))
cov = diag(c(1,1))


reg.K2 <- lapply(rep(0.5,1000),gaussian_mixture_bandit,t=1000,mix_weight=mix_weight,mu=mu,cov=cov)
mu=c(1,1,1.1,1,1)
cov=rep(1,5)
reg.K1 <- lapply(rep(0.5,1000),gaussian_bandit_mix_post,t=1000,mu=mu,cov=cov)

mix_weight2 = t(array(c(0.5,0.5,0.5,0.5),dim=c(2,2)))
mu = t(array(c(2,3,1,0),dim=c(2,2)))
reg.dif <- lapply(rep(0.5,1000),gaussian_mixture_bandit,t=1000,mix_weight=mix_weight,mu=mu,cov=cov)

library(Matrix)
plot(Reduce('+', reg.I)/1000, type="l")
lines(Reduce('+', reg.dif)/1000)
