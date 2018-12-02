mix_weight = t(array(c(0.5,0.5,0.3,0.7),dim=c(2,2)))
mu = t(array(c(1,2,0,3),dim=c(2,2)))
cov = diag(c(1,1))

reg.I <- lapply(rep(0.5,1000),gaussian_mixture_bandit,t=1000,mix_weight=mix_weight,mu=mu,cov=cov)
