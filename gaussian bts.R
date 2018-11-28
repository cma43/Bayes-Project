# alpha is parameter defining prior
# X is data. Observations from each arm are in a single row
# i.e. k arms with p_k obs is kxp_k matrix
# Need a measure equivalent to regret since sampling distribution is not specified
# mu is a vector of means of the sampling distribution
# cov is a covariance matrix for the sampling distribution

library(MASS)
gaussian_bandit2 <- function(alpha,t,mu,cov){
  k <- length(mu)  
  # initializing n since prior is improper
  n <- max( c(2,2-ceiling(2*alpha)))
  # initial rewards to get a proper posterior.
  
  # observe n rewards from each arm. First column is first arm. Second column is second arm
  init <-mvrnorm(n,mu, diag(cov))
  rewards <- matrix(NA,nrow=t,ncol=k)
  rewards[1:n,1:k] <- init
  n = rep(n,k)
  regret <- c(0)
  xbar <- apply(rewards, MARGIN=2,mean,na.rm=TRUE)
  s <- apply(rewards,MARGIN=2,sd,na.rm=TRUE)
  add = 0
  j.add = 0
  r = c(0)
  for(i in 1:t){
    
    # sampling from posterior for each of the k arms
    theta.hat <- sapply(1:k, function(l) rt(1, df=n[l]+2*alpha-1)
                        *sqrt(s[l]/(n[l]*(n[l]+2*alpha-1)))+xbar[l] )
    
  # finding best arm and randomly breaking ties
    j <- which.max(theta.hat)
    # all theoretical rewards
    obs.rewards <- mvrnorm(1,mu,diag(cov))
    obs.rewards[j] <- obs.rewards[j]+ add*(j==j.add)
    obs.rewards[j.add] <- obs.rewards[j.add] + add*(j !=j.add)
    # Use observed maximum instead of expected maximum to avoid computation error
    max.reward <- max(obs.rewards)
    # actual observed reward for best arm at time t
    obs.reward <- obs.rewards[j]
    #update n
    n[j] = n[j]+1
    
    # updating posterior stuff
    rewards[n[j],j] <- obs.reward
    xbar[j] <- mean(rewards[,j], na.rm=TRUE)
    s[j] <- sd(rewards[,j],na.rm=TRUE)
    regret[i+1] <- regret[i] + max.reward-obs.reward
    add <- 2*rbeta(1,1/2,1/2)
    j.add = j+1
    if(j==length(mu)){
      j.add = 1
    }
    r[i] = obs.reward
  }
  return(rewards)
  #return(sum(r))
}



