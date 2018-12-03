library(MASS)

# mix_weight is k by C matrix of C component weights for each arm
# mu is k by C matrix of C component means for each arm
gaussian_mixture_bandit <- function(alpha,t,mix_weight,mu,cov){
  
  k <- dim(mix_weight)[1]
  C <- dim(mix_weight)[2]
  # initializing n since prior is improper
  
  # Prior parameters
  V_inv0 = 0.01
  gam0=1
  gam = array(rep(1,k*C),dim=c(k,C))
  V_inv = array(rep(V_inv0,k*C),dim=c(k,C))
  u_0 = seq(C)
  u = array(rep(u_0,k),dim=c(k,C))
  alpha_0 = 0.01
  alpha = array(rep(alpha_0,k*C),dim=c(k,C))
  beta_0 = 0.01
  beta = array(rep(beta_0,k*C),dim=c(k,C))
  
  sample_from_comp <- function(a,k){
    w = rnorm(1,u[a,k],sqrt(1/V_inv[a,k]))
    om = rgamma(1,alpha[a,k],beta[a,k])
    return(rnorm(1,w,sqrt(1/om)))
  }
  
  init = NULL
  # observe n rewards from each arm. First column is first arm. Second column is second arm
  components = sapply(seq(C),(function(i){sample.int(C,size=n,replace=TRUE,prob=mix_weight[i,])}))
  for(i in 1:n){
    m = sapply(1:k,(function(j){mu[j,components[j]]}))
    init <- rbind(init,mvrnorm(1,m,cov))
  }
  
  resp <- rep(0,C)
  all_resp <- array(rep(0,t*k*C),dim=c(t,k,C))
  rewards <- matrix(NA,nrow=t,ncol=k)
  rewards[1:n,1:k] <- init
  n = rep(n,k)
  regret <- c(0)
  xbar <- apply(rewards, MARGIN=2,mean,na.rm=TRUE)
  s <- apply(rewards,MARGIN=2,sd,na.rm=TRUE)
  for(i in 1:t){
    
    # sampling from posterior for each of the k arms
    theta.hat <- sapply(1:k,function(a){sum(gam[a,]/sum(gam[a,])*sapply(1:C,function(m)sample_from_comp(a,m)))})
    #theta.hat <- sapply(1:k, function(l) rt(1, df=n[l]+2*alpha-1)
    #                    *sqrt(s[l]/(n[l]*(n[l]+2*alpha-1)))+xbar[l])
    
    # finding best arm and randomly breaking ties
    j <- sample(which(theta.hat ==max(theta.hat)),1)
    # all theoretical rewards
    components = sapply(seq(C),(function(i){sample.int(C,size=1,prob=mix_weight[i,])}))
    m = sapply(seq(k),(function(j){mu[j,components[j]]}))
    obs.rewards <- mvrnorm(1,m,cov)
    
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
    niter = 0
    tot_resp=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]))
    tot_resp_rew=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]*rewards[1:n[j],j]))
    tot_resp_rew2=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]*rewards[1:n[j],j]^2))
    print(tot_resp)
    print(tot_resp_rew)
    # Posterior update
    old_resp=c(10,10)
    print(alpha)
    print(beta)
    print(u)
    print(V_inv)
    print(gam)
    
    while(niter < 50 & max(abs(old_resp-resp))>0.01){
      old_resp = resp
      for(comp in 1:C){
        resp[comp] <- -0.5*(log(beta[j,comp])-digamma(alpha[j,comp])-1/V_inv[j,comp]+(obs.reward-u[j,comp])^2*alpha[j,comp]/beta[j,comp])+digamma(gam[j,comp])-digamma(sum(gam[j,]))
      }
      
      #Normalizing step
      resp = resp-min(resp)
      Z = sum(exp(resp))
      resp = resp-log(Z)
      print("----")
      print(resp)
      
      for(comp in 1:C){
        gam[j,comp]<-gam0+tot_resp[comp]+exp(resp[comp])
        #V_inv_old=V_inv[j,comp]
        V_inv[j,comp]<-V_inv0+tot_resp[comp]+exp(resp[comp])
        #u_old = u[j,comp]
        u[j,comp] <- 1/V_inv[j,comp]*(tot_resp_rew[comp]+exp(resp[comp])*obs.reward+V_inv0*u_0[comp])
        alpha[j,comp]<-alpha_0+0.5*(tot_resp[comp]+exp(resp[comp]))
        beta[j,comp]<-beta_0+0.5*(tot_resp_rew2[comp]+obs.reward^2*exp(resp[comp]))+0.5*u_0[comp]^2*V_inv0-0.5*u[j,comp]^2*V_inv[j,comp]
      }
      niter <- niter + 1
      #print(niter)
    }
    all_resp[n[j],j,]=exp(resp)
  }
  return(rewards)
}

# alpha is parameter defining prior
# X is data. Observations from each arm are in a single row
# i.e. k arms with p_k obs is kxp_k matrix
# Need a measure equivalent to regret since sampling distribution is not specified
# mu is a vector of means of the sampling distribution
# cov is a covariance matrix for the sampling distribution

gaussian_bandit_mix_post <- function(alpha,t,mu,cov){
  k <- length(mu)
  # initializing n since prior is improper
  n <- max( c(2,2-ceiling(2*alpha)))
  C=2
  
  # Prior parameters
  V_inv0 = 0.01
  gam0=1
  gam= array(rep(gam0,k*C),dim=c(k,C))
  V_inv = array(rep(V_inv0,k*C),dim=c(k,C))
  u_0 = seq(C)
  u = array(rep(u_0,k),dim=c(k,C))
  alpha_0 = 0.01
  alpha = array(rep(alpha_0,k*C),dim=c(k,C))
  beta_0 = 0.01
  beta = array(rep(beta_0,k*C),dim=c(k,C))
  
  sample_from_comp <- function(a,k){
    w = rnorm(1,u[a,k],sqrt(1/V_inv[a,k]))
    om = rgamma(1,alpha[a,k],beta[a,k])
    return(rnorm(1,w,sqrt(1/om)))
  }
  
  # initial rewards to get a proper posterior.
  init <-mvrnorm(n,mu, diag(cov))
  all_resp <- array(rep(0,t*k*C),dim=c(t,k,C))
  resp <- rep(0,C)
  rewards <- matrix(NA,nrow=t,ncol=k)
  rewards[1:n,1:k] <- init
  n = rep(n,k)
  regret <- c(0)
  xbar <- apply(rewards, MARGIN=2,mean,na.rm=TRUE)
  s <- apply(rewards,MARGIN=2,sd,na.rm=TRUE)
  for(i in 1:t){
    
    # sampling from posterior for each of the k arms
    theta.hat <- sapply(1:k,function(a){sum(gam[a,]/sum(gam[a,])*sapply(1:C,function(m)sample_from_comp(a,m)))})
    #theta.hat <- sapply(1:k, function(l) rt(1, df=n[l]+2*alpha-1)
    #                    *sqrt(s[l]/(n[l]*(n[l]+2*alpha-1)))+xbar[l])
    
    # finding best arm and randomly breaking ties
    j <- sample(which(theta.hat ==max(theta.hat)),1)
    # all theoretical rewards
    obs.rewards <- mvrnorm(1,mu,diag(cov))
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
    niter = 0
    tot_resp=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]))
    tot_resp_rew=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]*rewards[1:n[j],j]))
    tot_resp_rew2=sapply(seq(C),function(m)sum(all_resp[1:n[j],j,m]*rewards[1:n[j],j]^2))
    print(tot_resp)
    print(tot_resp_rew)
    # Posterior update
    old_resp=c(10,10)
    print(alpha)
    print(beta)
    print(gam)
    #print(u)
    #print(V_inv)
    while(niter < 50 & max(abs(old_resp-resp))>0.01){
      old_resp = resp
      for(comp in 1:C){
        resp[comp] <- -0.5*(log(beta[j,comp])-digamma(alpha[j,comp])-1/V_inv[j,comp]+(obs.reward-u[j,comp])^2*alpha[j,comp]/beta[j,comp])+digamma(gam[j,comp])-digamma(sum(gam[j,]))
      }
      
      #Normalizing step
      resp = resp-min(resp)
      Z = sum(exp(resp))
      resp = resp-log(Z)
      #print("----")
      #print(resp)
      
      for(comp in 1:C){
        gam[j,comp]<-gam0+tot_resp[comp]+exp(resp[comp])
        #V_inv_old=V_inv[j,comp]
        V_inv[j,comp]<-V_inv0+tot_resp[comp]+exp(resp[comp])
        #u_old = u[j,comp]
        u[j,comp] <- 1/V_inv[j,comp]*(tot_resp_rew[comp]+exp(resp[comp])*obs.reward+V_inv0*u_0[comp])
        alpha[j,comp]<-alpha_0+0.5*(tot_resp[comp]+exp(resp[comp]))
        beta[j,comp]<-beta_0+0.5*(tot_resp_rew2[comp]+obs.reward^2*exp(resp[comp]))+0.5*u_0[comp]^2*V_inv0-0.5*u[j,comp]^2*V_inv[j,comp]
      }
      niter <- niter + 1
      #print(niter)
    }
    all_resp[n[j],j,]=exp(resp)
  }
  return(rewards)
}

