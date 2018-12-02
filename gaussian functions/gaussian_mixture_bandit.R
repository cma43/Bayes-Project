library(MASS)



# mix_weight is k by C matrix of C component weights for each arm
# mu is k by C matrix of C component means for each arm
gaussian_mixture_bandit <- function(alpha,t,mix_weight,mu,cov){
  
  k <- dim(mix_weight)[1]
  C <- dim(mix_weight)[2]
  # initializing n since prior is improper
  
  # Prior parameters
  V_inv0 = 0.01
  gam = array(rep(1,k*C),dim=c(k,C))
  V_inv = array(rep(V_inv0,k*C),dim=c(k,C))
  u_0 = seq(C)
  u = t(array(rep(u_0,k),dim=c(k,C)))
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
    m = sapply(1:k,(function(j){mu[j,components[i]]}))
    init <- rbind(init,mvrnorm(1,m,cov))
  }
  
  resp <- rep(0,C)
  rewards <- matrix(NA,nrow=t,ncol=k)
  rewards[1:n,1:k] <- init
  n = rep(n,k)
  regret <- c(0)
  xbar <- apply(rewards, MARGIN=2,mean,na.rm=TRUE)
  s <- apply(rewards,MARGIN=2,sd,na.rm=TRUE)
  for(i in 1:t){
    
    # sampling from posterior for each of the k arms
    theta.hat <- sapply(1:k,function(a){sum(gam[a,]*sapply(1:C,function(m)sample_from_comp(a,m)))/sum(gam[a,])})
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
    # Posterior update
    old_resp=c(10,10)
    while(niter < 50 & max(abs(old_resp-resp))>0.01){
      old_resp = resp
      print(alpha)
      print(beta)
      print(u)
      print(V_inv)
      print(gam)
      for(comp in 1:C){
        resp[comp] <- -0.5*(log(beta[j,comp]-digamma(alpha[j,comp])-1/V_inv[j,comp]+(obs.reward-u[j,comp])^2*alpha[j,comp]/beta[j,comp]))+digamma(gam[j,comp])-digamma(sum(gam[j,]))
      }
      
      #Normalizing step
      resp = resp-min(resp)
      Z = sum(exp(resp))
      resp = resp-log(Z)
      print("----")
      print(resp)
      
      for(comp in 1:C){
        gam[j,comp]<-gam[j,comp]+exp(resp[comp])
        #V_inv_old=V_inv[j,comp]
        V_inv[j,comp]<-V_inv0+exp(resp[comp])
        #u_old = u[j,comp]
        u[j,comp] <- 1/V_inv[j,comp]*(exp(resp[comp])*obs.reward+V_inv0*u_0[comp])
        alpha[j,comp]<-alpha_0+0.5*exp(resp[comp])
        beta[j,comp]<-beta_0+0.5*obs.reward^2*exp(resp[comp])+0.5*u_0[comp]^2*V_inv0-0.5*u[j,comp]^2*V_inv[j,comp]
      }
      niter <- niter + 1
      #print(niter)
    }
  }
  return(rewards)
}
