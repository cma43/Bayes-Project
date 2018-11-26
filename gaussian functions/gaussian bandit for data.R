# alpha is parameter defining prior
# X is data. Observations from each arm are in a single column

# Need a measure equivalent to regret since sampling distribution is not specified

gaussian_bandit <- function(alpha,X,t,k){

# initializing n since prior is improper
n <- max( c(2,2-ceiling(2*alpha)))
# initial rewards to get a proper posterior.
init <- apply(X,2, sample, size=n, replace=FALSE)

rewards <- matrix(NA,nrow=t,ncol=k)
rewards[1:n,1:k] <- init

xbar <- apply(rewards, MARGIN=2,mean,na.rm=TRUE)
s <- apply(rewards,MARGIN=2,sd,na.rm=TRUE)

for(i in 1:t){
# sampling from each arm
theta.hat <- sapply(1:k, function(l) rt(1,df=n[l]+2*alpha-1)
                    *sqrt(s[l]/(n[l]*(n[l]+2*alpha-1)))+xbar[l])
# finding best arm and randomly breaking ties
j <- sample(which(theta.hat ==max(theta.hat)),1)

#update everything about posterior
n[j] = n[j]+1
rewards[n[j],j] <- obs.reward
xbar[j] <- mean(rewards[,j], na.rm=TRUE)
s[j] <- sd(rewards[,j],na.rm=TRUE)
}
# What we want to return depends on analysis. Add this later.
}