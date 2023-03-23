#############################################################################################
##################### METHODS TO CONSTRUCT CONFIDENCE INTERVALS #############################
#############################################################################################

# based on https://arxiv.org/pdf/1807.03582.pdf?fbclid=IwAR0RJA-eIuLOx3QRyKVszBFFPAnpE6Gh7h9MZtgvN_7c58PEiZb0EQrjY2c

# 3.1 Frequentist interval for p_hat

#--------------------------------------------------------------------------------
# exact Clopper-Pearson confidence interval for the relative frequency estimator
# for p in a Binomial model
#--------------------------------------------------------------------------------

frequentist.ci.binom <- function(n, k, alpha) {
  
  if (k == 0) {
    pl <- 0.0
    pu <- 1 - (alpha/2)^(1/n)
  }
  
  else if (k == n) {
    pl <- (alpha/2)^(1/n)
    pu <- 1.0
    }
  
  else {
    
    helper <- function(p, k, n, val) {
    return (pbinom(k, n, p) - val)
      }
    
    r <- uniroot(helper, k=(k-1),n=n, val=1-alpha/2,interval=c(0,1))
    pu <- r$rootr <- uniroot(helper, k=k,n=n, val=alpha/2,interval=c(0,1))
    pl <- r$root

    }
  
  return (data.frame(pl=pl, pu=pu))
}

# example 

CI <- frequentist.ci.binom(100, 5, 0.05) # 32 successes in 50 trials, we want a 95% CI
c( CI$pl, CI$pu.root)

CI$pl
#-----------------------------------------------------------------------
# Wilson approximate confidence interval for p when n is large
#-----------------------------------------------------------------------

approximate.ci.binom <- function(n, k, alpha) {
  
  if (k == 0) {
    pl <- 0.0
    sd <- sqrt(   ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pu <- ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  (  ( qnorm(1-alpha/2) / (2*n)) + qnorm(1-alpha/2) * sd )
  }
  
  else if (k == n) {
    sd <- sqrt(   ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pl <- ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) - qnorm(1-alpha/2) * sd )
    pu <- 1.0
  }
  
  else {
    
    sd <- sqrt( (( (k/n)* (1 - (k/n))  )/n) +  ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pl <- max(0, ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) - qnorm(1-alpha/2) * sd ) )
    pu <- min(1, ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) + qnorm(1-alpha/2) * sd ) )
    
  }
  
  return (data.frame(pl=pl, pu=pu))
}

# example

CI2 <-approximate.ci.binom(50, 32, 0.05)  # 32 successes in 50 trials, we want a 95% CI
CI2

# comparison of both frequentist CI when n is large and p is about the center of [0,1]

CI <-frequentist.ci.binom(10000, 4800, 0.05)  # 32 successes in 50 trials, we want a 95% CI
c( CI$p2, CI$p1.root)

CI2 <-approximate.ci.binom(10000, 4800, 0.05)  # 32 successes in 50 trials, we want a 95% CI
CI2



# 3.2 Likelihood ratio for p_hat

#--------------------------------------------------------------------
#  likelihood  ratio support  interval  for  the  relative  frequency
#--------------------------------------------------------------------


likratio.ci.binom <- function(n, k, K) {
  
  helper <- function(p, n, k, K) {
    return (p**k*(1-p)**(n-k) /((k/n)**k*(1-k/n)**(n-k))- 1/K)
  }
  
  pl <- rep(0,length(k))
  pu <- pl 
  
  if (k==0) {
    
    p1 <- 0
    
  } else {
    
    r <- uniroot(helper,n=n,k=k,K=K,interval=c(0,k/n))
    pl <- r$root
  }
  
  if (k==n) 
    {
    pu <- 1
    
    } else {
      
    r <- uniroot(helper,n=n,k=k,K=K,interval=c(k/n,1))
    pu <- r$root
  }
  
  return (data.frame(pl=pl, pu=pu))
  
}


# example

SCI <- likratio.ci.binom(50, 32, 8)  # 32 successes in 50 trials when the threshold is 8 (usual)




# 3.3 (1 - alpha) highest posterior density interval for p_hat

#--------------------------------------------------------------------
#  Highest posterior density  interval  for  the  relative  frequency
#--------------------------------------------------------------------

library(HDInterval)

CI <- hdi(qbeta, 1-alpha,shape1=(k+1),shape2=(n-k+1))  # n=50, k=32, alpha=0.05

pl <- CI[[1]] ; pu <- CI[[2]]

CI <- c(pl, pu)
CI




## 7 Coverage probability for the relative frequency - simulations
# 7.1 simulate samples of size n and compute an estimate p_hat of p
# 7.2 compute a confidence interval for p_hat
# 7.3 count the number of times the true value p falls into the given interval

set.seed(1986)
data <- rbinom(20, 1, 0.5)

# maximum likelihood estimation
mle = optimize(function(theta){sum(dbinom(x = data, size = 1, prob = theta, log = TRUE))},
               interval = c(0, 1),
               maximum = TRUE,
               tol = .Machine$double.eps^0.5)

p_mle = mle$maximum
p_mle # [1] 0.45

n = 20 # sample size
k = 9 # number of successes  (n*p_mle)
alpha = 0.05  
p = k/n   # true value = 0.45
p_hat = numeric(1000) # Number of replications: 1,000
ci_lower = numeric(1000)
ci_upper = numeric(1000)
ci = numeric(1000)
cr = 0

set.seed(1986)
for(i in 1:1000){
  
    p_hat[i] = mean(rbinom(n, 1, p))  # mle
    
    # compute Wilson approximate CI (based on normal distribution)
    ci_lower[i] = approximate.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pl
    ci_upper[i] = approximate.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pu
    
    if (p >= ci_lower[i] & p <= ci_upper[i])  cr=cr+1
}

cr/1000  # [1] 0.959
mean(ci_lower) # [1] 0.2458382
mean(ci_upper) # [1] 0.6705551


set.seed(1986)
for(i in 1:1000){
  
  p_hat[i] = mean(rbinom(n, 1, p))  # mle
  
  # compute exact Clopper-Pearson CI 
  ci_lower[i] = frequentist.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pl
  ci_upper[i] = frequentist.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pu.root
  
  if (p >= ci_lower[i] & p <= ci_upper[i])  cr=cr+1
}

cr/1000  # 1] 0.959
mean(ci_lower) # [1] 0.2390293
mean(ci_upper) # [1] 0.6833242


set.seed(1986)
for(i in 1:1000){
  
  p_hat[i] = mean(rbinom(n, 1, p))  # mle
  
  # compute high posterior density interval
  ci_lower[i] =  as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i])+1),shape2=(n-(n*p_hat[i])+1)))[1]
  ci_upper[i] =  as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i])+1),shape2=(n-(n*p_hat[i])+1)))[2]
  
  if (p >= ci_lower[i] & p <= ci_upper[i])  cr=cr+1
}

cr/1000  # 1] 0.959
mean(ci_lower) # [1] 0.2628083
mean(ci_upper) # [1] 0.6558001








######################################################################################
####### for R forum

## 7 Coverage probability for the relative frequency - simulations
# 7.1 simulate samples of size n and compute an estimate p_hat of p
# 7.2 compute a confidence interval for p_hat
# 7.3 count the number of times the true value p falls into the given interval

set.seed(1986)
data <- rbinom(20, 1, 0.5)

# maximum likelihood estimation
mle = optimize(function(theta){sum(dbinom(x = data, size = 1, prob = theta, log = TRUE))},
               interval = c(0, 1),
               maximum = TRUE,
               tol = .Machine$double.eps^0.5)

p_mle = mle$maximum
p_mle # [1] 0.45

n = 20 # sample size
k = 9 # number of successes  (n*p_mle)
alpha = 0.05  
p = k/n   # true value = 0.45
p_hat = numeric(1000) # Number of replications: 1,000
ci_lower = numeric(1000)
ci_upper = numeric(1000)
ci = numeric(1000)
cr = 0

set.seed(1986)
for(i in 1:1000){
  
  p_hat[i] = mean(rbinom(n, 1, p)) # mle
  
  # compute Wilson approximate CI (based on normal distribution)
  ci_lower[i] = approximate.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pl
  ci_upper[i] = approximate.ci.binom(n=n, k= n*p_hat[i], alpha=alpha)$pu
  
  if (p >= ci_lower[i] & p <= ci_upper[i])  cr=cr+1
}

cr/1000  # [1] 0.959
mean(ci_lower) # [1] 0.2458382
mean(ci_upper) # [1] 0.6705551

# plot
qplot(ci_lower, geom="histogram", fill=..count.., bins = 10) +
  theme_light()

n = 20 # sample size
k = 9 # number of successes  (n*p_mle)
alpha = 0.05  
p = k/n   # true value = 0.45
p_hat = numeric(1000) # Number of replications: 1,000
ci_lower = numeric(1000)
ci_upper = numeric(1000)
ci = numeric(1000)
cr = 0


set.seed(1986)
for(i in 1:1000){
  
  p_hat[i] = mean(rbinom(n, 1, p)) # mle
  
  # compute high posterior density interval
  ci_lower[i] =  as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i])+1),shape2=(n-(n*p_hat[i])+1)))[1]
  ci_upper[i] =  as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i])+1),shape2=(n-(n*p_hat[i])+1)))[2]
  
  if (p >= ci_lower[i] & p <= ci_upper[i])  cr=cr+1
}

cr/1000  # 1] 0.959
mean(ci_lower) # [1] 0.2628083
mean(ci_upper) # [1] 0.6558001

# plot
qplot(ci_upper, geom="histogram", fill=..count.., bins = 10) +
  theme_light()





#---------------------------
# my study # experiment 1
#---------------------------

library(ggplot2)
library(HDInterval)

## 7 Coverage probability for the relative frequency - simulations
# 7.1 simulate samples of size n and compute an estimate p_hat of p
# 7.2 compute a confidence interval for p_hat
# 7.3 count the number of times the true value p falls into the given interval

# function to compute Wilson confidence intervals
approximate.ci.binom <- function(n, k, alpha) {
  
  if (k == 0) {
    pl <- 0.0
    sd <- sqrt(   ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pu <- ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  (  ( qnorm(1-alpha/2) / (2*n)) + qnorm(1-alpha/2) * sd )
  }
  else if (k == n) {
    sd <- sqrt(   ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pl <- ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) - qnorm(1-alpha/2) * sd )
    pu <- 1.0
  }
  else {
    sd <- sqrt( (( (k/n)* (1 - (k/n))  )/n) +  ( ((qnorm(1-alpha/2))^2) / (4*n^2)  )  )
    pl <- max(0, ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) - qnorm(1-alpha/2) * sd ) )
    pu <- min(1, ( 1 / (1 + (qnorm(1-alpha/2)  / n) ) )  *  ( (k/n)    + ( qnorm(1-alpha/2) / (2*n)) + qnorm(1-alpha/2) * sd ) )
  }
  return (data.frame(pl=pl, pu=pu))
}


## ECR Wilson approximate interval

n = 100 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerW = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperW = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crW = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}
  p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

  for(i in 1:10000){
    for(j in 1:(n+1)) {
      
  # compute Wilson approximate CI (based on normal distribution)
  ci_lowerW[i,j] = approximate.ci.binom(n = n, k = n*p_hat[i,j], alpha=alpha)$pl
  ci_upperW[i,j] = approximate.ci.binom(n = n, k = n*p_hat[i,j], alpha=alpha)$pu
  
    if (p[i,j] >= ci_lowerW[i,j] & p[i,j] <= ci_upperW[i,j])  crW[j]=crW[j]+1 
    }
}

# coverage rate
crW/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crW/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Wilson approximate Confidence Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()



# average length of the confidence intervals
CI_lengthW=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthW[a, ] = apply(cbind(ci_lowerW[,a], ci_upperW[,a]),2, mean)
}

CI_lengthW

av_lengthW <- CI_lengthW[,2] - CI_lengthW[,1]
min(av_lengthW)
max(av_lengthW)
mean(av_lengthW)  # [1] 0.1534113

# Plot of average length
qplot(av_lengthW, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Histogram of average length interval - Approximate Wilson CI (0.1534)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.01), 
                     labels = seq(0.02,0.22, by=0.01)) + 
  theme_light()


## ECR Bayesian HPD interval

n = 100 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerB = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperB = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crB = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute high posterior density interval
    ci_lowerB[i,j] = as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i,j])+1),shape2=(n-(n*p_hat[i,j])+1)))[1]
    ci_upperB[i,j] = as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i,j])+1),shape2=(n-(n*p_hat[i,j])+1)))[2]
    
    if (p[i,j] >= ci_lowerB[i,j] & p[i,j] <= ci_upperB[i,j])  crB[j]=crB[j]+1 
  }
}

# coverage rate
crB/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crB/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Bayesian HPD Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthB=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthB[a, ] = apply(cbind(ci_lowerB[,a], ci_upperB[,a]),2, mean)
}

CI_lengthB

av_lengthB <- CI_lengthB[,2] - CI_lengthB[,1]
min(av_lengthB)
max(av_lengthB)
mean(av_lengthB) # [1] 0.14999

# Plot of average length
qplot(av_lengthB, geom="histogram", fill=..count.., bins = length(seq(0.01,0.24, by=0.01))) +
  ggtitle("Histogram of average length interval - Bayesian HPD Interval (0.1499)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.01,0.22, by=0.01), 
                     labels = seq(0.01,0.22, by=0.01)) + 
  theme_light()


## ECR Likelihood Ratio Support Interval

n = 100 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerL = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperL = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crL = numeric((n+1))
K=8 #

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute likelihood support ratio (K=8) support interval
    ci_lowerL[i,j] = likratio.ci.binom(n=n, k=n*p_hat[i,j], K=K)$pl
    ci_upperL[i,j] = likratio.ci.binom(n=n, k=n*p_hat[i,j], K=K)$pu
    
    if (p[i,j] >= ci_lowerL[i,j] & p[i,j] <= ci_upperL[i,j])  crL[j]=crL[j]+1 
  }
}

# coverage rate
crL/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crL/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Likelihood Ratio Support Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthL=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthL[a, ] = apply(cbind(ci_lowerL[,a], ci_upperL[,a]),2, mean)
}

CI_lengthL

av_lengthL <- CI_lengthL[,2] - CI_lengthL[,1]
min(av_lengthL)
max(av_lengthL)
mean(av_lengthL) # [1] 0.1563925

# Plot of average length
qplot(av_lengthL, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Histogram of average length interval - Likelihood Ratio Support Interval (0.1563)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.01), 
                     labels = seq(0.02,0.22, by=0.01)) + 
  theme_light()



## ECR Exact Clopper-Pearson (frequentist) CI

n = 100 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerE = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperE = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crE = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute exact frequentist CI
    ci_lowerE[i,j] = as.numeric(frequentist.ci.binom(n=n, k=n*p_hat[i,j], alpha=alpha))[1]
    ci_upperE[i,j] = as.numeric(frequentist.ci.binom(n=n, k=n*p_hat[i,j], alpha=alpha))[2]
    
    if (p[i,j] >= ci_lowerE[i,j] & p[i,j] <= ci_upperE[i,j])  crE[j]=crE[j]+1 
  }
}

# coverage rate
crE/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crE/rep(10000,(n+1))
qplot(x,y,geom='smooth',span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Exact Clopper-Pearson Confidence Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthE=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthE[a, ] = apply(cbind(ci_lowerE[,a], ci_upperE[,a]),2, mean)
}

av_lengthE <- CI_lengthE[,2] - CI_lengthE[,1]
min(av_lengthE)
max(av_lengthE)
mean(av_lengthE) # [1] 0.1601739

# Plot of average length
qplot(av_lengthE, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Histogram of average length interval - Exact Clopper-Pearson CI (0.1601) ") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.01), 
                     labels = seq(0.02,0.22, by=0.01)) + 
  theme_light()



############################
## multiple plots in ggplot2
############################

# ECR

x = k/n

g1 <- qplot(x,crE/rep(10000,(n+1)),geom='smooth',span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Exact Clopper-Pearson CI") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g2 <- qplot(x, crW/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Wilson approximate CI") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g3 <- qplot(x,crL/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Likelihood Ratio Support Interval (0.1601)") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g4 <- qplot(x,crB/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Bayesian HPD Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()

library(gridExtra)
grid.arrange(g1, g2, g3, g4, nrow = 2)



# AIL

p1 <- qplot(av_lengthE, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Exact Clopper-Pearson CI  (0.1601)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p2 <- qplot(av_lengthW, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Approximate Wilson CI (0.1534)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p3 <- qplot(av_lengthL, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Likelihood Ratio Support Interval (0.1563)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p4 <-qplot(av_lengthB, geom="histogram", fill=..count.., bins = length(seq(0.01,0.24, by=0.01))) +
  ggtitle("Average length - Bayesian HPD Interval (0.1499)") +
  scale_x_continuous(name = "Average length", breaks = seq(0.01,0.22, by=0.02), 
                     labels = seq(0.01,0.22, by=0.02)) + 
  theme_light()

grid.arrange(p1, p2, p3, p4, nrow = 2)




#---------------------------
# my study :) # experiment 2
#---------------------------

## 7 Coverage probability for the relative frequency - simulations
# 7.1 simulate samples of size n and compute an estimate p_hat of p
# 7.2 compute a confidence interval for p_hat
# 7.3 count the number of times the true value p falls into the given interval



## ECR Wilson approximate interval

n = 15 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerW = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperW = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crW = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}
p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute Wilson approximate CI (based on normal distribution)
    ci_lowerW[i,j] = approximate.ci.binom(n = n, k = n*p_hat[i,j], alpha=alpha)$pl
    ci_upperW[i,j] = approximate.ci.binom(n = n, k = n*p_hat[i,j], alpha=alpha)$pu
    
    if (p[i,j] >= ci_lowerW[i,j] & p[i,j] <= ci_upperW[i,j])  crW[j]=crW[j]+1 
  }
}

# coverage rate
crW/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crW/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Wilson approximate Confidence Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()



# average length of the confidence intervals
CI_lengthW=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthW[a, ] = apply(cbind(ci_lowerW[,a], ci_upperW[,a]),2, mean)
}

CI_lengthW

av_lengthW <- CI_lengthW[,2] - CI_lengthW[,1]
min(av_lengthW)
max(av_lengthW)
mean(av_lengthW)  # [1] 0.3806957

# Plot of average length
qplot(av_lengthW, geom="histogram", fill=..count.., bins = length(seq(min(av_lengthW),max(av_lengthW), by=0.05)+2)) +
  ggtitle("Histogram of average length interval - Approximate Wilson CI") +
  scale_x_continuous(name = "Average length", breaks = seq(min(av_lengthW),max(av_lengthW), by=0.05), 
                     labels = seq(min(av_lengthW),max(av_lengthW), by=0.05)) + 
  theme_light()


## ECR Bayesian HPD interval

n = 20 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerB = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperB = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crB = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute high posterior density interval
    ci_lowerB[i,j] = as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i,j])+1),shape2=(n-(n*p_hat[i,j])+1)))[1]
    ci_upperB[i,j] = as.numeric(hdi(qbeta, 1-alpha,shape1=((n*p_hat[i,j])+1),shape2=(n-(n*p_hat[i,j])+1)))[2]
    
    if (p[i,j] >= ci_lowerB[i,j] & p[i,j] <= ci_upperB[i,j])  crB[j]=crB[j]+1 
  }
}

# coverage rate
crB/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crB/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Bayesian HPD Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthB=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthB[a, ] = apply(cbind(ci_lowerB[,a], ci_upperB[,a]),2, mean)
}

CI_lengthB

av_lengthB <- CI_lengthB[,2] - CI_lengthB[,1]
min(av_lengthB)
max(av_lengthB)
mean(av_lengthB) # [1] 0.3092311

# Plot of average length
qplot(av_lengthB, geom="histogram", fill=..count.., bins = length(seq(0.01,0.24, by=0.01))) +
  ggtitle("Histogram of average length interval - Bayesian HPD Interval") +
  scale_x_continuous(name = "Average length", breaks = seq(0.01,0.22, by=0.01), 
                     labels = seq(0.01,0.22, by=0.01)) + 
  theme_light()


## ECR Likelihood Ratio Support Interval

n = 20 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerL = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperL = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crL = numeric((n+1))
K=8 #

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute likelihood support ratio (K=8) support interval
    ci_lowerL[i,j] = likratio.ci.binom(n=n, k=n*p_hat[i,j], K=K)$pl
    ci_upperL[i,j] = likratio.ci.binom(n=n, k=n*p_hat[i,j], K=K)$pu
    
    if (p[i,j] >= ci_lowerL[i,j] & p[i,j] <= ci_upperL[i,j])  crL[j]=crL[j]+1 
  }
}

# coverage rate
crL/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crL/rep(10000,(n+1))
qplot(x,y, geom='smooth', span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Likelihood Ratio Support Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthL=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthL[a, ] = apply(cbind(ci_lowerL[,a], ci_upperL[,a]),2, mean)
}

CI_lengthL

av_lengthL <- CI_lengthL[,2] - CI_lengthL[,1]
min(av_lengthL)
max(av_lengthL)
mean(av_lengthL) # [1] 0.322461

# Plot of average length
qplot(av_lengthL, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Histogram of average length interval - Likelihood Ratio Support Interval") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.01), 
                     labels = seq(0.02,0.22, by=0.01)) + 
  theme_light()



## ECR Exact Clopper-Pearson (frequentist) CI

n = 20 # sample size
k = 0:n # number of successes  
alpha = 0.05  
p = matrix(rep(k/n,10000), ncol=length(k/n), byrow = TRUE)   
p_hat = matrix(rep(0,10000*(n+1)), ncol=(n+1)) # Number of replications: 1,000
ci_lowerE = matrix(rep(0,10000*(n+1)), ncol=(n+1))
ci_upperE = matrix(rep(0,10000*(n+1)), ncol=(n+1))
crE = numeric((n+1))

set.seed(1986)
for(i in 1:10000){
  for(j in 0:n) {
    p_hat[i,j] = mean(rbinom(n, 1, j/n)) # mle
  } 
}

p_hat = cbind(p_hat[,(n+1)], p_hat[,1:n]) # reorder columns

for(i in 1:10000){
  for(j in 1:(n+1)) {
    
    # compute exact frequentist CI
    ci_lowerE[i,j] = as.numeric(frequentist.ci.binom(n=n, k=n*p_hat[i,j], alpha=alpha))[1]
    ci_upperE[i,j] = as.numeric(frequentist.ci.binom(n=n, k=n*p_hat[i,j], alpha=alpha))[2]
    
    if (p[i,j] >= ci_lowerE[i,j] & p[i,j] <= ci_upperE[i,j])  crE[j]=crE[j]+1 
  }
}

# coverage rate
crE/rep(10000,(n+1))  

# plot with ggplot2
x = k/n
y = crE/rep(10000,(n+1))
qplot(x,y,geom='smooth',span =0.15, color="darkred") +
  ggtitle("Empirical Coverage Rate behaviour as a function of p - Exact Clopper-Pearson Confidence Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


# average length of the confidence intervals
CI_lengthE=matrix(rep(0, 2*(n+1)), ncol=2)

for(a in 1:(n+1)) {
  CI_lengthE[a, ] = apply(cbind(ci_lowerE[,a], ci_upperE[,a]),2, mean)
}

av_lengthE <- CI_lengthE[,2] - CI_lengthE[,1]
min(av_lengthE)
max(av_lengthE)
mean(av_lengthE) # [1] 0.356273

# Plot of average length
qplot(av_lengthE, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Histogram of average length interval - Exact Clopper-Pearson CI") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.01), 
                     labels = seq(0.02,0.22, by=0.01)) + 
  theme_light()



############################
## multiple plots in ggplot2
############################

# ECR

x = k/n

g1 <- qplot(x,crE/rep(10000,(n+1)),geom='smooth',span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Exact Clopper-Pearson CI") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g2 <- qplot(x, crW/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Wilson approximate CI") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g3 <- qplot(x,crL/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Likelihood Ratio Support Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()


g4 <- qplot(x,crB/rep(10000,(n+1)), geom='smooth', span =0.1, color="darkred") +
  ggtitle("ECR behaviour - Bayesian HPD Interval") +
  xlab("p") + ylab("Empirical Coverage Rate") +
  ylim(0.9, 1) +
  theme_light()

library(gridExtra)
grid.arrange(g1, g2, g3, g4, nrow = 2)



# AIL

p1 <- qplot(av_lengthE, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Exact Clopper-Pearson CI") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p2 <- qplot(av_lengthW, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Approximate Wilson CI") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p3 <- qplot(av_lengthL, geom="histogram", fill=..count.., bins = length(seq(0.02,0.24, by=0.01)+2)) +
  ggtitle("Average length - Likelihood Ratio Support Interval") +
  scale_x_continuous(name = "Average length", breaks = seq(0.02,0.22, by=0.02), 
                     labels = seq(0.02,0.22, by=0.02)) + 
  theme_light()

p4 <-qplot(av_lengthB, geom="histogram", fill=..count.., bins = length(seq(0.01,0.24, by=0.01))) +
  ggtitle("Average length - Bayesian HPD Interval") +
  scale_x_continuous(name = "Average length", breaks = seq(0.01,0.22, by=0.02), 
                     labels = seq(0.01,0.22, by=0.02)) + 
  theme_light()

grid.arrange(p1, p2, p3, p4, nrow = 2)