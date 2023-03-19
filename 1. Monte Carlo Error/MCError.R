#-----------------------------------------------------------------------
# On the assessment of MC error in simulation-based statistical analyses
# partial reproductio
#-----------------------------------------------------------------------

------------------------------------
###### section 2 - MONTE CARLO ERROR
------------------------------------
  
  
### experiment 1

# we want to have an estimation of 3 operating caracteristics (and their standard deviation
# refered to as MCE) for the estimate of beta1, MLE estimate in logistic regression.

# the 3 operating caracteristics are percent bias, coverage probability of 95% CI and statistical power
# each MLE estimate and the operating caracteristics are computed on a sample of size N = 100
# R (number of replicates) = 100, 500, 1000, 2500, 5000, 10000
# I use M = 1,000 for convenience (time)

install.packages("arm")
library(arm) # for extraction of standard error estimates of MLE estimates 
betaA=numeric(100)
betaB=numeric(500)
betaC=numeric(1000)
betaD=numeric(2500)
betaE=numeric(5000)
betaF=numeric(10000)

sehatbetaA=numeric(100)
sehatbetaB=numeric(500)
sehatbetaA=numeric(100)
sehatbetaC=numeric(1000)
sehatbetaD=numeric(2500)
sehatbetaE=numeric(5000)
sehatbetaF=numeric(10000)

# 1. Number of replications: 100

set.seed(3)
phi_hat_c_r100=numeric(1000)
phi_hat_b_r100=numeric(1000)
phi_hat_p_r100=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:100)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaA[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaA[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaA) - qnorm(0.975) * mean(sehatbetaA)) 
    ci_upperR1= (mean(betaA) + qnorm(0.975) * mean(sehatbetaA))
  }
  phi_hat_c_r100[j] = ((1/100) * sum(betaA >= rep(ci_lowerR1,100) & betaA <= rep(ci_upperR1,100))) *100
  phi_hat_b_r100[j] <- ((1/100) * sum((betaA - (betax*rep(1, 100))) / (betax*rep(1, 100)))) *100
  phi_hat_p_r100[j] = ((1/100) * sum( (betaA/sehatbetaA) > qnorm(0.975))) *100
}


min(phi_hat_b_r100)   # Percent bias - min
# [1] -22.47584
max(phi_hat_b_r100)   # Percent bias - max
# [1] 24.42148
mean(phi_hat_b_r100)  # Percent bias - mean
# [1] 1.025418
sd(phi_hat_b_r100)    # Percent bias - MCE
# [1] 6.819096

min(phi_hat_c_r100)   # Coverage rate - min
# [1] 85
max(phi_hat_c_r100)   # Coverage rate - max
# [1] 100
mean(phi_hat_c_r100)  # Coverage rate- mean
# [1] 94.65
sd(phi_hat_c_r100)    # Coverage rate - MCE
# [1] 2.227209

min(phi_hat_p_r100)   # Power - min
# [1] 20
max(phi_hat_p_r100)   # Power - max
# [1] 50
mean(phi_hat_p_r100)  # Power- mean
# [1] 33.156
sd(phi_hat_p_r100)    # Power - MCE
# [1] 4.675631


# 2. Number of replications: 500

set.seed(3)
phi_hat_c_r500=numeric(1000)
phi_hat_b_r500=numeric(1000)
phi_hat_p_r500=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:500)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaB[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaB[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaB) - qnorm(0.975) * mean(sehatbetaB)) 
    ci_upperR1= (mean(betaB) + qnorm(0.975) * mean(sehatbetaB))
  }
  phi_hat_c_r500[j] =  ((1/500) * sum(betaB >= rep(ci_lowerR1,500) & betaB <= rep(ci_upperR1,500))) *100
  phi_hat_b_r500[j] <- ((1/500) * sum((betaB - (betax*rep(1, 500))) / (betax*rep(1, 500)))) *100
  phi_hat_p_r500[j] =  ((1/500) * sum( (betaB/sehatbetaB) > qnorm(0.975))) *100
  mceCLT1001[j]<-(1/500)* sqrt(sum((phi_hat_b_r100-((1/500) * sum((betaB - (betax*rep(1, 100))) / (betax*rep(1, 100)))) *100)^2))
}


min(phi_hat_b_r500)   # Percent bias - min
# [1] -10.0406
max(phi_hat_b_r500)   # Percent bias - max
# [1] 11.75664
mean(phi_hat_b_r500)  # Percent bias - mean
# [1] 0.9248017
sd(phi_hat_b_r500)    # Percent bias - MCE
# [1] 3.213109

min(phi_hat_c_r500)   # Coverage rate - min
# [1] 90.8
max(phi_hat_c_r500)   # Coverage rate - max
# [1] 97.4
mean(phi_hat_c_r500)  # Coverage rate- mean
# [1] 94.515
sd(phi_hat_c_r500)    # Coverage rate - MCE
# [1] 0.9924717

min(phi_hat_p_r500)   # Power - min
# [1] 26.2
max(phi_hat_p_r500)   # Power - max
# [1] 39.4
mean(phi_hat_p_r500)  # Power- mean
# [1] 32.9866
sd(phi_hat_p_r500)    # Power - MCE
# [1] 2.162549


# 3. Number of replications: 1000

set.seed(3)
phi_hat_c_r1000=numeric(1000)
phi_hat_b_r1000=numeric(1000)
phi_hat_p_r1000=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:1000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaC[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaC[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaC) - qnorm(0.975) * mean(sehatbetaC)) 
    ci_upperR1= (mean(betaC) + qnorm(0.975) * mean(sehatbetaC))
  }
  phi_hat_c_r1000[j] = ((1/1000) * sum(betaC >= rep(ci_lowerR1,1000) & betaC <= rep(ci_upperR1,1000))) *100
  phi_hat_b_r1000[j] <- ((1/1000) * sum((betaC - (betax*rep(1, 1000))) / (betax*rep(1, 1000)))) *100
  phi_hat_p_r1000[j] = ((1/1000) * sum( (betaC/sehatbetaC) > qnorm(0.975))) *100
}


min(phi_hat_b_r1000)   # Percent bias - min
# [1] -6.913045
max(phi_hat_b_r1000)   # Percent bias - max
# [1] 7.740247
mean(phi_hat_b_r1000)  # Percent bias - mean
# [1] 0.9809707
sd(phi_hat_b_r1000)    # Percent bias - MCE
# [1] 2.320492

min(phi_hat_c_r1000)   # Coverage rate - min
# [1] 92.3
max(phi_hat_c_r1000)   # Coverage rate - max
# [1] 96.5
mean(phi_hat_c_r1000)  # Coverage rate- mean
# [1] 94.5522
sd(phi_hat_c_r1000)    # Coverage rate - MCE
# [1] 0.7119987

min(phi_hat_p_r1000)   # Power - min
# [1] 28.4
max(phi_hat_p_r1000)   # Power - max
# [1] 37.7
mean(phi_hat_p_r1000)  # Power- mean
# [1] 33.0423
sd(phi_hat_p_r1000)    # Power - MCE
# [1] 1.548709


# 4. Number of replications: 2500

set.seed(3)
phi_hat_c_r2500=numeric(1000)
phi_hat_b_r2500=numeric(1000)
phi_hat_p_r2500=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:2500)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaD[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaD[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaD) - qnorm(0.975) * mean(sehatbetaD)) 
    ci_upperR1= (mean(betaD) + qnorm(0.975) * mean(sehatbetaD))
  }
  phi_hat_c_r2500[j] = ((1/2500) * sum(betaD >= rep(ci_lowerR1,2500) & betaD <= rep(ci_upperR1,2500))) *100
  phi_hat_b_r2500[j] <- ((1/2500) * sum((betaD - (betax*rep(1, 2500))) / (betax*rep(1, 2500)))) *100
  phi_hat_p_r2500[j] = ((1/2500) * sum( (betaD/sehatbetaD) > qnorm(0.975))) *100
}


min(phi_hat_b_r2500)   # Percent bias - min
# [1] -3.706412
max(phi_hat_b_r2500)   # Percent bias - max
# [1] 4.798646
mean(phi_hat_b_r2500)  # Percent bias - mean
# [1] 0.9572071
sd(phi_hat_b_r2500)    # Percent bias - MCE
# [1] 1.360274

min(phi_hat_c_r2500)   # Coverage rate - min
# [1] 93.04
max(phi_hat_c_r2500)   # Coverage rate - max
# [1] 96
mean(phi_hat_c_r2500)  # Coverage rate- mean
# [1] 94.52676
sd(phi_hat_c_r2500)    # Coverage rate - MCE
# [1] 0.4794499

min(phi_hat_p_r2500)   # Power - min
# [1] 29.32
max(phi_hat_p_r2500)   # Power - max
# [1] 36.52
mean(phi_hat_p_r2500)  # Power- mean
# [1] 33.02532
sd(phi_hat_p_r2500)    # Power - MCE
# [1] 0.9555203

# 5. Number of replications: 5000

set.seed(3)
phi_hat_c_r5000=numeric(1000)
phi_hat_b_r5000=numeric(1000)
phi_hat_p_r5000=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:5000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaE[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaE[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaE) - qnorm(0.975) * mean(sehatbetaE)) 
    ci_upperR1= (mean(betaE) + qnorm(0.975) * mean(sehatbetaE))
  }
  phi_hat_c_r5000[j] = ((1/5000) * sum(betaE >= rep(ci_lowerR1,5000) & betaE <= rep(ci_upperR1,5000))) *100
  phi_hat_b_r5000[j] <- ((1/5000) * sum((betaE - (betax*rep(1, 5000))) / (betax*rep(1, 5000)))) *100
  phi_hat_p_r5000[j] = ((1/5000) * sum( (betaE/sehatbetaE) > qnorm(0.975))) *100
}


min(phi_hat_b_r5000)   # Percent bias - min
# [1] -2.378006
max(phi_hat_b_r5000)   # Percent bias - max
# [1] 3.627774
mean(phi_hat_b_r5000)  # Percent bias - mean
# [1] 0.9044886
sd(phi_hat_b_r5000)    # Percent bias - MCE
# [1] 0.9539884

min(phi_hat_c_r5000)   # Coverage rate - min
# [1] 93.22
max(phi_hat_c_r5000)   # Coverage rate - max
# [1] 95.52
mean(phi_hat_c_r5000)  # Coverage rate- mean
# [1] 94.49948
sd(phi_hat_c_r5000)    # Coverage rate - MCE
# [1] 0.3616309

min(phi_hat_p_r5000)   # Power - min
# [1] 30.7
max(phi_hat_p_r5000)   # Power - max
# [1] 35.26
mean(phi_hat_p_r5000)  # Power- mean
# [1] 32.98502
sd(phi_hat_p_r5000)    # Power - MCE
# [1] 0.6623857



# 6. Number of replications: 10000

set.seed(3)
phi_hat_c_r10000=numeric(1000)
phi_hat_b_r10000=numeric(1000)
phi_hat_p_r10000=numeric(1000)
for(j in 1:1000)  {
  for(i in 1:10000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    Px1=0.3  # P(X=1)
    Px0=0.7   # P(X=0)
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaF[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
    sehatbetaF[i]<-se.coef(glm( y1~X,family="binomial"))[2]
    ci_lowerR1= (mean(betaF) - qnorm(0.975) * mean(sehatbetaF)) 
    ci_upperR1= (mean(betaF) + qnorm(0.975) * mean(sehatbetaF))
  }
  phi_hat_c_r10000[j] = ((1/10000) * sum(betaF >= rep(ci_lowerR1,10000) & betaF <= rep(ci_upperR1,10000))) *100
  phi_hat_b_r10000[j] <- ((1/10000) * sum((betaF - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100
  phi_hat_p_r10000[j] = ((1/10000) * sum( (betaF/sehatbetaF) > qnorm(0.975))) *100
}


min(phi_hat_b_r10000)   # Percent bias - min
# [1] -1.521861
max(phi_hat_b_r10000)   # Percent bias - max
# [1] 3.087407
mean(phi_hat_b_r10000)  # Percent bias - mean
# [1] 0.898388
sd(phi_hat_b_r10000)    # Percent bias - MCE
# [1] 0.6571217

min(phi_hat_c_r10000)   # Coverage rate - min
# [1] 93.57
max(phi_hat_c_r10000)   # Coverage rate - max
# [1] 99.36
mean(phi_hat_c_r10000)  # Coverage rate- mean
# [1] 94.50169
sd(phi_hat_c_r10000)    # Coverage rate - MCE
# [1] 0.3211775

min(phi_hat_p_r10000)   # Power - min
# [1] 31.44
max(phi_hat_p_r10000)   # Power - max
# [1] 34.39
mean(phi_hat_p_r10000)  # Power- mean
# [1] 32.98765
sd(phi_hat_p_r10000)    # Power - MCE
# [1] 0.4596255

#---------
# figure 1
#---------

beta11=numeric(10000)
beta21=numeric(10000)
beta31=numeric(10000)
beta41=numeric(10000)
beta51=numeric(10000)

set.seed(45) 
seeds<-sample(1:1000, 5, replace=TRUE)
seeds
S1 <-seeds[1]
S2 <-seeds[2]
S3 <-seeds[3]
S4 <-seeds[4]
S5 <-seeds[5]

beta0= -1     # known parameter for logistic regression
betax = log(2)   # [1] 0.6931472  parameter for logistic regression
Px1=0.3  # P(X=1)
Px0=0.7   # P(X=0)
levels=c(1,0)
N1=100
X=c(rep(1,30),rep(0,70))  # fixed

set.seed(S1)
for(i in 1:10000)
{
  z1 = beta0 + betax *X      # linear combination with a bias
  Py1 = 1/(1+exp(-z1))   # pass through an inv-logit function
  y1 <- rbinom(N1,1,Py1)
  beta11[i] <- glm(y1 ~ X, family="binomial"(link='logit'))$coeff[2]
}

summary(beta11)  #
# mean percent bias
((1/10000) * sum((beta11 - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100

set.seed(S2)
for(i in 1:10000)
{
  z1 = beta0 + betax *X      # linear combination
  Py1 = 1/(1+exp(-z1))   # inv-logit function
  y1 <- rbinom(N1,1,Py1)
  beta21[i] <- glm(y1 ~ X, family="binomial"(link='logit'))$coeff[2]
}

summary(beta21)  
# mean percent bias
((1/10000) * sum((beta21 - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100

set.seed(S3)
for(i in 1:10000)
{
  z1 = beta0 + betax *X      # linear combination
  Py1 = 1/(1+exp(-z1))   # inv-logit function
  y1 <- rbinom(N1,1,Py1)
  beta31[i] <- glm(y1 ~ X, family="binomial"(link='logit'))$coeff[2]
}

summary(beta31)  
# mean percent bias
((1/10000) * sum((beta31 - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100

set.seed(S4)
for(i in 1:10000)
{
  z1 = beta0 + betax *X      # linear combination 
  Py1 = 1/(1+exp(-z1))   # inv-logit function
  y1 <- rbinom(N1,1,Py1)
  beta41[i] <- glm(y1 ~ X, family="binomial"(link='logit'))$coeff[2]
}

summary(beta41)  # 
# mean percent bias
((1/10000) * sum((beta41 - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100


set.seed(S5)
for(i in 1:10000)
{
  z1 = beta0 + betax *X      # linear combination 
  Py1 = 1/(1+exp(-z1))   # inv-logit function
  y1 <- rbinom(N1,1,Py1)
  beta51[i] <- glm(y1 ~ X, family="binomial"(link='logit'))$coeff[2]
}

summary(beta51)  
# mean percent bias
((1/10000) * sum((beta51 - (betax*rep(1, 10000))) / (betax*rep(1, 10000)))) *100


# figure 1

den=1:10000
hplotR61=cumsum(round((round(beta11,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den
hplotR62=cumsum(round((round(beta21,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den
hplotR63=cumsum(round((round(beta31,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den
hplotR64=cumsum(round((round(beta41,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den
hplotR65=cumsum(round((round(beta51,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den
plot(hplotR61,type="l", xlim=c(1,10000), ylim=c(-5,5), main="MC estimates of percent bias for the MLE"
     ,xlab="Number of replicates, R", ylab="Percent bias")
lines(hplotR62)
lines(hplotR63)
lines(hplotR64)
lines(hplotR65)

# max, min, range
beta_max<- max((cumsum(round((round(beta11,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta21,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta31,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta41,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta51,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000])

beta_min<- min((cumsum(round((round(beta11,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta21,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta31,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta41,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000]
               ,(cumsum(round((round(beta51,2) - (betax*rep(1, 10000))) / (betax*rep(1, 10000)),2)*100) / den)[10000])

beta_max - beta_min  # range
# [1] 2.2058


#-----------------------------------------------------
###### section 4 - QUANTIFICATION OF MONTE CARLO ERROR
#-----------------------------------------------------

### experiment 2

###############################################################################
###############################################################################
## table 3 -  MCE associated with estimation of characteristics of MC sampling
## distribution for percent bias and mean

betaA2=numeric(100)
betaB2=numeric(500)
betaC2=numeric(1000)
betaD2=numeric(2500)
betaE2=numeric(5000)
betaF2=numeric(10000)
betamean100=numeric(1000)
betamean500=numeric(1000)
betamean1000=numeric(1000)
betamean2500=numeric(1000)
betamean5000=numeric(1000)
betamean10000=numeric(1000)

# R=100
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:100)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaA2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean100[j] = mean(betaA2)
}

# R=500
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:500)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaB2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean500[j] = mean(betaB2)
}

# R=1000
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:1000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaC2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean1000[j] = mean(betaC2)
}

# R=2500
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:2500)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaD2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean2500[j] = mean(betaD2)
}

# # R=5000
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:5000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaE2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean5000[j] = mean(betaE2)
}

# R=10000
set.seed(3)
for(j in 1:1000)  {
  for(i in 1:10000)
  { 
    beta0= -1     # known parameter for logistic regression
    betax = log(2)   # [1] 0.6931472  parameter for logistic regression
    levels=c(1,0)
    X=c(rep(1,30),rep(0,70))  # fixed
    z1 = beta0 + (betax *X)      # linear combination 
    Py1 = 1/(1+exp(-z1))   # inv-logit function
    y1 <- as.numeric(rbinom(100,1,Py1))
    betaF2[i] <- glm(y1 ~ X, family="binomial"(link=logit))$coeff[2]
  }
  betamean10000[j] = mean(betaF2)
}


# to observe if the results (I would have obtained if I used M=500,000 instead of 1000)
# are consistent with those seen in the study described
# in the paper, I simply replicate the means of MLE 500 times to obtain
# a vector betamean100_2 of length 500,000
# "if I had this quantity of this quality" principle

betamean100_2<-rep(betamean100,500) 
betamean500_2<-rep(betamean500,500)
betamean1000_2<-rep(betamean1000,500)
betamean2500_2<-rep(betamean2500,500)
betamean5000_2<-rep(betamean5000,500)
betamean10000_2<-rep(betamean10000,500)

# function to compute the mean percent bias for MCEclt
percentbiasfunction<-function(x, target=NULL)  {
  value=sd((x - target) / target * 100) / sqrt(length(x))
  return(value)
}

# MCE clt (M= 500,000)
percentbiasfunction(betamean100_2, betax)
# [1] 0.009638844
percentbiasfunction(betamean500_2, betax)
# [1] 0.004541754
percentbiasfunction(betamean1000_2, betax)
# [1] 0.003280034
percentbiasfunction(betamean2500_2, betax)
# [1] 0.001922758
percentbiasfunction(betamean5000_2, betax)
# [1] 0.00134847
percentbiasfunction(betamean10000_2, betax)
# [1] 0.0009288466

# I also compute mean percent bias mce clt for the data I obtain, namely 
# betamean100 with M=1,000. I anticipate larger mce clt.

# MCE clt (M= 1,000)
percentbiasfunction(betamean100, betax)
# [1] 0.2156387
percentbiasfunction(betamean500, betax)
# [1] 0.1016074
percentbiasfunction(betamean1000, betax)
# [1] 0.07338041
percentbiasfunction(betamean2500, betax)
# [1] 0.04301565
percentbiasfunction(betamean5000, betax)
# [1] 0.03016776
percentbiasfunction(betamean10000, betax)
# [1] 0.02078001

## parametric bootstrap

# 1. for Percent Bias
# selected functions available in the package MCE
# remove from R (archived), available here https://cran.r-project.org/web/packages/MCE/index.html
calcPB <-
  function(data, index, truth) (mean(data[index]) - truth) / truth * 100

calcSE <-
  function(data, index) sd(data[index])


mceBoot <-
  function(data, B, type="", truth=NULL)
  {
    ##
    R <- switch(type,
                "mean"=length(data),
                "PB"=length(data),
                "SE"=length(data),
                "RE"=nrow(data))
    
    ##
    subsets <- matrix(sample(1:R, R*B, replace=TRUE), ncol=B)
    bootRes <- rep(NA, B)
    for(b in 1:B)
    {
      if(type == "mean") bootRes[b] <- calcEX(data, subsets[,b])
      if(type == "PB")   bootRes[b] <- calcPB(data, subsets[,b], truth=truth)
      if(type == "SE")   bootRes[b] <- calcSE(data, subsets[,b])
      if(type == "RE")   bootRes[b] <- calcRE(data, subsets[,b])
    }
    value <- sd(bootRes)
    return(value)
  }

## for  M=500,000
# R=100
set.seed(3)
mceBoot(betamean100_2, B=100,type="PB", truth=betax)
# [1] 0.009481879
mceBoot(betamean100_2, B=200,type="PB", truth=betax)
# [1] 0.01030212
mceBoot(betamean100_2, B=500,type="PB", truth=betax)
# [1] 0.009391162

# R=500
set.seed(3)
mceBoot(betamean500_2, B=100,type="PB", truth=betax)
# [1] 0.003722718
mceBoot(betamean500_2, B=200,type="PB", truth=betax)
# [1] 0.004728553
mceBoot(betamean500_2, B=500,type="PB", truth=betax)
# [1] 0.00460928

# R=1000
set.seed(3)
mceBoot(betamean1000_2, B=100,type="PB", truth=betax)
# [1] 0.003515548
mceBoot(betamean1000_2, B=200,type="PB", truth=betax)
# [1] 0.00315945
mceBoot(betamean1000_2, B=500,type="PB", truth=betax)
# [1] 0.003183787

# R=2500
set.seed(3)
mceBoot(betamean2500_2, B=100,type="PB", truth=betax)
# [1] 0.002175779
mceBoot(betamean2500_2, B=200,type="PB", truth=betax)
# [1] 0.001860283
mceBoot(betamean2500_2, B=500,type="PB", truth=betax)
# [1] 0.002044577

# R=5000
set.seed(3)
mceBoot(betamean5000_2, B=100,type="PB", truth=betax)
# [1] 0.001351824
mceBoot(betamean5000_2, B=200,type="PB", truth=betax)
# [1] 0.001233471
mceBoot(betamean5000_2, B=500,type="PB", truth=betax)
# [1] 0.001323096

# R=10000
set.seed(3)
mceBoot(betamean10000_2, B=100,type="PB", truth=betax)
# [1] 0.0009799865
mceBoot(betamean10000_2, B=200,type="PB", truth=betax)
# [1] 0.0008644101
mceBoot(betamean10000_2, B=500,type="PB", truth=betax)
# [1] 0.0009861846


## for  M=1,000
# R=100
set.seed(3)
mceBoot(betamean100, B=100,type="PB", truth=betax)
# [1] 0.2210012
mceBoot(betamean100, B=200,type="PB", truth=betax)
# [1] 0.2268075
mceBoot(betamean100, B=500,type="PB", truth=betax)
# [1] 0.2128067

# R=500
set.seed(3)
mceBoot(betamean500, B=100,type="PB", truth=betax)
# [1] 0.09582288
mceBoot(betamean500, B=200,type="PB", truth=betax)
# [1] 0.1027374
mceBoot(betamean500, B=500,type="PB", truth=betax)
# [1] 0.1037522

# R=1000
set.seed(3)
mceBoot(betamean1000, B=100,type="PB", truth=betax)
# [1] 0.07237398
mceBoot(betamean1000, B=200,type="PB", truth=betax)
# [1] 0.07823361
mceBoot(betamean1000, B=500,type="PB", truth=betax)
# [1] 0.07268609

# R=2500
set.seed(3)
mceBoot(betamean2500, B=100,type="PB", truth=betax)
# [1] 0.04160736
mceBoot(betamean2500, B=200,type="PB", truth=betax)
# [1] 0.03973145
mceBoot(betamean2500, B=500,type="PB", truth=betax)
# [1] 0.04324723

# R=5000
set.seed(3)
mceBoot(betamean5000, B=100,type="PB", truth=betax)
# [1] 0.03083821
mceBoot(betamean5000, B=200,type="PB", truth=betax)
# [1] 0.03073943
mceBoot(betamean5000, B=500,type="PB", truth=betax)
# [1] 0.02968757

# R=10000
set.seed(3)
mceBoot(betamean10000, B=100,type="PB", truth=betax)
# [1] 0.02089425
mceBoot(betamean10000, B=200,type="PB", truth=betax)
# [1] 0.02239513
mceBoot(betamean10000, B=500,type="PB", truth=betax)
# [1] 0.02092768


# 2. for Monte Carlo Error MCE

phi_hat_b_r100_2 <- rep(phi_hat_b_r100, 500)
phi_hat_b_r500_2 <-rep(phi_hat_b_r500,500)
phi_hat_b_r1000_2 <-rep(phi_hat_b_r1000,500)
phi_hat_b_r2500_2 <-rep(phi_hat_b_r2500,500)
phi_hat_b_r5000_2<-rep(phi_hat_b_r5000,500)
phi_hat_b_r10000_2<-rep(phi_hat_b_r10000,500)

## M=500,000

# R=100
set.seed(3)
mceBoot(phi_hat_b_r100_2, B=100,type="SE", truth=NULL)
# [1] 0.006674453
mceBoot(phi_hat_b_r100_2, B=200,type="SE", truth=NULL)
# [1] 0.007113146
mceBoot(phi_hat_b_r100_2, B=500,type="SE", truth=NULL)
# [1] 0.006994765

# R=500
set.seed(3)
mceBoot(phi_hat_b_r500_2, B=100,type="SE", truth=NULL)
# [1] 0.002912512
mceBoot(phi_hat_b_r500_2, B=200,type="SE", truth=NULL)
# [1] 0.003158602
mceBoot(phi_hat_b_r500_2, B=500,type="SE", truth=NULL)
# [1] 0.003046988

# R=1000
set.seed(3)
mceBoot(phi_hat_b_r1000_2, B=100,type="SE", truth=NULL)
# [1] 0.002178902
mceBoot(phi_hat_b_r1000_2, B=200,type="SE", truth=NULL)
# [1] 0.002248723
mceBoot(phi_hat_b_r1000_2, B=500,type="SE", truth=NULL)
# [1] 0.002345778

# R=2500
set.seed(3)
mceBoot(phi_hat_b_r2500_2, B=100,type="SE", truth=NULL)
# [1] 0.001209567
mceBoot(phi_hat_b_r2500_2, B=200,type="SE", truth=NULL)
# [1] 0.001475144
mceBoot(phi_hat_b_r2500_2, B=500,type="SE", truth=NULL)
# [1] 0.00130233

# R=5000
set.seed(3)
mceBoot(phi_hat_b_r5000_2, B=100,type="SE", truth=NULL)
# [1] 0.0009052758
mceBoot(phi_hat_b_r5000_2, B=200,type="SE", truth=NULL)
# [1] 0.0008801176
mceBoot(phi_hat_b_r5000_2, B=500,type="SE", truth=NULL)
# [1] 0.0008987816

# R=10000
set.seed(3)
mceBoot(phi_hat_b_r10000_2, B=100,type="SE", truth=NULL)
# [1] 0.000614785
mceBoot(phi_hat_b_r10000_2, B=200,type="SE", truth=NULL)
# [1] 0.0006543925
mceBoot(phi_hat_b_r10000_2, B=500,type="SE", truth=NULL)
# [1] 0.0006535844


## M=1,000

# R=100
set.seed(3)
mceBoot(phi_hat_b_r100, B=100,type="SE", truth=NULL)
# [1] 0.1600673
mceBoot(phi_hat_b_r100, B=200,type="SE", truth=NULL)
# [1] 0.1547736
mceBoot(phi_hat_b_r100, B=500,type="SE", truth=NULL)
# [1] 0.1599613

# R=500
set.seed(3)
mceBoot(phi_hat_b_r500, B=100,type="SE", truth=NULL)
# [1] 0.06775524
mceBoot(phi_hat_b_r500, B=200,type="SE", truth=NULL)
# [1] 0.06603143
mceBoot(phi_hat_b_r500, B=500,type="SE", truth=NULL)
# [1] 0.06993108

# R=1000
set.seed(3)
mceBoot(phi_hat_b_r1000, B=100,type="SE", truth=NULL)
# [1] 0.0481949
mceBoot(phi_hat_b_r1000, B=200,type="SE", truth=NULL)
# [1] 0.05243203
mceBoot(phi_hat_b_r1000, B=500,type="SE", truth=NULL)
# [1] 0.05009632

# R=2500
set.seed(3)
mceBoot(phi_hat_b_r2500, B=100,type="SE", truth=NULL)
# [1] 0.02907816
mceBoot(phi_hat_b_r2500, B=200,type="SE", truth=NULL)
# [1] 0.02785746
mceBoot(phi_hat_b_r2500, B=500,type="SE", truth=NULL)
# [1] 0.02928191

# R=5000
set.seed(3)
mceBoot(phi_hat_b_r5000, B=100,type="SE", truth=NULL)
# [1] 0.02105426
mceBoot(phi_hat_b_r5000, B=200,type="SE", truth=NULL)
# [1] 0.01911916
mceBoot(phi_hat_b_r5000, B=500,type="SE", truth=NULL)
# [1] 0.02010327

# R=10000
set.seed(3)
mceBoot(phi_hat_b_r10000, B=100,type="SE", truth=NULL)
# [1] 0.01462418
mceBoot(phi_hat_b_r10000, B=200,type="SE", truth=NULL)
# [1] 0.01275701
mceBoot(phi_hat_b_r10000, B=500,type="SE", truth=NULL)
# [1] 0.0142794

# With regards to the results obtained with the two methods, i.e.central limit theorem
# and bootstrap, we observe a reduction of MCE of more than 90% in using R=10000 
# compared to R=100, whatever the number of simulations M.

(0.01462418 - 0.1600673) / 0.1600673   # M=1,000
# [1] -0.9086373
(0.000614785 - 0.006674453) / 0.006674453  # M=500,000
# [1] -0.9078898

# In this respect, my results are consistent with those from the study.

#----
# end
#----
