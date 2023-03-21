#------------------------------------------------
# Bootstrap-based methods for Quantile regression
#------------------------------------------------

library("quantreg")

# 1. introductory example
n = 30
set.seed(1986)
x = rnorm(n)
beta0 = 5
beta1 = 2
y = beta0 + beta1*x + rnorm(n)
x = c(x, 3)
y = c(y, 4)

# Linear regression model (OLS)
model.1 = lm(y~x)

plot(x,y, pch = 18)
range = seq(-(1+min(x)),(1+max(x)),0.01)
lines(x,model.1$fitted.values, col = "black", lwd = 2)

# Quantile regression model
model.2 = rq(y~x)
lines(x,model.2$fitted.values, col = "firebrick2", lwd = 2)
title("Linear regression vs Quantile regression")
legend("bottomright", c("Linear regression (OLS) estimator","Quantile regression estimator"), lwd=c(2.5,2.5), 
        col =  c("black", "firebrick2"), cex = 0.8)

# 2. code example (1) (2)
  
# bootstrap
start.time <- Sys.time()

set.seed(1986)
tau = 0.5
MC_sim = 1000
MC_length_CG = MC_length_BN = matrix(rep(0, MC_sim*3),MC_sim, 3)
MC_coverage_CG = MC_coverage_BN <- matrix(rep(0, MC_sim*3),MC_sim, 3)

for(j in 1:MC_sim){
  n = 50
  b0 = b1 = b2 = 1
  e = rt(n,3)
  x1 = rlnorm(n)
  x2 <- numeric(n)
  x2[1:(0.8*n)] = 1
  x2[((0.8*n)+1):n] = 0
  
  y <- numeric(n)
  y = b0 + b1*x1 + b2*x2 + 3^(-1/2)*(2+(1+(x1-8)^2 + x2)/10)*e
  M1 = rq(y~x1+x2)
  
  boot_sim = 1000
  coef_boot_CG = coef_boot_BN <- matrix(rep(0, boot_sim*3),boot_sim, 3)
  w = w_CG = w_BN <- numeric(n)
  
  #Bootstrap
  
  #CG
  B_CG = boot.rq(cbind(rep(1,length(x1)), x1,x2), y, R = 1000, bsmethod = "wild")
  coef_boot_CG = B_CG$B
  
  #BN
  for(i in 1:boot_sim){
    u = runif(n)
    w_BN[u<=0.5] = 1
    w_BN[u>0.5] = -1
    # Correction for residuals in finite sample
    x = cbind(rep(1,length(x1)), x1,x2)
    r = M1$residuals
    f0 <- akj(r, z = 0)$dens
    r <- r + hat(x) * (tau - I(r < 0))/f0
    
    y_boot_BN = M1$fitted.values + w_BN * abs(r)
    coef_boot_BN[i,] = rq(y_boot_BN~x1+x2)$coefficients
  }
  
  # CG
  MC_quant = quantile(coef_boot_CG[,1], c(0.05,0.95))
  MC_length_CG[j,1] = MC_quant[2] - MC_quant[1]
  MC_coverage_CG[j,1] = (MC_quant[2]>=b0 && MC_quant[1]<=b0)
  
  MC_quant = quantile(coef_boot_CG[,2], c(0.05,0.95))
  MC_length_CG[j,2] = MC_quant[2] - MC_quant[1]
  MC_coverage_CG[j,2] = (MC_quant[2]>=b1 && MC_quant[1]<=b1)
  
  MC_quant = quantile(coef_boot_CG[,3], c(0.05,0.95))
  MC_length_CG[j,3] = MC_quant[2] - MC_quant[1]
  MC_coverage_CG[j,3] = (MC_quant[2]>=b2 && MC_quant[1]<=b2)
  
  # BN
  MC_quant = quantile(coef_boot_BN[,1], c(0.05,0.95))
  MC_length_BN[j,1] = MC_quant[2] - MC_quant[1]
  MC_coverage_BN[j,1] = (MC_quant[2]>=b0 && MC_quant[1]<=b0)
  
  MC_quant = quantile(coef_boot_BN[,2], c(0.05,0.95))
  MC_length_BN[j,2] = MC_quant[2] - MC_quant[1]
  MC_coverage_BN[j,2] = (MC_quant[2]>=b1 && MC_quant[1]<=b1)
  
  MC_quant = quantile(coef_boot_BN[,3], c(0.05,0.95))
  MC_length_BN[j,3] = MC_quant[2] - MC_quant[1]
  MC_coverage_BN[j,3] = (MC_quant[2]>=b2 && MC_quant[1]<=b2)
}
AGG_MC_length_CG = colMeans(MC_length_CG)
AGG_MC_length_BN = colMeans(MC_length_BN)

AGG_MC_coverage_CG = colMeans(MC_coverage_CG)
AGG_MC_coverage_BN = colMeans(MC_coverage_BN)

Results = c("CG", round(AGG_MC_coverage_CG[1]*100,1), round(AGG_MC_length_CG[1],1), round(sd(MC_length_CG[,1])/sqrt(n),2),
            round(AGG_MC_coverage_CG[2]*100,1), round(AGG_MC_length_CG[2],1), round(sd(MC_length_CG[,2])/sqrt(n),2),
            round(AGG_MC_coverage_CG[3]*100,1), round(AGG_MC_length_CG[3],1), round(sd(MC_length_CG[,3])/sqrt(n),2)
)

Results = rbind(Results,c("BN", round(AGG_MC_coverage_BN[1]*100,1), round(AGG_MC_length_BN[1],1), round(sd(MC_length_BN[,1])/sqrt(n),2),
                          round(AGG_MC_coverage_BN[2]*100,1), round(AGG_MC_length_BN[2],1), round(sd(MC_length_BN[,2])/sqrt(n),2),
                          round(AGG_MC_coverage_BN[3]*100,1), round(AGG_MC_length_BN[3],1), round(sd(MC_length_BN[,3])/sqrt(n),2)
))

colnames(Results) = c("Method","b0 Cover","b0 length", "b0 sd","b1 Cover","b1 length", "b1 sd","b2 Cover","b2 length", "b2 sd") 
Results

end.time <- Sys.time()
time <- end.time - start.time
time
  
# 3. practical example: low birthweight

data(birthwt, package = 'MASS')
data <- birthwt
head(data)
#    low age lwt race smoke ptl ht ui ftv  bwt
# 85   0  19 182    2     0   0  0  1   0 2523
# 86   0  33 155    3     0   0  0  0   3 2551
# 87   0  20 105    1     1   0  0  0   1 2557
# 88   0  21 108    1     1   0  0  1   2 2594
# 89   0  18 107    1     1   0  0  1   0 2600
# 91   0  21 124    3     0   0  0  0   0 2622
dim(data)
# [1] 189  10
class(data)
# [1] "data.frame"

quantile.model <- rq(bwt ~ age + race + smoke, 
                     tau = seq(0.05, 0.95, by = 0.05),
                     data = data)
summary(quantile.model)
# tau: [1] 0.5
# Coefficients:
#           coefficients     lower bd   upper bd  
# (Intercept) 3248.20000   2620.25015 4011.99351
# age           17.60000    -19.56983   33.77936
# race        -283.20000   -362.99101  -42.80247
# smoke       -512.80000   -701.82717 -241.68587

# wild bootstrap
quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                 tau = 0.5, R = 1000, bsmethod = "wild")$B[,1], c(0.025, 0.975))
#       2.5%    97.5% 
#   2509.913 4006.367 

# wild bootstrap loops
attach(data)
set.seed(1986)

# compute vector of lower quantiles
birthweight.CG.lower <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
birthweight.CG.lower[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,1], 0.025)[[1]]
}
# compute vector of lower quantiles
birthweight.CG.upper <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.upper[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,1], 0.975)[[1]]
}
# data frame for beta0
dfbeta0 <- data.frame(quantile = seq(0.05, 0.95, by = 0.05),
                      beta0 = quantile.model$coefficients[1,],
                      lower = birthweight.CG.lower,
                      upper = birthweight.CG.upper)

# compute vector of lower quantiles
birthweight.CG.lower <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.lower[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,2], 0.025)[[1]]
}
# compute vector of lower quantiles
birthweight.CG.upper <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.upper[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,2], 0.975)[[1]]
}
# data frame for beta1
dfbeta1 <- data.frame(quantile = seq(0.05, 0.95, by = 0.05),
                      beta1 = quantile.model$coefficients[2,],
                      lower1 = birthweight.CG.lower,
                      upper1 = birthweight.CG.upper)

# compute vector of lower quantiles
birthweight.CG.lower <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.lower[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,3], 0.025)[[1]]
}
# compute vector of lower quantiles
birthweight.CG.upper <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.upper[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,3], 0.975)[[1]]
}
# data frame for beta2
dfbeta2 <- data.frame(quantile = seq(0.05, 0.95, by = 0.05),
                      beta2 = quantile.model$coefficients[3,],
                      lower2 = birthweight.CG.lower,
                      upper2 = birthweight.CG.upper)

# compute vector of lower quantiles
birthweight.CG.lower <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.lower[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,4], 0.025)[[1]]
}
# compute vector of lower quantiles
birthweight.CG.upper <- numeric(length(seq(0.05, 0.95, by = 0.05)))
for(i in 1:length(seq(0.05, 0.95, by = 0.05))) {
  birthweight.CG.upper[i] = quantile(boot.rq(cbind(rep(1,length(bwt)),age,race,smoke), bwt,
                                             tau = i/20, R = 1000, bsmethod = "wild")$B[,4], 0.975)[[1]]
}
# data frame for beta3
dfbeta3 <- data.frame(quantile = seq(0.05, 0.95, by = 0.05),
                      beta3 = quantile.model$coefficients[4,],
                      lower3 = birthweight.CG.lower,
                      upper3 = birthweight.CG.upper)
# plot
library(ggplot2)
library(gridExtra)
library(grid)
p1<-ggplot(data=dfbeta0, aes(x=quantile, y=beta0)) + 
  geom_point() + 
  geom_line() +
  ylab('intercept') +
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, col = 'firebrick2', alpha=0.1) +
  theme_minimal()
p2<-ggplot(data=dfbeta1, aes(x=quantile, y=beta1)) + 
  geom_point() + 
  geom_line() +
  ylab('age') +
  geom_ribbon(aes(ymin=lower1, ymax=upper1), linetype=2, col = 'firebrick2', alpha=0.1) +
  theme_minimal()
p3<-ggplot(data=dfbeta2, aes(x=quantile, y=beta2)) + 
  geom_point() + 
  geom_line() +
  ylab('race') +
  geom_ribbon(aes(ymin=lower2, ymax=upper2), linetype=2, col = 'firebrick2', alpha=0.1) +
  theme_minimal()
p4<-ggplot(data=dfbeta3, aes(x=quantile, y=beta3)) + 
  geom_point() + 
  geom_line() +
  ylab('smoke') +
  geom_ribbon(aes(ymin=lower3, ymax=upper3), linetype=2, col = 'firebrick2', alpha=0.1) +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, nrow = 2,
             top = textGrob("Wild bootstrap CI for parameter estimates",gp=gpar(fontsize=20,font=1)))


#######
# end #
#######