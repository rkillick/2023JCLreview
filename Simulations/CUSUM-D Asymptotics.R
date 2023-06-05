## Compute the quantitles CUSUM-D Distribution
## Used to perform CUSUM-D test (of a single changepoint) 
## Under the assumption that model has a fixed trend and AR(1) errors
## Author: Xueheng Shi
## Based on "Changepoint Detection in Climate Time Series with Long-Term Trends
## by COLIN GALLAGHER AND ROBERT LUND
## Date: 03/21/2021

library(e1071)

n = 10000
t = seq(0,1,length=n)
Super.Bhz = NULL

set.seed(123)

for(k in 1:10^6){ 
  ##Brownian Motion and Brownian Bridge
  bm = c(0, cumsum(rnorm(n - 1,0,sqrt(1/n))))
  BB = bm - t*rep(bm[length(bm)], length.out = length(bm))
  #BBz=rbridge(end = 1, frequency = 10000)
  inteBB = sum(BB)/n
  
  #Super.Bz = max( abs(BB-6*t*(1-t)*inteBBZ) )
  Super.Bz = max( abs(BB-6*t*(1-t)*inteBB) )
  Super.Bhz = rbind(Super.Bhz, Super.Bz)
  
  #print(k)
}

##Based on 10^6 replications
quantile(Super.Bhz, probs = c(0.9,0.95,0.975,0.99))
#90%       95%     97.5%       99% 
#  0.8321218 0.9021685 0.9656287 1.0426011 

pvalue = sum(Super.Bhz>0.928964)/length(Super.Bhz)
pvalue
##[1] 0.0377
