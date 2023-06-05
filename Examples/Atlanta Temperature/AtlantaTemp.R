## GA Search and Binary Search
## Review Paper for Journal of Climate 
## Author: Xueheng Shi
## Date: 07/08/2022

rm(list = ls())

library(latex2exp)
library(wbs)

source("GASearch.R")

##Differenced Series
##based on Shi's Autocovariance Paper(2022) and using theta=-1
##Yule-Walker Estimator at lag1, i.e., gamma(1)/gamma(0)
phi.yw = function(x){
  ## x is input time series
  acv.h = acf(x, lag.max = 1, type = "covariance", plot = F)$acf
  acv.h[2]/acv.h[1] 
}

phi.mom = function(x){
  y = diff(x)
  return(2*phi.yw(y)+1)
}


##Atlanta Example
ATL.dat = read.csv(file.choose(), head = T) #Load Atlanta_YrAvg_C3.csv
#View(ATL.dat)
ATL.Temp = ts(ATL.dat[ ,2], frequency = 1, start=1879, end=2012) 

xt = ATL.dat[ ,2]
N = length(xt)


##Run BIC-GA
BICGA = GA::ga(type="binary", fitness = BinSearch.BIC,
               nBits = N, maxiter = 20000, run = 4000,
               popSize = 200, monitor = T)

BICsol = (BICGA@solution)[1, ] ##can be 2 sols since X[1] is free, take 1st sol
BICsol[1] = 0  ##1st data point is not a changepoint in our definition
mBIC = sum(BICsol) #Number of detected changepoints

tBIC = as.vector( which(BICsol == 1) )
tBIC ##changepoint times (indices)
#43  82 106
tBIC+1878 ## changepoint times (years)
##[1] 1921 1960 1984

##Plot the model
pdf("Atlanta-BIC.pdf", width = 10, height = 5)
par(oma=c(0,2,0,0))
plot(ATL.Temp, xlab="Year", ylab=TeX("Average Temperature (^oC)"), main="Atlanta")
segments(1879, mean(xt[1:42]), 1921, mean(xt[1:42]),  col= "red", lwd=3)
segments(1921, mean(xt[43:81]), 1960, mean(xt[43:81]),  col= "red", lwd=3)
segments(1960, mean(xt[82:105]), 1984,mean(xt[82:105]),  col= "red", lwd=3)
segments(1984, mean(xt[106:134]), 2012,mean(xt[106:134]),  col= "red", lwd=3)
abline(v=c(1921, 1960, 1984), lty="dashed", col="red", lwd=2) ##GA-BIC
dev.off()



##Binary Segmentation 
##Changepoint detection on the AR(1) residuals
phi.hat = phi.mom(xt)
Zt = xt[-1] - phi.hat*xt[-N] ##AR(1) residuals

##Binary Segmentation, also called sbs
bs = wbs::sbs(Zt)
bscpt = wbs::changepoints(bs)
mBS = bscpt$no.cpt.th #number of changepoints
tBS = as.vector(unlist(bscpt$cpt.th)) 
tBS ##changepoint times(indices)
#104
tBS+1878 ##changepoint years
#1982

##Plot the model
pdf("Atlanta-BS.pdf", width = 10, height = 5)
par(oma=c(0,2,0,0))
plot(ATL.Temp, xlab="Year", ylab=TeX("Average Temperature (^oC)"), main="Atlanta")
segments(1879, mean(xt[1:103]), 1982, mean(xt[1:103]),  col= "red", lwd=3)
segments(1982, mean(xt[104:134]), 2012,mean(xt[104:134]),  col= "red", lwd=3)
abline(v=c(1982), lty="dashed", col="red", lwd=2) ##GA-BIC
dev.off()
