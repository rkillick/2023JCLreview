## Sea Ice example
## Review Paper for Journal of Climate 
## Author: Xueheng Shi
## Date: 09/16/2022


rm(list = ls())

library(latex2exp)
library(wbs)

#################### Load data Sea Ice###################
raw = read.csv(file.choose(), header=T)
View(raw)


##Plot the data
sea.ice = ts(raw[ ,2], frequency = 1, start=1979, end=2021) 

##Binary Segmentation, also called sbs
xt = raw[,2]
N = length(xt)
bs = wbs::sbs(xt)
bscpt = wbs::changepoints(bs)
mBS = bscpt$no.cpt.th #number of CPT
tBS = as.vector(unlist(bscpt$cpt.th)) #CPT location
#t.BS = ifelse(mBS==0, "NA", toString(tBS) )
tBS ##changepoint times
#16 23 37
tBS+1978 ##changepoint years
#2001 1994 2015


pdf("SeaIce-BS-AR1.pdf", width = 10, height = 5)
par(oma=c(0,2,0,0))
plot(sea.ice, xlab="Year", ylab=TeX("Sea Ice $(\\times 10^6\\, km^2)$"))
#segments(1900, mu1.hat, 1987, mu1.hat, col= "red", lwd=3)
#segments(1988, mu2.hat, 2021, mu2.hat, col= "red", lwd=3)
segments(1979, mean(sea.ice[1:15]), 1994, mean(sea.ice[1:15]),  col= "red", lwd=3)
segments(1994, mean(sea.ice[16:22]), 2001, mean(sea.ice[16:22]),  col= "red", lwd=3)
segments(2001, mean(sea.ice[23:36]), 2015,mean(sea.ice[23:36]),  col= "red", lwd=3)
segments(2015, mean(sea.ice[37:43]), 2021,mean(sea.ice[37:43]),  col= "red", lwd=3)
abline(v=c(2001, 1994, 2015), lty="dashed", col="red", lwd=2)
dev.off()



##Fit Mean shifts + AR1 errors
source("GASearch.R")
##Run BIC-GA
BICGA = GA::ga(type="binary", fitness = BinSearch.BIC,
               nBits = N, maxiter = 10000, run = 1000,
               popSize = 50, monitor = T)

BICsol = (BICGA@solution)[1, ] ##can be 2 sols since X[1] is free, take 1st sol
BICsol[1] = 0 
mBIC = sum(BICsol)

tBIC = as.vector( which(BICsol == 1) )
tBIC
#17 28 38 39
tBIC+1978
#1995 2006 2016 2017

pdf("SeaIce-BIC-AR1.pdf", width = 10, height = 5)
par(oma=c(0,2,0,0))
plot(sea.ice, xlab="Year", ylab=TeX("Sea Ice $(\\times 10^6\\, km^2)$ "))
segments(1979, mean(sea.ice[1:16]), 1995, mean(sea.ice[1:16]),  col= "red", lwd=3)
segments(1995, mean(sea.ice[18:27]), 2006, mean(sea.ice[18:27]),  col= "red", lwd=3)
segments(2006, mean(sea.ice[28:37]), 2016,mean(sea.ice[28:37]),  col= "red", lwd=3)
segments(2016, mean(sea.ice[38:38]), 2017,mean(sea.ice[38:38]),  col= "red", lwd=3)
segments(2017, mean(sea.ice[39:43]), 2021,mean(sea.ice[39:43]),  col= "red", lwd=3)
abline(v=c(1995, 2006, 2016, 2017), lty="dashed", col="red", lwd=2)
dev.off()


##Fit Mean shifts + Fixed Trend + AR1 errors
source("MeanCPT_Trend_AR1Errors.R")

BIC.out = GA::ga(type="binary", fitness = MeanChange.BIC,
                 nBits = 14, maxiter = 10000, run=2000,
                 popSize = 100, monitor = T)
summary(BIC.out)
BIC.trend = (BIC.out@solution)[1,] ##can be 2 sols
BIC.trend[1]=0
mBIC2 = sum(BIC.trend)
tBIC2 = as.vector(which(BIC.trend == 1))
tBIC2*3 ##the search uses a mini-spacing =3
##Result is Zero Changepoint

## Estimate the trend for plotting
## Zero CPT so the design matrix D is easy to compute
Y = matrix(xt, nrow=N, byrow=T) #Y: Dependent Variable
D = matrix(c(rep(1,N), seq(1,N)), nrow=N, byrow=F) #design matrix 

beta.hat = solve((t(D)%*%D))%*%t(D)%*%Y ##Regression Coefficients
#beta.hat
#[1,] 11.59910299
#[2,] -0.05316219

Rt = Y - D%*%beta.hat  ##AR(1) residuals
phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
#[1] 0.05547652

##Endpoints of the segments, for plotting
x1=1979
y1=11.59910299
x2=2021
y2=11.59910299-0.05316219*(2021-1979)

pdf("SeaIce-BIC-Trend+AR1.pdf", width = 10, height = 5)
par(oma=c(0,2,0,0))
plot(sea.ice, xlab="Year", ylab=TeX("Sea Ice $(\\times 10^6\\, km^2)$"))
#segments(1900, mu1.hat, 1987, mu1.hat, col= "red", lwd=3)
segments(x1, y1, x2, y2, lty="dashed", col= "red", lwd=2)
abline(a=11.59910299,b=-0.05316219, lty="dashed", col="red", lwd=2)
dev.off()
