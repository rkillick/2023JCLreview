## Central England Temperature
## Review Paper for Journal of Climate 
## Author: Xueheng Shi
## Date: 05/02/2022
## Update: 07/01/2022

rm(list = ls())

library(latex2exp)

####################Functions########################
## compute phi for AR(1) series
phi = function(x){
  ## x is input time series
  n = length(x)
  sum(x[-1]*x[-n])/sum(x^2)
}



#################### Data Analysis ###################
CET = read.csv(file.choose(),header=T) ##Load CentralEnglandT1659-2020_yearly.csv

##use the data 1900-2020
CET120 = CET[242:362, ]
Temp = as.vector( CET120[,2] )


#################### Data Analysis ###################


#CentralEngland.Yr.Temp = read.csv(file='CentralEnglandT1659-2020_yearly.csv', header=T)
#View(CentralEngland.Yr.Temp)

Yr.Temp = ts(CET120[ ,2], frequency = 1, start=1900, end=2020) 
Xt = Temp
n = length(Xt)

##CUSUM Test on Xt, assuming IID errors
abs.CUSUMx = abs(cumsum(Xt)-sum(Xt)*(1:n)/n)/sqrt(n)
k=which.max(abs.CUSUMx) ##Changepoint time
max.CUSUMx =  max(abs.CUSUMx) ##Uncropped CUSUMx
k
##Under H0
sigma2.hat = sum((Xt-mean(Xt))^2)/n
Cusum.Score = max.CUSUMx/sqrt(sigma2.hat)
Cusum.Score>1.358099 ##CUSUMx Test result, if TRUE then reject H0:no changepoint

##SCUSUM developed for IID by Kirch
SCUSUM.K = (1/n)*sum(abs.CUSUMx^2)/sigma2.hat
SCUSUM.K 


##CUSUMz Test for AR1 errors (based on residuals of AR1 series)
phi.hat = phi(Xt)
muhat = mean(Xt)
Zt = c(Xt[1], (Xt[-1] - muhat) - phi.hat*(Xt[-n] - muhat) )
sigma.sq = var(Zt)
absCUSUMz = abs(cumsum(Zt)-sum(Zt)*(1:n)/n)/sqrt(n)
maxCUSUMz =  max(absCUSUMz)  ##CUSUMz Statistic
maxCUSUMz/sqrt(sigma.sq)>1.358099
k2=which.max(absCUSUMz)
k2

##SCUSUMz Test (based on residuals of AR1 series)
(1/n)*sum(absCUSUMz^2)/sigma.sq
##0.1798982




## Single Mean Shift Test, assuming fixed Trend and AR1 errors
tt = 1:n
beta.hat = 12*sum( Xt*(tt-mean(tt)) )/(n*(n+1)*(n-1))
##[1] 0.009229982
mu.hat = mean(Xt)-beta.hat*mean(tt)
et = Xt-(mu.hat+beta.hat*tt) ##autocorrelated residuals

##Compute under H0
phi2 = phi(et)
#[1] 0.1936949
muhat = mean(et)
Zt = c(et[1], (et[-1] - muhat) - phi2*(et[-n] - muhat) )
sigma2 = sum( (Zt-mean(Zt))^2 )/(n-2)

absDk = abs(cumsum(Zt))/sqrt(n)
maxDk =  max(absDk)  ##eqn 3.3-3.4 of Lund (2012)
maxDk/sqrt(sigma2)>0.906 ##critical value from CUSUM-D Asymptotic
which.max(absDk)
##88
#> maxDk/sqrt(sigma2)
#[1] 0.928964


acf(Xt, type = "correlation", plot=F)
##0.425



##Plots of the single changepoint tests
##Plots for model fitting

## IID Fitting Result
##starting year = 1990, which index is 1
mu1.hat = mean(Yr.Temp[1:88])
mu2.hat = mean(Yr.Temp[89:121])

pdf("CET1900-2020IID-SCPT.pdf", width = 10, height = 5)
plot(Yr.Temp, xlab="Year", ylab=TeX("Average Temperature (^oC)"))
segments(1900, mu1.hat, 1987, mu1.hat, col= "red", lwd=3)
segments(1988, mu2.hat, 2021, mu2.hat, col= "red", lwd=3)
abline(v=1988, lty="dashed", col="red", lwd=2)
dev.off()


## Fixed Slope+AR1 Errors Model
## Estimate the coefficients for segments
MeanCptTrendEstimate = function(Xt, tau){
  
  N = length(Xt) #length of Xt
  loc.ind = rep(0, N)
  loc.ind[tau] = 1
  m = sum(loc.ind) #Number of CPTs in configuration
  
  Y = matrix(Xt, nrow=N, byrow=T) #Y: Dependent Variable
  D = matrix(c(rep(1,N), seq(1,N)), nrow=N, byrow=F) #design mat for zero cpt
  
  if(m == 0){ #Zero CPT
    beta.hat = solve((t(D)%*%D))%*%t(D)%*%Y ##Regression Coefficients
    Rt = Y - D%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
  }
  else{
    #tau.vec = loc.ind*(1:N) #CPT index
    tau = which(loc.ind == 1) #CPT locations
    tau.ext = c(1, tau, (N+1)) #CPT boundary 1 and N+1
    #seg.len = diff(tau.ext)  #length of each segment
    
    Dm = matrix(0, nrow=N, ncol=m+1, byrow=T) #Design mat for mean shifts
    for(i in seq(1,m+1) ){
      #print(i)
      Dm[seq(tau.ext[i], tau.ext[i+1]-1), i]=1
    }
    M = cbind(Dm,seq(1,N))
    beta.hat = solve((t(M)%*%M))%*%t(M)%*%Y ##Regression Coefs
    Rt = Y - M%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
    
  }
  return( c(phi.hat, beta.hat) )
}

MeanCptTrendEstimate(Xt, c(89))
##[1] 0.092790150 9.304804346 9.924446415 0.003135142

##Endpoints for plotting segments
p0 = c(1900, 9.304804346)
p1 = c(1987, 9.304804346+0.003135142*(1987-1900))
p2 = c(1988, 9.924446415+0.003135142*(1988-1900))
p3 = c(2020, 9.924446415+0.003135142*(2020-1900))
pp= matrix(rbind(p0,p1,p2,p3), nrow=4)


#par(cex=1.5)
pdf("CET1900-2020FixedSlope-SCPT.pdf", width = 10, height = 5)
plot(Yr.Temp, xlab="Year", ylab=TeX("Average Temperature (^oC)"))
segments(pp[1,1], pp[1,2], pp[1+1,1], pp[1+1,2], col= "red",lwd=2)
segments(pp[3,1], pp[3,2], pp[3+1,1], pp[3+1,2], col= "red",lwd=2)
abline(v=1988, lty="dashed", col="red", lwd=2)
dev.off()
