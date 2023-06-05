##Compute the Critical Value of SCUSUM Test
##Based the dist of integral of squared Brownian Bridges
##Based on Tolmatz 2002

##Author: Xueheng Shi
##Date: 09/11/2019

##Results for Table 1 in JCL Review

#N=100; lambda=0.45
# Dp = function(x){
#   (exp(-z^2/4)/gamma(0.5) )*(integrate(integrand, 0, Inf)$value)
# }


## CDF of integral of squared Brownian Bridges
##Based on eqn(2) in Tolmatz 2002
Fcdf = function(lambda){
  #lambda = 0.6
  N = 500
  item = NULL
  for(k in 1:N){
    integrand = function(x){ exp(-0.5*lambda*x^2)/sqrt(-x*sin(x)) }
    item[k]= (integrate(integrand, (2*k-1)*pi, (2*k)*pi)$value)
  }
  1-(2/pi)*sum(item)  #is used to search root
}


#Fcdf(0.4614)
# 0.9500124

##Computer Quantile
##change q to compute different quantiles
Fcdf.quantile = function(lambda, q = 0.90){
  N = 500
  item = NULL
  for(k in 1:N){
    integrand = function(x){ exp(-0.5*lambda*x^2)/sqrt(-x*sin(x)) }
    item[k]= (integrate(integrand, (2*k-1)*pi, (2*k)*pi)$value)
  }
  1-(2/pi)*sum(item)-q  #is used to search root
}

uniroot(Fcdf.quantile, c(0,1), tol = 0.0001)$root
# $root
# [1] 0.4613719