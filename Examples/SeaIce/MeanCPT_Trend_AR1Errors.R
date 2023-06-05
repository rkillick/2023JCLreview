
## Mean Shifts+Constant Trend+Global AR(1)Errors
MeanChange.BIC = function(loc.minseg, Xt=xt){

  N = length(Xt) #length of Xt
  loc.ind = rep(0, N)
  loc.ind[seq(1,N-1) %% 3 == 0] = loc.minseg
  
  Y = matrix(Xt, nrow=N, byrow=T) #Y: dependent variable
  
  m = sum(loc.ind) #Num of CPTs in a configuration
  D = matrix(c(rep(1,N), seq(1,N)), nrow=N, byrow=F) #design mat for zero cpt
  
  if(m == 0){ #Zero CPT
    beta.hat = solve((t(D)%*%D))%*%t(D)%*%Y ##Regression Coefs
    Rt = Y - D%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
    Rt.hat = c(0, phi.hat*Rt[-N])
    sigma2.hat = sum( (Rt-Rt.hat)^2 )/N
    BIC = N*log(sigma2.hat)+N+N*log(2*pi) + 4*log(N)
  }
  else{
    #tau.vec = loc.ind*(1:N) #CPT index
    tau = which(loc.ind == 1) #CPT locations
    tau.ext = c(1, tau, (N+1)) #CPT boundary 1 and N+1
    seg.len = diff(tau.ext)  #length of each segment
    
    Dm = matrix(0, nrow=N, ncol=m+1, byrow=T) #Design mat for mean shifts
    for(i in seq(1,m+1) ){
      #print(i)
      Dm[seq(tau.ext[i], tau.ext[i+1]-1), i]=1
    }
    M = cbind(Dm,seq(1,N))
    beta.hat = solve((t(M)%*%M))%*%t(M)%*%Y ##Regression Coefs
    Rt = Y - M%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
    Rt.hat = c(0, phi.hat*Rt[-N])
    sigma2.hat = sum( (Rt-Rt.hat)^2 )/N
    BIC = N*log(sigma2.hat)+N+N*log(2*pi)+(2*m+4)*log(N)
  }
  return(-BIC)
}






MeanCptTrendAR1.LogLKHD = function(tau, Xt=xt){

  N = length(Xt) #length of Xt
  loc.ind = rep(0, N)
  loc.ind[tau] = 1
  
  Y = matrix(Xt, nrow=N, byrow=T) #Y: dependent variable
  
  m = sum(loc.ind) #Num of CPTs in a configuration
  D = matrix(c(rep(1,N), seq(1,N)), nrow=N, byrow=F) #design mat for zero cpt
  
  if(m == 0){ #Zero CPT
    beta.hat = solve((t(D)%*%D))%*%t(D)%*%Y ##Regression Coefs
    Rt = Y - D%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
    Rt.hat = c(0, phi.hat*Rt[-N])
    sigma2.hat = sum( (Rt-Rt.hat)^2 )/N
    LogLKHD =  -0.5*N*log(sigma2.hat) -0.5*N-0.5*N*log(2*pi)
  }
  else{
    #tau.vec = loc.ind*(1:N) #CPT index
    tau = which(loc.ind == 1) #CPT locations
    tau.ext = c(1, tau, (N+1)) #CPT boundary 1 and N+1
    seg.len = diff(tau.ext)  #length of each segment
    
    Dm = matrix(0, nrow=N, ncol=m+1, byrow=T) #Design mat for mean shifts
    for(i in seq(1,m+1) ){
      #print(i)
      Dm[seq(tau.ext[i], tau.ext[i+1]-1), i]=1
    }
    M = cbind(Dm,seq(1,N))
    beta.hat = solve((t(M)%*%M))%*%t(M)%*%Y ##Regression Coefs
    Rt = Y - M%*%beta.hat  ##AR(1) residuals
    phi.hat = sum( Rt[-N]*Rt[-1] )/sum( Rt[-1]^2 )
    Rt.hat = c(0, phi.hat*Rt[-N])
    sigma2.hat = sum( (Rt-Rt.hat)^2 )/N
    LogLKHD =  -0.5*N*log(sigma2.hat) -0.5*N-0.5*N*log(2*pi)
  }
  return(c(sigma2.hat,LogLKHD))
}

#MeanCptTrendAR1.LogLKHD(c(219), Xt=xt)
##[1]    0.3249955 -310.2218880