## Use statistics (obtained by QQ Lu) from Mott series to generate data
## The data contain a seasonal component
## Binary segmentation was tested on it
## Author: Xueheng Shi
## Date: 04/21/2021


#library(TSA)
#library(readr)
library(wbs)


##Monthly average temperatures, computed by QQ Lu
mon.avg = c(-11.210286,  -7.933000,  -2.340857,   5.606232,  12.198143,
17.282754,  21.080143,  19.994000,  13.798143,   7.156377,
-2.196000,  -8.200286)

##Monthly standard errors, computed by QQ Lu
mon.sd = c(4.782172, 4.623396, 3.341025, 2.068479, 2.005662, 1.901481,
           1.847303, 1.761586, 1.987364, 1.696544, 2.744698, 3.635875)

##Number of days in 12 months
mon.days = c(31,28,31,30,31,30,31,31,30,31,30,31)
#sum(mon.days) ##365


##Store results
m.vec = NULL
tau.vec = NULL


## Simulated data has two periods, T=365
for(k in 1:100){
  
  dailytemp1 = NULL
  for(i in 1:12){
    #i = (k-1)%%12 + 1
    temp = rnorm(mon.days[i], mean=mon.avg[i],sd=mon.sd[i] )
    dailytemp1 = c(dailytemp1, temp)
  }
  
  
  dailytemp2 = NULL
  for(i in 1:12){
    #i = (k-1)%%12 + 1
    temp = rnorm(mon.days[i], mean=mon.avg[i],sd=mon.sd[i] )
    dailytemp2 = c(dailytemp2, temp)
  }
  

  xt = c(dailytemp1, dailytemp2)
  
  bs = wbs::sbs(xt)
  bscpt = wbs::changepoints(bs)
  mBS = bscpt$no.cpt.th #number of CPT
  tBS = as.vector(unlist(bscpt$cpt.th)) #CPT location
  #t.BS = ifelse(mBS==0, "NA", toString(tBS) )
  m.vec = rbind(m.vec, mBS) 
  #tau.vec = rbind(tau.vec, tBS)
  print(k)
}

pdf("Simulated-Mott.pdf", width = 10, height = 5)
plot(as.ts(xt), ylab="Daily Temperautre in Degree of Celsius")
op = par(mfrow=c(1,1),mgp=c(2.5, 0.8,0))
plot(c(1, 750), y=c(-25, 30),type='n',xlab="Day",ylab="Daily Temperautres(Degrees Celsius)", axes=F)
axis(1, at=c(0,150, 300, 450, 600, 750), labels=c(0,150, 300, 450, 600, 750), tck=-0.02)
#axis(1, at=seq(1850, 2010, 10), labels=F,cex=0.1, tck=-0.01)
axis(2, at=seq(-20, 20, 10),labels=seq(-20, 20, 10), tck=-0.02)
#axis(2, at=seq(0, 110, 10), labels=F,cex=0.1, tck=-0.01)
box()
lines(1:730, as.ts(xt))
dev.off()

## Summarize Binary Segmentation Detection
boxplot(m.vec)
summary(m.vec)

#par(mfrow=c(1,2), mar=c(2,2,2,2))
#plot(as.ts(xt), ylab="Daily Temperautre in Degree of Celsius")

pdf("Simulated-BinSeg-Mott.pdf", width = 3, height = 6)
boxplot(m.vec)
dev.off()






