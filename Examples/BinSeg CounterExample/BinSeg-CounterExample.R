##  JCL Review Paper
##  Binary segmentation counterexample
##  Author: Dr. Xueheng Shi
##  Date: 01/15/2022

library(ggplot2)


##Simulated time series with 3 alternating mean shifts
tsdat = data.frame(tt=seq(1,500), 
                Xt=rnorm(500)+c(rep(0,125),rep(1,125), rep(0,125), rep(1,125)))

mean.struct = data.frame(tt=seq(1,500), 
                         mean.shifts=c(rep(0,125),rep(1,125), rep(0,125), rep(1,125)))


##plot simulated data
#ggsave("BSCountExample.pdf", width = 20, height = 10, units = "cm")
p = ggplot() + 
  geom_line(data = tsdat, aes(x = tt, y = Xt), color = "blue") +
  geom_line(data = mean.struct , aes(x = tt, y = mean.shifts), color = "red") +
  xlab('t') +
  ylab('Xt')
p


##The simulation was performed on Top500 due to heavy computational cost of GA
##results were stored in 3MCPT_ar1_N=500_phi=0_Merged.csv
result = read.csv(file.choose(), header=T) #Load 3MCPT_ar1_N=500_phi=0_Merged.csv
#View(result)

##Number of repeats
n = length(result$phi.hat)


##Distributions (Table 3 added for Review paper)

round(table(result$mBIC)/n,3)
#3            4            5            6 
#0.9235294118 0.0568627451 0.0186274510 0.0009803922 

round(table(result$mMBIC)/n,3)
##3   4   5 
#985  25  10 

round(table(result$mMDL)/n,3)
###3   4   5   6   7   8   9 
#894  39  56  19   9   2   1 

round(table(result$mBS)/n,3)
###  0   1   2   3   4 
#####9 255   3 738  15 

round(table(result$mWBS)/n,3)
##  2   3   4   5   6   7 
####2 967  36  13   1   1  



##Absolute distance between estimated and true changepoint times
#ggsave("BSCounterExample_Dist_boxplots.pdf", width = 20, height = 10, units = "cm")
dist = data.frame(Methods = c(rep("BS",n),rep("GA-BIC",n),rep("GA-MDL",n),rep("GA-mBIC",n)),
                  Distance = c(result$dBS, result$dBIC, result$dMDL, result$dMBIC))                                                
dist$Methods = factor(dist$Methods , levels=c("BS", "GA-BIC", "GA-MDL", "GA-mBIC"))

dcpt = ggplot(dist, aes(x=Methods, y=Distance)) + 
  geom_boxplot()+
  #stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")+
  #ggtitle("Boxplots of the Distance between Estimated and True Changepoint Configuratons") + 
  ylab("Absolute Distance between estimated \n and true changepoint configuratons")
  #stat_summary(fun=mean, geom="line", aes(group=1), color="blue") 
dcpt




## number of changepoints
numb = data.frame(Methods = c(rep("BinSeg",n),rep("BIC",n),rep("MDL",n),rep("mBIC",n)),
                  Number = c(result$mBS, result$mBIC, result$mMDL, result$mMBIC))                                                
numb$Methods = factor(numb$Methods , levels=c("BinSeg", "BIC", "MDL", "mBIC"))

#ggsave("BSCounterExample_Num_boxplots.pdf", width = 20, height = 10, units = "cm")
ncpt = ggplot(numb, aes(x=Methods, y=Number)) + 
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")+ 
  geom_hline(yintercept=3,  linetype="dashed", color = "blue")
#stat_summary(fun=mean, geom="line", aes(group=1), color="blue") 
ncpt

