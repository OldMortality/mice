##
##  R2jags
##  
##  Binomial model
##
##  This one works

library(R2jags)
library(stringr)
library(ggplot2)
library(reshape2)

#
# usually called from do_mice.r
#
#

female <- grepl('female',mice$Animal_type)
PLT <- grepl('PLT',mice$Animal_type)
Pound <- grepl('Pound',mice$Animal_type)

{
for (i in 1:length(female)) {
  if (female[i]) {
    female[i] <- 1
  } else {
    female[i] <- 0
  }
}
for (i in 1:length(PLT)) {
  if (PLT[i]) {
    PLT[i] <- 1
  } else {
    PLT[i] <- 0
  }
}
for (i in 1:length(Pound)) {
  if (Pound[i]) {
    Pound[i] <- 1
  } else {
    Pound[i] <- 0
  }
}

gtt <- data.frame(female,PLT,Pound)
}

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24))


# impute missing data
{
## impute missing data. 
#
# replace each NA in N_matrix, R_matrix with the mean value of 
##    mice of the same type, female and period. If this is not

##    available, take the previous period.
length(which(is.na(N_matrix)) )
for (i in 1:dim(mice)[1]) {
  for (j in 1:5) {
    if (is.na(N_matrix[i,j])) {
      this_female <- female[i]
      this_Pound <- Pound[i]
      # which mice are of the same type and gender and this mouse?
      this_PLT <- PLT[i]
      similar <- which(gtt$female==this_female & gtt$PLT==this_PLT & gtt$Pound==this_Pound)
      #length(N_matrix[similar,k])
      N_imp <- round(mean(N_matrix[similar,j],na.rm=T))
      R_imp <- round(mean(R_matrix[similar,j],na.rm=T))
       
      if (is.nan(N_imp)) {
        # if that did not work, we take mean of previous time period
        N_imp <- round(mean(N_matrix[similar,j-1],na.rm=T))
        R_imp <- round(mean(R_matrix[similar,j-1],na.rm=T))
      }
      if (is.nan(N_imp)) {
        # one period further back
        N_imp <- round(mean(N_matrix[similar,j-2],na.rm=T))
        R_imp <- round(mean(R_matrix[similar,j-2],na.rm=T))
      }
      if (is.nan(N_imp)) {
        # this will do
        N_imp <- round(mean(N_matrix[similar,j-3],na.rm=T))
        R_imp <- round(mean(R_matrix[similar,j-3],na.rm=T))
         
      }
      if (is.nan(N_imp)) {
        print(str_c('error imputing data',i,j))
      }
      N_matrix[i,j] <- N_imp
      R_matrix[i,j] <- R_imp
      
    }
  }
 
}
which(is.na(N_matrix)) 
}
 

#compute some stuff
{
  # create dummy variables 
  time_3  <- c(0,1,0,0,0,0)
  time_10 <- c(0,0,1,0,0,0)
  time_17 <- c(0,0,0,1,0,0)
  time_24 <- c(0,0,0,0,1,0)

  
P_matrix <- R_matrix/N_matrix
d <- data.frame(female,PLT,Pound, P_matrix)
d$group <- NA
d[which(d$female & ( d$PLT==F & d$Pound==F)),"group"] <- 1
d[which(!d$female & ( d$PLT==F & d$Pound==F)),"group"] <- 2
d[which(d$female & d$PLT),"group"] <- 3
d[which(!d$female & d$PLT),"group"] <- 4
d[which(d$female & d$Pound),"group"] <- 5
d[which(!d$female & d$Pound),"group"] <- 6
}

 
## jags models
source("/Users/micheldelange/Documents/mice/rscripts/mousemodels.r")

parameters <- c("alpha0", "alpha1","alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8",
                "beta1","beta2",
                "gamma1","gamma2","gamma3","gamma4","gamma5",
                "delta1","delta2","delta3","delta4","delta5",
                "delta6","delta7","delta8","delta9","delta10",
                "sigma","p","b0","b","Tobs","Trep","eta_hat",  
                "chirep","chi","phi")

runjags <- function(theJagsModel) {
   samples <- jags(model = theJagsModel,
                data = list('R' = R_matrix,
                            'N' = N_matrix,
                            'female' = female,
                            'PLT' = PLT,
                            'Pound' = Pound,
                            'time_3' =  time_3,
                            'time_10' = time_10,
                            'time_17' = time_17,
                            'time_24' = time_24,
                            'N_mouse' = dim(R_matrix)[1]),
                parameters.to.save= parameters,
                n.chains=4,
                n.iter=5000,
                n.burnin=1000,
                DIC=T)
   return(samples) 
}




# make sure you run the right one
# this is the one with sex
samples <-runjags(mouseModel)

par(mfrow=c(2,1))
coda=samples$BUGSoutput$sims.matrix
#colnames(coda)
Tobs=coda[,1]
Trep=coda[,2]
bpv=mean(Trep>Tobs)
bpv
fletch <- Tobs/Trep
sqrt(var(fletch))
hist(fletch,30)
plot(Trep ~ Tobs)
abline(0,1,col='red')
#gelman.diag(as.mcmc(samples))
#gelman.plot(as.mcmc(samples)[,3:11])
#traceplot(as.mcmc(samples)[,],ask=F)


# get posteriors
source("/Users/micheldelange/Documents/mice/rscripts/get_posteriors.r") 


#print(str_c("DIC: ",round(100*DIC)/100))



# so you can get a credible region just by sampling q for your favourite  mouse
# male PLT mouse at time 1. This eta_hat does not include the stroke size, so it
# is for your average stroke size.
q <- expit(eta_hat)



# show posterior for p, and actual p
#  they are very tight
#{
#par(mfrow=c(3,2))
#for (mouse in 13:24) {
#  for (i in 1:5) {
#    main = str_c('p for mouse ',mouse,' period ',i)
#    {
#    hist(p[,mouse,i],30,probability=T,main=main,xlim=c(0,1))
#    abline(v=P_matrix[mouse,i],col='red')
##    }
#  }
#}
#}

# plot posteriors
{
par(mfrow=c(1,2))
plot_posterior <- function(s,title,redline) {
  hist(s,30,probability=T,main=title)
  abline(v=redline,col='red')
  # 95% CR in blue
  lower <- quantile(s,0.025)
  upper <- quantile(s,0.975)
  segments(x0 =lower,y0=0.1,x1=upper,y1=0.1,col="blue")
  return(c(lower,upper))

}


plot_posteriors<- function(s, title) {
  plot_posterior(s,str_c(title,', parameter'),redline=0)
  cr <- plot_posterior(exp(s),str_c(title,', odds-ratio'),redline=1)
  print(cr)
  return(cr)
}


# plot posteriors
{
#par(mfrow=c(1,1))
#cr <- plot_posteriors(alpha0,'a0')
#cr <- plot_posteriors(alpha1,'a1, female')
#cr <- plot_posteriors(alpha2,'a2, PLT')
#cr <- plot_posteriors(alpha3,'a3, Pound')
#cr <- plot_posteriors(alpha4,'a4, time 3')
#cr <- plot_posteriors(alpha5,'a5, time 10')
#cr <- plot_posteriors(alpha6,'a6, time 17')
#cr <- plot_posteriors(alpha7,'a7, time 24')
#cr <- plot_posteriors(alpha8,'a8, time 31')

##cr <- plot_posteriors(beta1,'beta1, female*PLT')
#cr <- plot_posteriors(beta2,'beta2, female*Pound')

#cr <- plot_posteriors(gamma1,'gamma1, female*Pound')
#cr <- plot_posteriors(gamma2,'gamma1, female*Pound')
#cr <- plot_posteriors(gamma3,'gamma1, female*Pound')
#cr <- plot_posteriors(gamma4,'gamma1, female*Pound')
#cr <- plot_posteriors(gamma5,'gamma1, female*Pound')
}

#par(mfrow=c(2,1))
#hist(expit(alpha0),30,main="male WT t0 mean stroke size",probability = T)
#low <- quantile(expit(alpha0),0.025)
#upp <- quantile(expit(alpha0),0.975)
#segments(low,0,upp,0,col='red')
#hist(expit(alpha0 + phi),30,main="male WT t0 +1sd stroke size",probability = T)
#low <- quantile(expit(alpha0+ phi),0.025)
#upp <- quantile(expit(alpha0+ phi),0.975)
#segments(low,0,upp,0,col='red')
#
}




 

###
###
### plot CR and data
###
###


getTypicalMouse <- function(sex,type) {
  if (sex=="female" & type=="WT") { return(23) }
  if (sex=="male" & type=="WT") { return(20) }
  if (sex=="female" & type=="PLT") { return(10) }
  if (sex=="male" & type=="PLT") { return(1) }
  if (sex=="female" & type=="Pound") { return(11) }
  if (sex=="male" & type=="Pound") { return(9) }
}






plot_CR_data <- function() {
  
  # for plotting the data
  s1 <- subset(d,female==T & Pound==F & PLT==F)
  s2 <- subset(d,female==F & Pound==F & PLT==F)
  s3 <- subset(d,female==T & PLT==T)
  s4 <- subset(d,female==F & PLT==T)
  s5 <- subset(d,female==T & Pound==T)
  s6 <- subset(d,female==F & Pound==T) 
  
  par(mfrow=c(1,1))
  title <- cell_type
  plot(0,0,pch='',xlim=c(0,32),ylim=c(0,1),xaxt = "n",
       main=title,xlab="time",ylab="proportion",grid(NA, 5, lwd = 2))
  axis(1, at=c(0,3,10,17,24,31), labels=c("0","3","10","17","24","31"))
  legend(0,0.9,c("WT","PLT","Pound"),
         lty=c(1,1,1),ncol=3,
         lwd=c(2.5,2.5),col=c("black","blue","red"))
  
  
  for (time in c(1,2,3,4,5)) {
    if (time == 1) { t = 0}
    if (time == 2) { t = 3}
    if (time == 3) { t = 10}
    if (time == 4) { t = 17}
    if (time == 5) { t = 24}
    if (time == 6) { t = 31}
  
    for (type in c("WT","PLT","Pound")) {
      
      for (sex in c("female","male")) {
        
        if (type == "WT") { col="black"}
        if (type == "PLT") { col="blue"}
        if (type == "Pound") { col="red"} 
        
          mouse <- getTypicalMouse(sex,type)
         
          segments(t,quantile(q[,mouse,time],0.025),t,quantile(q[,mouse,time ],0.975),col=col)
          
          ##
          ## plot the data
          ##
          if (sex=="female" & type == "WT") {
            points(rep(t,dim(s1)[1]),s1[,3+time],col=col)
          }
          if (sex=="male" & type == "WT") {
            points(rep(t,dim(s2)[1]),s2[,3+time],col=col)
          }
          if (sex=="female" & type == "PLT") {
            points(rep(t,dim(s3)[1]),s3[,3+time],col=col)
          }
          if (sex=="male" & type == "PLT") {
            points(rep(t,dim(s4)[1]),s4[,3+time],col=col)
          }
          if (sex=="female" & type == "Pound") {
            points(rep(t,dim(s5)[1]),s5[,3+time],col=col)
          }
          if (sex=="male" & type == "Pound") {
            points(rep(t,dim(s6)[1]),s6[,3+time],col=col)
          }
          t <- t + 0.5
        }      
      }
    }
}  
plot_CR_data()  


###
### difference by sex:
###   3 plots in one: 1 plot for each type:
###
###
{
par(mfrow=c(1,3))
mf <- getTypicalMouse("female","WT")
mm <- getTypicalMouse("male","WT")
## WT
plot('',xlim=(c(1,5)),ylim=c(-3,3),xaxt = "n",
     main=str_c(cell_type, "female vs male WT"),
     xlab="time",ylab='logit(pfemale)-logit(pmale)')
axis(1, at=c(1,2,3,4,5), labels=c("0","3","10","17","24"))

grid(NA, 5, lwd = 2)
abline(h=0,col='red')
for (t in 1:5) {
    pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
    #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
    segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}
## PLT
plot('',xlim=(c(1,5)),ylim=c(-3,3),,xaxt = "n",
     main=str_c(cell_type, "female vs male PLT"),
     xlab="time",ylab='logit(pfemale)-logit(pmale)')
axis(1, at=c(1,2,3,4,5), labels=c("0","3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')
mf <- getTypicalMouse("female","PLT")
mm <- getTypicalMouse("male","PLT")
for (t in 1:5) {
  pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}
## Pound
plot('',xlim=(c(1,5)),ylim=c(-3,3),xaxt='n',
     main=str_c(cell_type, "female vs male Pound"),
     xlab="time",ylab='logit(pfemale)-logit(pmale)')
axis(1, at=c(1,2,3,4,5), labels=c("0","3","10","17","24"))

grid(NA, 5, lwd = 2)
abline(h=0,col='red')
mf <- getTypicalMouse("female","Pound")
mm <- getTypicalMouse("male","Pound")
for (t in 1:5) {
  pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}
}


 
#hist(fletch,30)
#par(mfrow=c(1,1))
#plot(Trep ~ Tobs)
#abline(0,1,col='red')
#gelman.diag(as.mcmc(samples))
#gelman.plot(as.mcmc(samples)[,3:11])
#traceplot(as.mcmc(samples)[,],ask=F)


# get posteriors
#source("/Users/micheldelange/Documents/mice/rscripts/get_posteriors.r") 
#q <- expit(eta_hat)


#
# compare types at each time, take a female mouse
#
{
n_tests <- 5
eps <- 0.05
    
if (buon) {
    eps <- 0.05 / n_tests
} 

Mwt <- getTypicalMouse("female","WT")
Mplt <- getTypicalMouse("female","PLT")
Mpound <- getTypicalMouse("female","Pound")

par(mfrow=c(1,2))
plot('',xlim=(c(1,5)),ylim=c(-3,3),xaxt='n',
     main=str_c(cell_type, "logit PLT vs WT"),
     xlab="time",ylab='diff (logit scale)')
axis(1, at=c(1,2,3,4,5), labels=c("0","3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')
for (t in 1:5) {
  pdiff <- eta_hat[,Mplt,t] - eta_hat[,Mwt,t]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,eps/2),t,quantile(pdiff,1-eps/2))
}

plot('',xlim=(c(1,5)),ylim=c(-3,3),xaxt='n',
     main=str_c(cell_type, " logit Pound vs WT"),
     xlab="time",ylab='diff (logit scale')
axis(1, at=c(1,2,3,4,5), labels=c("0","3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')
for (t in 1:5) {
  pdiff <- eta_hat[,Mpound,t] - eta_hat[,Mwt,t]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,eps/2),t,quantile(pdiff,1-eps/2))
}
}

# plots at probability scale.
{
par(mfrow=c(5,3))
for (t in 1:5) {
  # eta_hat is at the logit scale
  ewt <- eta_hat[,Mwt,t]
  eplt <- eta_hat[,Mplt,t]
  ep <- eta_hat[,Mpound,t]
  #hist(expit(ewt),30,xlim=c(0,1), main=str_c("Probablity of cell, type WT time: ",t))
  #hist(expit(eplt),30,xlim=c(0,1),main=str_c("Probablity of cell, type PLT time: ",t))
  #hist(expit(ep),30,xlim=c(0,1),  main=str_c("Probablity of cell, type Pound time: ",t))
}
}

## difference vs time 0.
{
n_tests <- 4
eps <- 0.05
if (buon) {
  eps <- 0.05 / n_tests
} 

par(mfrow=c(1,3))

plot('',xlim=(c(2,5)),ylim=c(-3,3),xaxt="n",
     main=str_c(cell_type, "WT time vs t0"),
     xlab="time",ylab='diff')
axis(1, at=c(2,3,4,5), labels=c("3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')


for (t in 2:5) {
  pdiff <- eta_hat[,Mwt,t] - eta_hat[,Mwt,1]
  #print(c(quantile(pdiff,eps/2),quantile(pdiff,1-eps/2)))
  segments(t,quantile(pdiff,eps/2),t,quantile(pdiff,1-eps/2))
}


plot('',xlim=(c(2,5)),ylim=c(-3,3),xaxt='n',
     main=str_c(cell_type, " PLT time vs t0"),
     xlab="time",ylab='diff')
axis(1, at=c(2,3,4,5), labels=c("3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')

for (t in 2:5) {
  pdiff <- eta_hat[,Mplt,t] - eta_hat[,Mplt,1]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,eps/2),t,quantile(pdiff,1-eps/2))
}

plot('',xlim=(c(2,5)),ylim=c(-3,3),xaxt="n",
     main=str_c(cell_type, "Pound vs t1"),
     xlab="time",ylab='diff')
axis(1, at=c(2,3,4,5), labels=c("3","10","17","24"))
grid(NA, 5, lwd = 2)
abline(h=0,col='red')

for (t in 2:5) {
  pdiff <- eta_hat[,Mpound,t] - eta_hat[,Mpound,1]
  #print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,eps/2),t,quantile(pdiff,1-eps/2))
}

}

