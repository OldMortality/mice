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
stroke_size
female <- grepl('female',mice$Animal_type)
PLT <- grepl('PLT',mice$Animal_type)
Pound <- grepl('Pound',mice$Animal_type)

for (i in length(female)) {
  if (female[i]) {
    female[i] <- 1
  } else {
    female[i] <- 0
  }
}
for (i in length(PLT)) {
  if (PLT[i]) {
    PLT[i] <- 1
  } else {
    PLT[i] <- 0
  }
}
for (i in length(Pound)) {
  if (Pound[i]) {
    Pound[i] <- 1
  } else {
    Pound[i] <- 0
  }
}

gtt <- data.frame(female,PLT,Pound)


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
 

#compute some stuff, including eta_bar
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
#table(d$female,d$group)
#table(d$PLT,d$group)
#table(d$Pound,d$group)


# averages of logit-p by group (rows), time (cols)
logit <- function(x) {
  return(x/(1-x))
}

eta_bar <- matrix(NA, nrow = 6, ncol = 6)
for (i in 1:6) {
  for (j in 1:6) {
    eta_bar[i,j] <- mean(logit(d[which(d$group==i),j+3]))
  }
}
eta_bar 
print('done')
}

init_chains <- function() {
# initialise the chains
i1 <- 0
i2 <- 0.1
i3 <- -0.1
i4 <- 0.245
my_inits <- list(
  list(alpha0 = i1, alpha1 = i1, alpha2 = i1, alpha3 = i1, alpha4 = i1, alpha5 = i1,
       alpha6 = i1, alpha7 = i1, alpha8 = i1, sigma0 = 50, beta1 = i1, beta2 = i2, sigma =50,mu=200),
  list(alpha0 = i2, alpha1 = i2, alpha2 = i2 ,alpha3 = i2, alpha4 = i2, alpha5 = i2,
       alpha6 = i2, alpha7 = i2, alpha8 = i2, sigma0 = 50, beta1 = i2, beta2 = i1, sigma =32,mu=300),
  list(alpha0 = i3, alpha1 = i3, alpha2 = i3, alpha3 = i3, alpha4 = i3, alpha5 = i3,
       alpha6 = i3, alpha7 = i3, alpha8 = i3, sigma0 = 50, beta1 = -1 * i1, beta2=i3, sigma =56,mu=400),
  list(alpha0 = i4, alpha1 = i4, alpha2 = i4, alpha3 = i4, alpha4 = i4, alpha5 = i4, 
       alpha6 = i4, alpha7 = i4, alpha8 = i4, sigma0 = 50, beta1 = i3, beta2 = -1 * i4, sigma =80,mu=200)
) 
return(my_inits)
}
my_inits <- init_chains()
 

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

runjags <- function() {
   samples <- jags(model = mouseModel,
                data = list('R' = R_matrix,
                            'N' = N_matrix,
                            'female' = female,
                            'PLT' = PLT,
                            'Pound' = Pound,
                            'time_3' =  time_3,
                            'time_10' = time_10,
                            'time_17' = time_17,
                            'time_24' = time_24, 
                            'stroke_size' = stroke_size,
                            'N_mouse' = dim(R_matrix)[1]),
                parameters.to.save= parameters,
                n.chains=4,
                n.iter=5000,
                n.burnin=1000,
                DIC=T)
   return(samples) 
}


### leave female out 
parameters <- c("alpha0", "alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8",
                "delta1","delta2","delta3","delta4","delta5",
                "delta6","delta7","delta8","delta9","delta10",
                "sigma","p","b0","b","Tobs","Trep","eta_hat",  
                "chirep","chi","phi")

runjags <- function() {
  samples <- jags(model = mouseModelA,
                  data = list('R' = R_matrix,
                              'N' = N_matrix,
                              'PLT' = PLT,
                              'Pound' = Pound,
                              'time_3' =  time_3,
                              'time_10' = time_10,
                              'time_17' = time_17,
                              'time_24' = time_24, 
                              'stroke_size' = stroke_size,
                              'N_mouse' = dim(R_matrix)[1]),
                  parameters.to.save= parameters,
                  n.chains=4,
                  n.iter=5000,
                  n.burnin=1000,
                  DIC=T)
  return(samples) 
}



samples <-runjags()

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
gelman.diag(as.mcmc(samples))
gelman.plot(as.mcmc(samples)[,3:11])
#traceplot(as.mcmc(samples)[,],ask=F)


# get posteriors
{
phi <- samples$BUGSoutput$sims.list$phi
alpha0 <- samples$BUGSoutput$sims.list$alpha0
alpha1 <- samples$BUGSoutput$sims.list$alpha1
alpha2 <- samples$BUGSoutput$sims.list$alpha2
alpha3 <- samples$BUGSoutput$sims.list$alpha3
alpha4 <- samples$BUGSoutput$sims.list$alpha4
alpha5 <- samples$BUGSoutput$sims.list$alpha5
alpha6 <- samples$BUGSoutput$sims.list$alpha6
alpha7 <- samples$BUGSoutput$sims.list$alpha7
alpha8 <- samples$BUGSoutput$sims.list$alpha8
alpha9 <- samples$BUGSoutput$sims.list$alpha9
beta1 <- samples$BUGSoutput$sims.list$beta1
beta2 <- samples$BUGSoutput$sims.list$beta2
gamma1 <- samples$BUGSoutput$sims.list$gamma1
gamma2 <- samples$BUGSoutput$sims.list$gamma2
gamma3 <- samples$BUGSoutput$sims.list$gamma3
gamma4 <- samples$BUGSoutput$sims.list$gamma4
gamma5 <- samples$BUGSoutput$sims.list$gamma5
delta1 <- samples$BUGSoutput$sims.list$delta1
delta2 <- samples$BUGSoutput$sims.list$delta2
delta3 <- samples$BUGSoutput$sims.list$delta3
delta4 <- samples$BUGSoutput$sims.list$delta4
delta5 <- samples$BUGSoutput$sims.list$delta5
delta6 <- samples$BUGSoutput$sims.list$delta6
delta7 <- samples$BUGSoutput$sims.list$delta7
delta8 <- samples$BUGSoutput$sims.list$delta8
delta9 <- samples$BUGSoutput$sims.list$delta9
delta10 <- samples$BUGSoutput$sims.list$delta10
eta_hat <- samples$BUGSoutput$sims.list$eta_hat
p <- samples$BUGSoutput$sims.list$p
p[1,1,1]
b0 <- samples$BUGSoutput$sims.list$b0
b <- samples$BUGSoutput$sims.list$b
sigma <- samples$BUGSoutput$sims.list$sigma
sigma2 <- samples$BUGSoutput$sims.list$sigma2
p <- samples$BUGSoutput$sims.list$p
b0 <- samples$BUGSoutput$sims.list$b0
b <- samples$BUGSoutput$sims.list$b
 
DIC <- samples$BUGSoutput$DIC
chi <- samples$BUGSoutput$sims.list$chi
chirep <- samples$BUGSoutput$sims.list$chirep
print("got posteriors")
}

print(str_c("DIC: ",round(100*DIC)/100))

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

{
dim(eta_hat)
# eta_bar business
# eta_hat comes from the mcmc. It is 4000 * mice * time = 4000 x 35 x 6
#   model predictions for logit(p) over the
#     mcmc chain, by group (type*sex) and time period
eta_hat <-  array(0,dim=c(4000,6,6))
for (m in 1: 4000) {
  for (group in 1:6) {
    for (t in 1:6) {
      eta_hat[m,group,t] <- mean(eta_hat[m,which(d$group==group),t])  
    }
  }
}

# for the normal model, we have counts in eta_hat, but proportions in eta_bar.
eta_hat <- eta_hat / 1000

## compute Mean Square Pred Error
mspe <- numeric(4000)

for (m in 1:dim(eta_hat)[1]) {
  sum <- 0
  for (group in 1:6) {
    for (t in 1:6) {
      
      term <- sqrt((1/36) * (eta_hat[m,group,t] - eta_bar[group,t])^2)
      sum <- sum + term
      
    }
  }
  mspe[m] <- sum
  
}
hist(mspe,30,main="MSPE, summed over 6x6 sex*type*time")
median(mspe)
}

#plot(eta_hat[,1,1])


# so you can get a credible region just by sampling q for your favourite  mouse
# male PLT mouse at time 1

q <- expit(eta_hat)


# mouse 1 (male PLT mouse) at time 1
c(quantile(q[,1,1],0.025),quantile(q[,1,1],0.975))
# female PLT mouse at time 1
female[22]
c(quantile(q[,22,1],0.025),quantile(q[,22,1],0.975))
# female - male PLT mouse
p_diff <- q[,1,1] - q[,22,1]
c(quantile(p_diff,0.025),quantile(p_diff,0.975))
p_diff <- q[,1,2] - q[,11,2]
c(quantile(p_diff,0.025),quantile(p_diff,0.975))


# show posterior for p, and actual p
#  they are very tight
{
par(mfrow=c(3,2))
for (mouse in 13:24) {

  for (i in 1:6) {
    main = str_c('p for mouse ',mouse,' period ',i)
    {
    hist(p[,mouse,i],30,probability=T,main=main,xlim=c(0,1))
    abline(v=P_matrix[mouse,i],col='red')
    }
  }
}
}

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

{
par(mfrow=c(1,1))
cr <- plot_posteriors(alpha0,'a0')
cr <- plot_posteriors(alpha1,'a1, female')
cr <- plot_posteriors(alpha2,'a2, PLT')
cr <- plot_posteriors(alpha3,'a3, Pound')
cr <- plot_posteriors(alpha4,'a4, time 3')
cr <- plot_posteriors(alpha5,'a5, time 10')
cr <- plot_posteriors(alpha6,'a6, time 17')
cr <- plot_posteriors(alpha7,'a7, time 24')
cr <- plot_posteriors(alpha8,'a8, time 31')

cr <- plot_posteriors(beta1,'beta1, female*PLT')
cr <- plot_posteriors(beta2,'beta2, female*Pound')

cr <- plot_posteriors(gamma1,'gamma1, female*Pound')
cr <- plot_posteriors(gamma2,'gamma1, female*Pound')
cr <- plot_posteriors(gamma3,'gamma1, female*Pound')
cr <- plot_posteriors(gamma4,'gamma1, female*Pound')
cr <- plot_posteriors(gamma5,'gamma1, female*Pound')
}


#hist(p[,1,1],30,probability=T,main='mouse 1, period 1')
#hist(p[,2,1],30,probability=T,main='mouse 2, period 1')
}

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}


###
# Compare 3 types of male mice at time 0
###
 
  

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
  title <- c(str_c("CD19pos_B220 CR and data"))
  plot(0,0,pch='',xlim=c(0,32),ylim=c(0,0.6),xaxt = "n",
       main=title,xlab="time",ylab="proportion",grid(NA, 5, lwd = 2))
  axis(1, at=c(0,3,10,17,24,31), labels=c("0","3","10","17","24","31"))
  legend(0,0.59,c("WT","PLT","Pound"),
         lty=c(1,1,1),ncol=3,
         lwd=c(2.5,2.5),col=c("black","blue","red"))
  
  
  for (time in c(1,2,3,4,5,6)) {
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
          print(mouse)
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


par(mfrow=c(2,1))
hist(q[,37,1],30,main="female")
hist(q[,6,1],30,main="male")

for (t in 1:6) {
  
  pdiff <- eta_hat[,37,t] - eta_hat[,6,t]
  pdiff <- q[,37,t] - q[,6,t]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
}

########## difference by sex ############

diff_sex <- function(type,time) {
  
  
  if (type=="WT") {
    pf <- q[,37,time]  
    pm <- q[,6,time]  
  }
  if (type=="PLT") {
    pf <- q[,33,time]  
    pm <- q[,1,time]  
  }
  if (type=="Pound") {
    pf <- q[,25,time]  
    pm <- q[,8,time]  
  }

  pf <- expit(pf)
  pm <- expit(pm)
  return(pf-pm)  
}

par(mfrow=c(2,3))

count <- 0
for (type in c("WT","PLT","Pound")) {
  print(type)
  for (t in c(1,2,3,4,5,6)) {
  
   
    pdiff <- diff_sex(type,t)
    
    # 95% central CR.
    lwr <- quantile(pdiff,0.025)
    upr <- quantile(pdiff,0.975)
    title = str_c(type,"",t)
    if ( !(lwr < 0 & upr > 0) ) {
      # 0 is not in the CR
      count <- count + 1
      title <- str_c(title," ***")
      print(str_c(type, ' ',t,' ',lwr,' ',upr,"  ***"))
    } else {
      title = str_c(type,"",t)
      print(str_c(type,' ',t,' ',lwr,' ',upr))
    }
    
    
    
    hist(pdiff,30,xlim=c(-0.5,0.5),main=title,xlab="pf - pm")
    abline(v=0,col='red')
    segments(lwr,0.01,upr,0.01,col='red',lw=2)
    # is it different from 0?
    
  }
}
count
par(mfrow=c(1,1))

########################### difference by time #############

p_wt <- expit(q[,6,1] )
p_pound <- expit(q[,8,1] )
p_plt <- expit(q[,15,1] )
c(quantile(p_plt-p_wt,0.025),quantile(p_plt-p_wt,0.975))
c(quantile(p_pound-p_wt,0.025),quantile(p_pound-p_wt,0.975))

p_wt <- expit(q[,6,3] )
p_pound <- expit(q[,8,3] )
p_plt <- expit(q[,15,3] )
c(quantile(p_plt-p_wt,0.025),quantile(p_plt-p_wt,0.975))
c(quantile(p_pound-p_wt,0.025),quantile(p_pound-p_wt,0.975))


### difference by sex:

par(mfrow=c(3,1))
mf <- getTypicalMouse("female","WT")
mm <- getTypicalMouse("male","WT")
  
plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="WT mice, logit(pfemale)-logit(pmale)",
     xlab="time",ylab='diff')
abline(h=0,col='red')
for (t in 1:5) {
    pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
    print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
    segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="PLT mice, logit(pfemale)-logit(pmale)",
     xlab="time",ylab='diff')

abline(h=0,col='red')
mf <- getTypicalMouse("female","PLT")
mm <- getTypicalMouse("male","PLT")
for (t in 1:5) {
  pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="Pound mice, logit(pfemale)-logit(pmale)",
     xlab="time",ylab='diff')

abline(h=0,col='red')
mf <- getTypicalMouse("female","Pound")
mm <- getTypicalMouse("male","Pound")
for (t in 1:5) {
  pdiff <- eta_hat[,mf,t] - eta_hat[,mm,t]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

# compare types at each time, doesn't matter of which sex
Mwt <- getTypicalMouse("female","WT")
Mplt <- getTypicalMouse("female","PLT")
Mpound <- getTypicalMouse("female","Pound")

par(mfrow=c(2,1))
plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="PLT vs WT, logit(pPLT)-logit(pWT)",
     xlab="time",ylab='diff')
abline(h=0,col='red')
for (t in 1:5) {
  pdiff <- eta_hat[,Mplt,t] - eta_hat[,Mwt,t]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="Pound vs WT, logit(pPound)-logit(pWT)",
     xlab="time",ylab='diff')
abline(h=0,col='red')
for (t in 1:5) {
  pdiff <- eta_hat[,Mpound,t] - eta_hat[,Mwt,t]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}
ep <- eta_hat[,Mpound,t]
ewt <- eta_hat[,Mwt,t]
par(mfrow=c(2,1))
hist(exp(ep))
hist(exp(ewt))


par(mfrow=c(3,1))

plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="WT time vs time 1 ((logitp-t)-logit(p-1)",
     xlab="time",ylab='diff')
abline(h=0,col='red')


for (t in 2:5) {
  pdiff <- eta_hat[,Mwt,t] - eta_hat[,Mwt,1]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="PLT time vs time 1 ((logitp-t)-logit(p-1)",
     xlab="time",ylab='diff')
abline(h=0,col='red')

for (t in 2:5) {
  pdiff <- eta_hat[,Mplt,t] - eta_hat[,Mplt,1]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}
plot('',xlim=(c(0,5)),ylim=c(-3,3),
     main="Pound time vs time 1 ((logitp-t)-logit(p-1)",
     xlab="time",ylab='diff')
abline(h=0,col='red')

for (t in 2:5) {
  pdiff <- eta_hat[,Mpound,t] - eta_hat[,Mpound,1]
  print(c(quantile(pdiff,0.025),quantile(pdiff,0.975)))
  segments(t,quantile(pdiff,0.025),t,quantile(pdiff,0.975))
}

