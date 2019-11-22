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
gtt <- data.frame(female,PLT,Pound)


N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))

##head(R_matrix)
##### is there life on earth? ######
#R_matrix[which(d$female),1] <- round(N_matrix[which(d$female),1] * 0.7)
###################################


## impute missing data. 
#
# replace each NA in N_matrix, R_matrix with the mean value of 
##    mice of the same type, female and period. If this is not

##    available, take the previous period.
length(which(is.na(N_matrix)) )
for (i in 1:40) {
  for (j in 1:6) {
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

# for printing only
P_matrix <- R_matrix/N_matrix



d <- data.frame(female,PLT,Pound, P_matrix)


{
par(mfrow=c(1,1))
plot(0,0,pch='',xlim=c(1,7),ylim=c(0,1),
     main='t=0,3: F wt, M wt; F plt  M plt  ; F Pound  M pound',
     xaxt = 'n')
axis(side = 1 ,at=c(1,2,3,5,6,7), 
labels=rep(c("WT","PLT","Pound"),2))
s1 <- subset(d,female==T & Pound==F & PLT==F)
s2 <- subset(d,female==F & Pound==F & PLT==F)
s3 <- subset(d,female==T & PLT==T)
s4 <- subset(d,female==F & PLT==T)
s5 <- subset(d,female==T & Pound==T)

s6 <- subset(d,female==F & Pound==T)
which(d$PLT & d$female )

# confusing, but R_.. are now the proportions
mean(s1$R_0)
mean(s2$R_0)
mean(s3$R_0)
mean(s4$R_0)
mean(s5$R_0)
mean(s6$R_0)

mean(s1$R_3)
mean(s2$R_3)
mean(s3$R_3)
mean(s4$R_3)
mean(s5$R_3)
mean(s6$R_3)

 
#t=0
points(rep(1.0,dim(s1)[1]),s1$R_0,col='red')
points(rep(1.1,dim(s2)[1]),s2$R_0)
points(rep(2.0,dim(s3)[1]),s3$R_0,col='red')
points(rep(2.1,dim(s4)[1]),s4$R_0)
points(rep(3.0,dim(s5)[1]),s5$R_0,col='red')
points(rep(3.1,dim(s6)[1]),s6$R_0)
text(2,0.8,"t=0",cex=2)

abline(v=4)
# t=3
points(rep(5,dim(s1)[1]),  s1$R_3,col='red')
points(rep(5.1,dim(s2)[1]),s2$R_3)
points(rep(6,dim(s3)[1]),s3$R_3,col='red')
points(rep(6.1,dim(s4)[1]),s4$R_3)
points(rep(7,dim(s5)[1]),s5$R_3,col='red')
points(rep(7.1,dim(s6)[1]),s6$R_3)
text(6,0.8,"t=1",cex=2)
}

####### End of plot ########

####### t tests ######################
{
## all types
t.test(R_0 ~ female,data=d)
t.test(R_3 ~ female,data=d)
t.test(R_10 ~ female,data=d)
t.test(R_17 ~ female,data=d)
t.test(R_24 ~ female,data=d)

d[which(d$Pound==F & d$PLT==F),]
t.test(R_0 ~ female,data=d[which(d$Pound==F & d$PLT==F),])
t.test(R_3 ~ female,data=d[which(d$Pound==F & d$PLT==F),])
t.test(R_10 ~ female,data=d[which(d$Pound==F & d$PLT==F),])
t.test(R_17 ~ female,data=d[which(d$Pound==F & d$PLT==F),])
t.test(R_24 ~ female,data=d[which(d$Pound==F & d$PLT==F),])
t.test(R_31 ~ female,data=d)


d[which(d$Pound),]
t.test(R_0 ~ female,data=d[which(d$Pound),])
t.test(R_3 ~ female,data=d[which(d$Pound),])
t.test(R_10 ~ female,data=d[which(d$Pound),])
t.test(R_17 ~ female,data=d[which(d$Pound),])
t.test(R_24 ~ female,data=d[which(d$Pound),])

d[which(d$PLT),]
t.test(R_0 ~ female,data=d[which(d$PLT),])
t.test(R_3 ~ female,data=d[which(d$PLT),])
t.test(R_10 ~ female,data=d[which(d$PLT),])
t.test(R_17 ~ female,data=d[which(d$PLT),])
t.test(R_24 ~ female,data=d[which(d$PLT),])
t.test(R_31 ~ female,data=d[which(d$PLT),])
t.test(R_31 ~ female,data=d[which(d$PLT),])
}


########################







# create dummy variables 
time_3  <- c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)



# initialise the chains
{
i1 <- 0
i2 <- 0.1
i3 <- -0.1
i4 <- 0.245
my_inits <- list(
  list(alpha0 = i1, alpha1 = i1, alpha2 = i1, alpha3 = i1, alpha4 = i1, alpha5 = i1,
       alpha6 = i1, alpha7 = i1, alpha8 = i1, sigma0 = 50, beta1 = i1, beta2 = i2, sigma =50),
  list(alpha0 = i2, alpha1 = i2, alpha2 = i2 ,alpha3 = i2, alpha4 = i2, alpha5 = i2,
       alpha6 = i2, alpha7 = i2, alpha8 = i2, sigma0 = 50, beta1 = i2, beta2 = i1, sigma =32),
  list(alpha0 = i3, alpha1 = i3, alpha2 = i3, alpha3 = i3, alpha4 = i3, alpha5 = i3,
       alpha6 = i3, alpha7 = i3, alpha8 = i3, sigma0 = 50, beta1 = -1 * i1, beta2=i3, sigma =56),
  list(alpha0 = i4, alpha1 = i4, alpha2 = i4, alpha3 = i4, alpha4 = i4, alpha5 = i4, 
       alpha6 = i4, alpha7 = i4, alpha8 = i4, sigma0 = 50, beta1 = i3, beta2 = -1 * i4, sigma =80)
)

} 
 
## jags models
source("/Users/micheldelange/Documents/mice/rscripts/mousemodels.r")




parameters <- c("alpha0", "alpha1","alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8","beta1","beta2",
                "gamma1","gamma2","gamma3","gamma4","gamma5",
                "sigma","sigma0","p","b0","b","Tobs","Trep","q")

runjags <- function() {

samples <- jags(model = mouseModel ,
                data = list('R' = R_matrix,
                            'N' = N_matrix,
                            'female' = female,
                            'PLT' = PLT,
                            'Pound' = Pound,
                            'time_3' =  time_3,
                            'time_10' = time_10,
                            'time_17' = time_17,
                            'time_24' = time_24,
                            'time_31' = time_31, 
                            'N_mouse' = 40 ),
                parameters,
                n.chains=4,
                n.iter=5000,
                n.burnin=1000,
                inits=my_inits,
                
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
 sgelman.diag(as.mcmc(samples))
gelman.plot(as.mcmc(samples)[,3:11])
traceplot(as.mcmc(samples)[,1:3],ask=F)

s <- as.mcmc(samples)


# get posteriors
alpha0_sample <- samples$BUGSoutput$sims.list$alpha0
alpha1_sample <- samples$BUGSoutput$sims.list$alpha1
alpha2_sample <- samples$BUGSoutput$sims.list$alpha2
alpha3_sample <- samples$BUGSoutput$sims.list$alpha3
alpha4_sample <- samples$BUGSoutput$sims.list$alpha4
alpha5_sample <- samples$BUGSoutput$sims.list$alpha5
alpha6_sample <- samples$BUGSoutput$sims.list$alpha6
alpha7_sample <- samples$BUGSoutput$sims.list$alpha7
alpha8_sample <- samples$BUGSoutput$sims.list$alpha8
beta1_sample <- samples$BUGSoutput$sims.list$beta1
beta2_sample <- samples$BUGSoutput$sims.list$beta2
gamma1_sample <- samples$BUGSoutput$sims.list$gamma1
gamma2_sample <- samples$BUGSoutput$sims.list$gamma2
gamma3_sample <- samples$BUGSoutput$sims.list$gamma3
gamma4_sample <- samples$BUGSoutput$sims.list$gamma4
gamma5_sample <- samples$BUGSoutput$sims.list$gamma5
p_sample <- samples$BUGSoutput$sims.list$p
p_sample[1,1,1]

sigma_sample <- samples$BUGSoutput$sims.list$sigma
sigma0_sample <- samples$BUGSoutput$sims.list$sigma0
p_sample <- samples$BUGSoutput$sims.list$p
b0_sample <- samples$BUGSoutput$sims.list$b0
b_sample <- samples$BUGSoutput$sims.list$b
q_sample <- samples$BUGSoutput$sims.list$q

# chain,mouse,period
for (mouse in 1:20) {
  print(q_sample[1,mouse,1])
}
# so you can get a credible region just by sampling q for your type of mouse
# male PLT mouse at time 1
c(quantile(q_sample[,1,1],0.025),quantile(q_sample[,1,1],0.975))
c(quantile(q_sample[,11,1],0.025),quantile(q_sample[,11,1],0.975))
# female - male PLT mouse
p_diff <- q_sample[,1,1] - q_sample[,11,1]
c(quantile(p_diff,0.025),quantile(p_diff,0.975))
# no difference between male, female PLT at t=0
p_diff <- q_sample[,1,2] - q_sample[,11,2]
c(quantile(p_diff,0.025),quantile(p_diff,0.975))
# no difference between male, female PLT at t=3
# male PLT definitely higher at time 3


# show posterior for p, and actual p
#  they are very tight
par(mfrow=c(3,2))
for (mouse in 13:24) {

  for (i in 1:6) {
    main = str_c('p for mouse ',mouse,' period ',i)
    {
    hist(p_sample[,mouse,i],30,probability=T,main=main,xlim=c(0,1))
    abline(v=P_matrix[mouse,i],col='red')
    }
  }
}


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


par(mfrow=c(1,1))
cr <- plot_posteriors(alpha0_sample,'a1, female')
cr <- plot_posteriors(alpha1_sample,'a1, female')
cr <- plot_posteriors(alpha2_sample,'a2, PLT')
cr <- plot_posteriors(alpha3_sample,'a3, Pound')
cr <- plot_posteriors(alpha4_sample,'a4, time 3')
cr <- plot_posteriors(alpha5_sample,'a5, time 10')
cr <- plot_posteriors(alpha6_sample,'a6, time 17')
cr <- plot_posteriors(alpha7_sample,'a7, time 24')
cr <- plot_posteriors(alpha8_sample,'a8, time 31')

cr <- plot_posteriors(beta1_sample,'beta1, female*PLT')
cr <- plot_posteriors(beta2_sample,'beta2, female*Pound')

cr <- plot_posteriors(gamma1_sample,'gamma1, female*Pound')
cr <- plot_posteriors(gamma2_sample,'gamma1, female*Pound')
cr <- plot_posteriors(gamma3_sample,'gamma1, female*Pound')
cr <- plot_posteriors(gamma4_sample,'gamma1, female*Pound')
cr <- plot_posteriors(gamma5_sample,'gamma1, female*Pound')



#hist(p_sample[,1,1],30,probability=T,main='mouse 1, period 1')
#hist(p_sample[,2,1],30,probability=T,main='mouse 2, period 1')


expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}


###
# Compare 3 types of male mice at time 0
###
 
  par(mfrow=c(1,1))
  logit_p_wild_type <- alpha0_sample
  logit_p_plt       <- alpha0_sample +   alpha2_sample 
  logit_p_pound     <- alpha0_sample +   alpha3_sample 

  p_wild_type <- expit(logit_p_wild_type)[,1]
  p_plt       <- expit(logit_p_plt)[,1]
  p_pound     <- expit(logit_p_pound)[,1]

  c(quantile(p_wild_type,0.025),quantile(p_wild_type,0.975))
  c(quantile(p_plt,0.025),quantile(p_plt,0.975))
  c(quantile(p_pound,0.025),quantile(p_pound,0.975))

  dat <- data.frame(p_wild_type,p_plt,p_pound)
  colnames(dat) <- c('WT',"PLT","POUND")
  dim(dat)
  data <- reshape2::melt(dat)
  dim(data)
  ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)

  ## add period 10
  logit_p_wild_type3 <- alpha0_sample + alpha1_sample * 1 +  alpha5_sample
  logit_p_plt3       <- alpha0_sample + alpha1_sample * 1 +  alpha2_sample + alpha5_sample 
  logit_p_pound3     <- alpha0_sample + alpha1_sample * 1 +  alpha3_sample + alpha5_sample

  p_wild_type3 <- expit(logit_p_wild_type3)[,1]
  p_plt3       <- expit(logit_p_plt3)[,1]
  p_pound3     <- expit(logit_p_pound3)[,1]

  dat <- data.frame(p_wild_type,
                  p_wild_type3,
                  p_plt,
                  p_plt3,
                  p_pound,
                  p_pound3)
  colnames(dat) <- c('WT',"WT3","PLT","PLT3","POUND","POUND3")
  dim(dat)
  data <- melt(dat)
  dim(data)
  ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5)
  table(data$variable)



 


###
###
### plot segments 
###
###


### You can use this function to get any quantile. For example
### for a female WT mouse at t=0, you would do:
###    parms <- c(1,0,0,0,0,0,0,0);
###    get_quantile(parms)
###
get_quantile <- function(parms) {
    p <- expit(alpha0_sample + 
               alpha1_sample * parms[1]+
               alpha2_sample * parms[2] +
               alpha3_sample * parms[3] +
               alpha4_sample * parms[4] +
               alpha5_sample * parms[5] +
               alpha6_sample * parms[6] +
               alpha7_sample * parms[7] +
               alpha8_sample * parms[8] + 
               beta1_sample * parms[9] + 
               beta2_sample * parms[10] +
               gamma1_sample * parms[4] * parms[1] +
               gamma2_sample * parms[5] * parms[1] +
               gamma3_sample * parms[6] * parms[1] +
               gamma4_sample * parms[7] * parms[1] + 
               gamma5_sample * parms[8] * parms[1]  
                 )
    return(c(quantile(p,0.025),quantile(p,0.975)))
}



plot_segments <- function(t,parms_wt,parms_plt,parms_pound) {
  q_wt <- get_quantile(parms_wt)
  segments(t,q_wt[1],t,q_wt[2])
  points(t,mean(q_wt),pch='.',cex=5)
  
  q_plt <- get_quantile(parms_plt)
  segments(t+delta,q_plt[1],t+delta,q_plt[2],col='blue')
  points(t+delta,mean(q_plt),col='blue',pch='.',cex=5)
  
  q_pound <- get_quantile(parms_pound)
  segments(t+delta2,q_pound[1],t+delta2,q_pound[2],col='red')
  points(t+delta2,mean(q_pound),col='red',pch='.',cex=5)
  
}

###
### female mice, 
##
{
par(mfrow=c(2,1))
title <- str_c(cell_type, ', female mice')
plot(0,0,pch='',xlim=c(0,32),ylim=c(0,0.5),main=title,xlab="time",ylab="proportion",
     grid(NA, 5, lwd = 2))
legend(15,0.5,c("WT","PLT","Pound"),
      lty=c(1,1,1),
      lwd=c(2.5,2.5),col=c("black","blue","red"))
delta <- 0.2
delta2 <- 2 * delta


### t=0
parms_wt  <-   c(1,0,0,0,0,0,0,0,0,0)
parms_plt <-   c(1,1,0,0,0,0,0,0,1,0)
parms_pound <- c(1,0,1,0,0,0,0,0,0,1)
plot_segments(t=0,parms_wt,parms_plt,parms_pound)

### t=3
parms_wt  <-   c(1,0,0,1,0,0,0,0,0,0)
parms_plt <-   c(1,1,0,1,0,0,0,0,1,0)
parms_pound <- c(1,0,1,1,0,0,0,0,0,1)
plot_segments(t=3,parms_wt,parms_plt,parms_pound)

### t=10
parms_wt  <-   c(1,0,0,0,1,0,0,0,0,0)
parms_plt <-   c(1,1,0,0,1,0,0,0,1,0)
parms_pound <- c(1,0,1,0,1,0,0,0,0,1)
plot_segments(t=10,parms_wt,parms_plt,parms_pound)

### t=17
parms_wt  <-   c(1,0,0,0,0,1,0,0,0,0)
parms_plt <-   c(1,1,0,0,0,1,0,0,1,0)
parms_pound <- c(1,0,1,0,0,1,0,0,0,1)
plot_segments(t=17,parms_wt,parms_plt,parms_pound)

### t=24
parms_wt  <-   c(1,0,0,0,0,0,1,0,0,0)
parms_plt <-   c(1,1,0,0,0,0,1,0,1,0)
parms_pound <- c(1,0,1,0,0,0,1,0,0,1)
plot_segments(t=24,parms_wt,parms_plt,parms_pound)

### t=31
parms_wt  <-   c(1,0,0,0,0,0,0,1,0,0)
parms_plt <-   c(1,1,0,0,0,0,0,1,1,0)
parms_pound <- c(1,0,1,0,0,0,0,1,0,1)
plot_segments(t=31,parms_wt,parms_plt,parms_pound)

###
### male mice
###

title <- str_c(cell_type,', male mice')
plot(0,0,pch='',xlim=c(0,32),ylim=c(0,0.5),main=title,xlab="time",ylab="proportion",grid(NA, 5, lwd = 2))

#legend(15,0.5,c("WT","PLT","Pound"),
#       lty=c(1,1,1),
#       lwd=c(2.5,2.5),col=c("black","blue","red"))
delta <- 0.2
delta2 <- 2 * delta


### t=0
parms_wt  <-   c(0,0,0,0,0,0,0,0,0,0)
parms_plt <-   c(0,1,0,0,0,0,0,0,0,0)
parms_pound <- c(0,0,1,0,0,0,0,0,0,0)
plot_segments(t=0,parms_wt,parms_plt,parms_pound)

### t=3
parms_wt  <-   c(0,0,0,1,0,0,0,0,0,0)
parms_plt <-   c(0,1,0,1,0,0,0,0,0,0)
parms_pound <- c(0,0,1,1,0,0,0,0,0,0)
plot_segments(t=3,parms_wt,parms_plt,parms_pound)

### t=10
parms_wt  <-   c(0,0,0,0,1,0,0,0,0,0)
parms_plt <-   c(0,1,0,0,1,0,0,0,0,0)
parms_pound <- c(0,0,1,0,1,0,0,0,0,0)
plot_segments(t=10,parms_wt,parms_plt,parms_pound)

### t=17
parms_wt  <-   c(0,0,0,0,0,1,0,0,0,0)
parms_plt <-   c(0,1,0,0,0,1,0,0,0,0)
parms_pound <- c(0,0,1,0,0,1,0,0,0,0)
plot_segments(t=17,parms_wt,parms_plt,parms_pound)

### t=24
parms_wt  <-   c(0,0,0,0,0,0,1,0,0,0)
parms_plt <-   c(0,1,0,0,0,0,1,0,0,0)
parms_pound <- c(0,0,1,0,0,0,1,0,0,0)
plot_segments(t=24,parms_wt,parms_plt,parms_pound)

### t=31
parms_wt  <-   c(0,0,0,0,0,0,0,1,0,0)
parms_plt <-   c(0,1,0,0,0,0,0,1,0,0)
parms_pound <- c(0,0,1,0,0,0,0,1,0,0)
plot_segments(t=31,parms_wt,parms_plt,parms_pound)
}


##### same idea, but now male vs female, for t=0
plot_segments_mf <- function(t) {
  
  par(mfrow=c(1,1))
  title <- str_c(cell_type,', male vs female mice')
  plot(0,0,pch='',xlim=c(0,10),ylim=c(0,1.0),main=title,xlab=c("time",0, " male on the left"),ylab="proportion",grid(NA, 5, lwd = 2))
  legend(8,0.2,c("WT","PLT","Pound"),
         lty=c(1,1,1),
         lwd=c(2.5,2.5),col=c("black","blue","red"))
  
  
  # male
  parms_wt  <-   c(0,0,0,0,0,0,0,0,0,0)
  parms_plt <-   c(0,1,0,0,0,0,0,0,0,0)
  parms_pound <- c(0,0,1,0,0,0,0,0,0,0)
  
  parms_wt[4] <- 1
  parms_plt[4] <- 1
  parms_pound[4] <- 1
  
  q_wt <- get_quantile(parms_wt)
  segments(t,q_wt[1],t,q_wt[2])
  points(t,mean(q_wt),pch='.',cex=5)
  
  q_plt <- get_quantile(parms_plt)
  segments(t+delta,q_plt[1],t+delta,q_plt[2],col='blue')
  points(t+delta,mean(q_plt),col='blue',pch='.',cex=5)
  
  q_pound <- get_quantile(parms_pound)
  segments(t+delta2,q_pound[1],t+delta2,q_pound[2],col='red')
  points(t+delta2,mean(q_pound),col='red',pch='.',cex=5)
  
  # female
  parms_wt  <-   c(1,0,0,0,0,0,0,0,0,0)
  parms_plt <-   c(1,1,0,0,0,0,0,0,1,0)
  parms_pound <- c(1,0,1,0,0,0,0,0,0,1)
  
  parms_wt[4] <- 1
  parms_plt[4] <- 1
  parms_pound[4] <- 1
  
  
  q_wt <- get_quantile(parms_wt)
  segments(5 + t,q_wt[1],5 + t,q_wt[2])
  points(5 +t,mean(q_wt),pch='.',cex=5)
  
  q_plt <- get_quantile(parms_plt)
  segments(5 +t+delta,q_plt[1],5 + t+delta,q_plt[2],col='blue')
  points(5 +t+delta,mean(q_plt),col='blue',pch='.',cex=5)
  
  q_pound <- get_quantile(parms_pound)
  segments(5 +t+delta2,q_pound[1],5 + t+delta2,q_pound[2],col='red')
  points(5 +t+delta2,mean(q_pound),col='red',pch='.',cex=5)
  
  ### now plot the data
  s1 <- subset(d,female==T & Pound==F & PLT==F)
  s2 <- subset(d,female==F & Pound==F & PLT==F)
  s3 <- subset(d,female==T & PLT==T)
  s4 <- subset(d,female==F & PLT==T)
  s5 <- subset(d,female==T & Pound==T)
  s6 <- subset(d,female==F & Pound==T)
  
  # male WT
  points(rep(2,0,dim(s2)[1]),s2$R_3,col='black')
  # male PLT
  points(rep(2.2,dim(s4)[1]),s4$R_3,col='blue')
  # male Pound
  points(rep(2.4,dim(s6)[1]),s6$R_3,col='red')
  
  # female WT
  points(rep(7.0,dim(s1)[1]),s1$R_3,col='black')
  # female PLT``
  points(rep(7.2,dim(s3)[1]),s3$R_3,col='blue')
  # female Pound
  points(rep(7.4,dim(s5)[1]),s5$R_3,col='red')
  
  
  text(1,0.6,'MALE')
  text(7,0.6,'FEMALE')
}
# go here
plot_segments_mf(0)
 

 



##########



 

logit_p_wild_type <- alpha0_sample 
logit_p_plt       <- alpha0_sample +   alpha2_sample 
logit_p_pound     <- alpha0_sample +   alpha3_sample 
p_wild_type <- expit(logit_p_wild_type)[,1]
p_plt       <- expit(logit_p_plt)[,1]
p_pound     <- expit(logit_p_pound)[,1]
c(quantile(p_wild_type,0.025),quantile(p_wild_type,0.975))
c(quantile(p_plt,0.025),quantile(p_plt,0.975))
c(quantile(p_pound,0.025),quantile(p_pound,0.975))


##

get_posterior_for_time <- function(time_period) {
  if (! time_period %in% c(0,3,10,17,24,31)) {
    print(c("error time period ",time_period))
  }
  if (time_period == 0) {
    return(rep(0,length(alpha4_sample)))
  }
  if (time_period == 3) {
    return(alpha4_sample)
  }
  if (time_period == 10) {
    return(alpha5_sample)
  }
  if (time_period == 17) {
    return(alpha6_sample)
  }
  if (time_period == 24) {
    return(alpha7_sample)
  }
  if (time_period == 31) {
    return(alpha8_sample)
  }
}
  



########## difference by sex #############

diff_sex <- function(type,time) {
  
  # female mice
  base <- alpha0_sample + alpha1_sample + get_posterior_for_time(time)
  
  if (type=="WT") {
    logit_pf <- base 
  }
  if (type=="PLT") {
    logit_pf <- base + alpha2_sample + beta1_sample
  }
  if (type=="Pound") {
    logit_pf <- base + alpha3_sample + beta2_sample
  }

  # male mice
  base <- alpha0_sample +  get_posterior_for_time(time)
  if (type=="WT") {
    logit_pm <-base 
  }
  if (type=="PLT") {
    logit_pm <- base + alpha2_sample 
  }
  if (type=="Pound") {
    logit_pm <- base + alpha3_sample
  }
  pf <- expit(logit_pf)
  pm <- expit(logit_pm)
  return(pf-pm)  
}

par(mfrow=c(3,6))

count <- 0
for (type in c("WT","PLT","Pound")) {
  # for (t in c(0,3,10,17,24,31)) {
  for (t in c(0)) {
    pdiff <- diff_sex(type,t)
    
    abline(v=0,col='red')
    lwr <- quantile(pdiff,0.025)
    upr <- quantile(pdiff,0.975)
    title = str_c(type,"",t)
    if ( !(lwr < 0 & upr > 0) ) {
      # 0 is not in the CR
      count <- count + 1
      title <- str_c(title," ***")
      print(str_c(type,'\t ',t,' ',lwr,' ',upr,"  ***"))
    } else {
      print(str_c(type,' ',t,' ',lwr,' ',upr))
    }
    
    
    
    hist(pdiff,30,xlim=c(-0.5,0.5),main=title,xlab="pf - pm")
    
    segments(lwr,0.01,upr,0.01,col='red',lw=2)
    # is it different from 0?
    
  }
}
count
par(mfrow=c(1,1))

########################### 

#verify these results

# female PLT
f_plt_0 <- alpha0_sample + alpha1_sample + alpha2_sample +  beta1_sample
# male PLT
m_plt_0 <- alpha0_sample +  alpha2_sample 


f_plt_0 <- expit(f_plt_0)
m_plt_0 <- expit(m_plt_0)
# female - male PLT

par(mfrow=c(2,1))
hist(f_plt_0,30)
hist(m_plt_0,30)

d_plt_0 <- f_plt_0 - m_plt_0
title = "p-female - p_male, PLT mice at time 0"
hist(d_plt_0,30,main=title)
abline(v=0,col='red')
segments(quantile(d_plt_0,0.025),0,quantile(d_plt_0,0.975),0,col='red',lw=2)


# female Pound
f_0 <- alpha1_sample + alpha3_sample +  beta2_sample
# male Pound
m_0 <- alpha1_sample + alpha3_sample 


f_0 <- expit(f_0)
m_0 <- expit(m_0)
# female - male PLT
d_0 <- f_0 - m_0
title = "p-female - p_male, Pound mice at time 0"
hist(d_0,30,main=title)
abline(v=0,col='red')
segments(quantile(d_0,0.025),0,quantile(d_0,0.975),0,col='red',lw=2)








