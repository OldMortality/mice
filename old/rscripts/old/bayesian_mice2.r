##
##  R2jags
##  Long data
##  Binomial model
##
##

library(R2jags)
library(stringr)
mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]

# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

N_0 <- mice$Live_cells_0
R_0 <- mice$CD19pos_B220pos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD19pos_B220pos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD19pos_B220pos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD19pos_B220pos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD19pos_B220pos_24
N_31 <- mice$Live_cells_31
R_31 <- mice$CD19pos_B220pos_31

gender <- grepl('female',mice$Animal_type)
type_1 <- grepl('PLT',mice$Animal_type)
type_2 <- grepl('Pound',mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))
P_matrix <- R_matrix/N_matrix


time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)


 
# long notation (i.e. observations as rows)
per_long <- vector()
mouse_long <- vector()
gender_long <- vector()
type_1_long <- vector()
type_2_long <- vector()
time_3_long <- vector()
time_10_long <- vector()
time_17_long <- vector()
time_24_long <- vector()
time_31_long <- vector()

N_long <- vector()
R_long <- vector()
 

counter = 0
skipped = 0
for (i in 1: N) {
  for (j in 1:6) {
    
    counter <- counter + 1
    if (is.na(N_matrix[i,j])) {
      skipped <- skipped + 1
      next
    }
    per_long    <- c(per_long,j)
    mouse_long  <- c(mouse_long,i)
    gender_long <- c(gender_long,gender[i])
    type_1_long <- c(type_1_long,type_1[i])
    type_2_long <- c(type_2_long, type_2[i])
    N_long      <- c(N_long,N_matrix[i,j])
    R_long      <- c(R_long,R_matrix[i,j])
    time_3_long <- c(time_3_long,time_3[j])
    time_10_long <-c(time_10_long,time_10[j])
    time_17_long <- c(time_17_long,time_17[j])
    time_24_long <- c(time_24_long,time_24[j])
    time_31_long <- c(time_31_long,time_31[j])
  }
  
}

long_data <- data.frame(mouse_long,
                   per_long,
                   gender_long,
                   type_1_long,
                   type_2_long,
                   N_long,
                   R_long,
                   time_3_long,
                   time_10_long,
                   time_17_long,
                   time_24_long,
                   time_31_long)

write.csv(long_data,"/Users/micheldelange/Documents/mice/data/long_data.csv")

mouseModel <- function()
{
  for ( i in 1 : N) {
      R_long[i] ~ dbin(p[mouse[i], period[i]],N_long[i])
      b[mouse[i], period[i]] ~ dnorm(0.0, tau); 
      logit(p[mouse[i], period[i]]) <- alpha0 + alpha1 * gender[i] +
        alpha2 * type_1[i] +
        alpha3 * type_2[i] +
        alpha4 * time3[i] +
        alpha5 * time10[i] +
        alpha6 * time17[i] +
        alpha7 * time24[i] +
        alpha8 * time31[i] +
        b[mouse[i],period[i]]
    }
  
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}


## just checking tau
r.debug <- runif(1000,0,100)
t.debug <- r.debug ^ -2
# very low precision
t.debug


time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)

# initialise the chains
i1 <- 0
i2 <- 0.1
i3 <- -0.1
i4 <- 0.245

my_inits <- list(
  list(alpha0 = i1, 
       alpha1 = i1,
       alpha2 = i1 ,
       alpha3 = i1,
       alpha4 = i1,
       alpha5 = i1,
       alpha6 = i1,
       alpha7 = i1,
       alpha8 = i1,
       sigma =50),
  list(alpha0 = i2, 
       alpha1 = i2,
       alpha2 = i2 ,
       alpha3 = i2,
       alpha4 = i2,
       alpha5 = i2,
       alpha6 = i2,
       alpha7 = i2,
       alpha8 = i2,
       sigma =32),
  list(alpha0 = i3, 
       alpha1 = i3,
       alpha2 = i3,
       alpha3 = i3,
       alpha4 = i3,
       alpha5 = i3,
       alpha6 = i3,
       alpha7 = i3,
       alpha8 = i3,
       sigma =56),
  list(alpha0 = i4, 
       alpha1 = i4,
       alpha2 = i4,
       alpha3 = i4,
       alpha4 = i4,
       alpha5 = i4,
       alpha6 = i4,
       alpha7 = i4,
       alpha8 = i4,
       sigma =80)
  
)


mouseModel2 <- function() 
{
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
  
  
  phi[mouse[i],period[i]]~ dgamma(1.0E-3, 1.0E-3)
  
  for (i in 1:N) {
    
    R_long[i] ~ dbinom(p[mouse[i], period[i]], N_long[i])
    
    p[mouse[i], period[i]] ~ dbeta(alpha[mouse[i], period[i]], beta[mouse[i], period[i]])
    
    alpha[mouse[i], period[i]] = phi[mouse[i], period[i]] * mu[mouse[i], period[i]]
    
    beta[mouse[i], period[i]] = phi[mouse[i], period[i]] * (1 - mu[mouse[i], period[i]])
    
    logit(mu[mouse[i], period[i]]) <- alpha0 + alpha1 * gender[i] +
      alpha2 * type_1[i] +
      alpha3 * type_2[i] +
      alpha4 * time3[i] +
      alpha5 * time10[i] +
      alpha6 * time17[i] +
      alpha7 * time24[i] +
      alpha8 * time31[i] +
      b[mouse[i],period[i]]
  }
  

  
}

parameters <- c("alpha0", "alpha1","alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8","sigma","p","b")
samples <- jags(model = mouseModel2() ,
                data = list('R_long' = R_long,
                            'N_long' = N_long,
                            'period' = per_long,
                            'mouse'  = mouse_long,
                            'gender' = gender_long,
                            'type_1' = type_1_long,
                            'type_2' = type_2_long,
                            'time3' =  time_3_long,
                            'time10' = time_10_long,
                            'time17' = time_17_long,
                            'time24' = time_24_long,
                            'time31' = time_31_long,
                            'N' = 167 ),
                parameters,
                n.chains=4,
                n.iter=5000,
                n.burnin=1000,
                inits=my_inits,
                n.thin=1,
                DIC=T)






parameters <- c("alpha0", "alpha1","alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8","sigma","p","b")
samples <- jags(model = mouseModel ,
                data = list('R_long' = R_long,
                            'N_long' = N_long,
                            'period' = per_long,
                            'mouse'  = mouse_long,
                            'gender' = gender_long,
                            'type_1' = type_1_long,
                            'type_2' = type_2_long,
                            'time3' =  time_3_long,
                            'time10' = time_10_long,
                            'time17' = time_17_long,
                            'time24' = time_24_long,
                            'time31' = time_31_long,
                            'N' = 167 ),
                parameters,
                n.chains=4,
                n.iter=5000,
                n.burnin=1000,
                inits=my_inits,
                n.thin=1,
                DIC=T)



#samples$BUGS$summ
# not sure how this samples from the 4 MC's. I think I am 
#   probably just sampling from the first one this way.
alpha0_sample <- samples$BUGSoutput$sims.list$alpha0
alpha1_sample <- samples$BUGSoutput$sims.list$alpha1
alpha2_sample <- samples$BUGSoutput$sims.list$alpha2
alpha3_sample <- samples$BUGSoutput$sims.list$alpha3
alpha4_sample <- samples$BUGSoutput$sims.list$alpha4
alpha5_sample <- samples$BUGSoutput$sims.list$alpha5
alpha6_sample <- samples$BUGSoutput$sims.list$alpha6
alpha7_sample <- samples$BUGSoutput$sims.list$alpha7
alpha8_sample <- samples$BUGSoutput$sims.list$alpha8
sigma_sample <- samples$BUGSoutput$sims.list$sigma
p_sample <- samples$BUGSoutput$sims.list$p
b_sample <- samples$BUGSoutput$sims.list$b




# show posterior for p, and actual p

par(mfrow=c(3,2))
for (mouse in 1:12) {

  for (i in 1:6) {
    main = str_c('p for mouse ',mouse,' period ',i)
    {
    hist(p_sample[,mouse,i],30,probability=T,main=main,xlim=c(0,1))
    abline(v=P_matrix[mouse,i],col='red')
    }
  }
}

## show random effect
par(mfrow=c(3,2))
for (mouse in 1:12) {
  
  for (i in 1:6) {
    main = str_c('random effect mouse ',mouse,' period ',i)
    {
      hist(b_sample[,mouse,i],30,probability=T,main=main,xlim=c(0,1))
      
    }
  }
}


hist(alpha0_sample,30,probability=T)
hist(alpha1_sample,30,probability=T)
hist(alpha2_sample,30,probability=T)
hist(alpha3_sample,30,probability=T)
hist(alpha4_sample,30,probability=T)
hist(alpha5_sample,30,probability=T)
hist(alpha6_sample,30,probability=T)
hist(alpha7_sample,30,probability=T)
hist(alpha8_sample,30,probability=T)
hist(sigma_sample,30,probability=T)

hist(p_sample[,1,1],30,probability=T,main='mouse 1, period 1')
hist(p_sample[,2,1],30,probability=T,main='mouse 2, period 1')


expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}


###
# Compare type_1 with wild type mice (type_1=type_2 == 0) for
# female mice at time 1 (all time_1 ==0)

par(mfrow=c(1,3))
logit_p_wild_type <- alpha0_sample + alpha1_sample * 1 
logit_p_type_1    <- alpha0_sample + alpha1_sample * 1 +  alpha2_sample 
p_wild_type <- expit(logit_p_wild_type)
p_type_1 <- expit(logit_p_type_1)
hist(p_wild_type,30,probability = T)
hist(p_type_1,30,probability=T)
p_diff <- p_wild_type - p_type_1
hist(p_diff,30,probability=T,main="p_wild_type - p_type_1")
lower <- quantile(p_diff,0.025)
upper <- quantile(p_diff,0.975)
segments(x0 =lower,y0=0.1,x1=upper,y1=0.1,col="red")


# Compare plt with wild type mice for
# female mice at time 1 (all time_1 ==0)

logit_p_wild_type <- alpha0_sample + alpha1_sample * 1 
logit_p_plt       <- alpha0_sample + alpha1_sample * 1 +  alpha3_sample 
p_wild_type <- expit(logit_p_wild_type)
p_plt <- expit(logit_p_plt)
hist(p_wild_type,30,probability = T)
hist(p_plt,30,probability=T)
p_diff <- p_wild_type - p_plt
hist(p_diff,30,probability=T,main='p_wt - p_type2')
lower <- quantile(p_diff,0.025)
upper <- quantile(p_diff,0.975)
segments(x0 =lower,y0=0.1,x1=upper,y1=0.1,col="red")


# Compare type_1 with wild type mice for
# male mice at time 1 (all time_1 ==0)

logit_p_wild_type <- alpha0_sample + alpha1_sample * 1 
logit_p_type_1    <- alpha0_sample + alpha1_sample * 1 +  alpha2_sample 
p_wild_type <- expit(logit_p_wild_type)
p_type_1 <- expit(logit_p_plt)
hist(p_wild_type,30,probability = T)
hist(p_type_1,30,probability=T)
p_diff <- p_wild_type - p_type_1
hist(p_diff,30,probability=T,main='p_wt - p_type1')
lower <- quantile(p_diff,0.025)
upper <- quantile(p_diff,0.975)
segments(x0 =lower,y0=0.1,x1=upper,y1=0.1,col="red")


## male wt vs male plt (type_2), period 1
logit_p_wild_type <- alpha0_sample 
logit_p_plt       <- alpha0_sample + alpha3_sample 
p1 <- expit(logit_p_wild_type)
p2 <- expit(logit_p_plt)
hist(p1,30,probability = T)
hist(p2,30,probability=T)
p_diff <- p1 - p2
hist(p_diff,30,probability=T,main='p_wt - p_plt')
lower <- quantile(p_diff,0.025)
upper <- quantile(p_diff,0.975)
segments(x0 =lower,y0=0.1,x1=upper,y1=0.1,col="red")

#### model diagnostics.

# compute D

# this will store D-values
ds <- numeric()
ds2 <- numeric()

# loop through the MC
counter <- 1
for (t in 1000:4000) {
  result <- 0
  result2 <- 0
   
  # loop over the observations
  for (i in 1:length(N_long)) {
      # loop over time period
      # observed value:
      obs <- R_long[i]
      # predicted value for p for this mouse and time period:
      pred_logit_p <-
        alpha0_sample[t] + 
        alpha1_sample[t] * gender_long[i] +
        alpha2_sample[t] * type_1_long[i] +
        alpha3_sample[t] * type_2_long[i] +
        alpha4_sample[t] * time_3_long[i] +
        alpha5_sample[t] * time_10_long[i] +
        alpha6_sample[t] * time_17_long[i] +
        alpha7_sample[t] * time_24_long[i] +
        alpha8_sample[t] * time_31_long[i] +
        b_sample[t,mouse_long[i],per_long[i]]
      pred_p <- expit(pred_logit_p)
      actual_p <- R_long[i]/N_long[i]
      # expected value of the binomial for this mouse and time period
      exp_obs <- N_long[i] * pred_p 
      result <- result + ((obs-exp_obs)^2)/(exp_obs)
      #print(result)
      # simulated from the model
      sim <- rbinom(1,size=N_long[i],prob=pred_p)
      result2 <- result2 + ((sim-exp_obs)^2)/(exp_obs) 
      #print(c(result,result2))
  }
  ds[counter] <- result
  ds2[counter] <- result2
  counter <- counter + 1
  
  
}


hist(ds)
hist(ds2)


## the moment of truth
plot(ds,ds2,
     xlim=c(min(ds,ds2),max(ds,ds2)),
     ylim=c(min(ds,ds2),max(ds,ds2)))

abline(0,1,col='red')
sum(ds2<ds)
p.value <- sum(ds2<ds)/length(ds)
print("==============")
p.value


#######






##########
#plot(alpha1_sample,pch='.')
traceplot(samples$BUGSoutput,ask=T)

turnout.fit <- samples
#turnout.mcmc <- as.mcmc(turnout.fit)
#turnout.mat <- as.matrix(turnout.mcmc)
#turnout.out <- as.data.frame(turnout.mat)
#p <- turnout.out[, grep("p[", colnames(turnout.out), fixed = T)]


############ GLM ###########

myglm <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long, family="quasibinomial" )

myglm2 <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long
             + gender_long * type_1_long + gender_long * type_2_long 
             , 
             family="quasibinomial" )


summary(myglm)

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

pred <- predict(myglm,se=TRUE)

max_plot = 167
plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  low <- expit(pred$fit[i] - 2 * pred$se.fit[i])
  upp <- expit(pred$fit[i] + 2 * pred$se.fit[i])
  segments(x0=i,y0=low,x1=i,y1=upp,col='black')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
}

library(lme4)
myglm3 <- glmer(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long + (1 | mouse_long), family="binomial")
summary(myglm3)
pre3 <- predict(myglm3)
max_plot = 167
plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  points(i,expit(pre3[i]),col='black',pch='x')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
  segments(i,R_long[i]/N_long[i],i,expit(pre3[i]))
}




library(lme4)
myglm3c <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")
myglm3d <- glm(cbind(R_long,N_long) ~ 1 + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")

summary(myglm3c)
summary(myglm3d)
#(763946-8*2)/5106 + 2 * (8+1)
#1] 167.6142
#> (754808-9*2)/5106 + 2 * (9+1)
#[1] 167.8241

myglm3a <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
                + time_17_long + time_24_long + time_31_long, family="quasibinomial")
myglm3b <- glm(cbind(R_long,N_long) ~ 1 + type_1_long + type_2_long + time_3_long + time_10_long
              + time_17_long + time_24_long + time_31_long, family="quasibinomial")


anova(myglm3a,myglm3b,test='F')


s <- rstudent(myglm3)
qqnorm(s) 
abline(0,5601^0.5,col='red')
summary(myglm3)
hist(s,probability = T,30)
plot(density(s))
### david

fits<-fitted(myglm)

plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  low <- qbinom(0.025,N_long[i],fits[i])/N_long[i]
  upp <- qbinom(0.975,N_long[i],fits[i])/N_long[i]
  segments(x0=i,y0=low,x1=i,y1=upp,col='black')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
}


### beta binomial Bayes


