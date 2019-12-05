#
#
# 2019. Took this from mice.df.rmd.
#       Fits Bayesian model. 
#
#
rm(list = ls())
library(stringr)
library(lme4)
setwd("~/Documents/mice")
source('createMiceDF.R')

mice.df <- createMiceDF(celltype= "CD19pos_B220")
  
s <- summary(m)

formula <- "y~ type + factor(period) * female +  (1|mouse)" 
m <- lmer( formula , data=mice.df)
s <- summary(m)
CI_lower <- s$coefficients[,2] - 1.96*s$coefficients[,2]
CI_upper <- s$coefficients[,2] + 1.96*s$coefficients[,2]
cbind(CI_lower,CI_upper)
  
  
mice.df <- mice.df[-which(mice.df$female),]
head(mice.df)
dim(mice.df)
m <- lmer(y~   type * factor(period) +  (1|mouse) , data=mice.df)
s <- summary(m)
CI_lower <- s$coefficients[,1] - 1.96*s$coefficients[,2]
CI_upper <- s$coefficients[,1] + 1.96*s$coefficients[,2]
round(cbind(CI_lower,CI_upper),2)

#round(confint(m,level=0.90),2)

#plot(m)
r <- residuals(m)

boxplot(r~mice.df$type)
boxplot(r~mice.df$period)



policePlot()

# police method
par(mfrow=c(2,5))
real.pos <- floor(runif(1,min=1,max=11))
for (i in 1:10) {
  y.sim <- unlist(simulate(m))
  m.sim <- lmer(y.sim~type * period + (1|mouse) , data=mice.df)
  if (i == real.pos) {
    hist(residuals(m),main="",15,xlim=c(-0.1,0.1))
  } else {
    hist(residuals(m.sim),main="",15,xlim=c(-0.1,0.1))
  }
}

#hist(r)
#qqplot(residuals(m,type="ss"))



plots <- function(d,p,r) {
  par(mfrow=c(2,2))
  col <- rep('black',dim(d)[1])
  col[which(d$type=='T2')] <- 'red'
  plot(r~p,main="residuals ~ fitted",col=col)
  abline(h=0,col='red')
#  plot(r~d$stroke_size)
  
  boxplot(r~d$type,main='residuals~type')
  boxplot(r~d$period,main='residuals~time')
}

r <- residuals(m)
p <- predict(m)
plots(d=mice.df,p=p,r=r)

qqnorm(residuals(m,scale=T))
abline(0,1,col='red')


# show fitted vs actual
f <- fitted(m)
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.4),main='fitted vs observed (red) Normal model')
grid()
points(1:38,mice.df$y,col='red')


###
### Beta regression
###
library(glmmTMB)

m.beta <- glmmTMB(y~type * factor(period)  + (1|mouse), 
                  data=mice.df, family=list(family="beta",link="logit"))

s <- summary(m.beta)
CI_lower <- s$coefficients$cond[,1] - 1.67*s$coefficients$cond[,2]
CI_upper <- s$coefficients$cond[,1] + 1.67*s$coefficients$cond[,2]
cbind(exp(CI_lower),exp(CI_upper))



# police method
par(mfrow=c(2,5))
real.pos <- floor(runif(1,min=1,max=11))
for (i in 1:10) {
  y.sim <- unlist(simulate(m.beta))
  m.sim <- lmer(y.sim~type * period + (1|mouse) , data=mice.df)
  if (i == real.pos) {
    hist(residuals(m.beta),main="",15,xlim=c(-0.1,0.1))
  } else {
    hist(residuals(m.sim),main="",15,xlim=c(-0.1,0.1))
  }
}


# show fitted vs actual
f.beta <- fitted(m.beta)
par(mfrow=c(1,1))
plot(1:38,f.beta,ylim=c(0,.4),main='fitted vs observed(red), beta-regression')
grid()
points(1:38,mice.df$y,col='red')
plots(d=mice.df,p=predict(m.beta),r=residuals(m.beta))




