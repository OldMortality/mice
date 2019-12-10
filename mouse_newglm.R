#
#
# 2019. Took this from mice.df.rmd.
#      
#
rm(list = ls())
library(stringr)
library(lme4)
setwd("~/Documents/mice")
source('createMiceDF.R')
source('waldInterval.R')

mice = readFile(filename='data/mice.csv') 
celltypes <- getcelltypes(mice)

# get data in long format
C = 5
celltypes[C]
# create dataset in the long format for the given cell type
mice.df <- createMiceDF(celltype = celltypes[C], mice )
head(mice.df)

mice=NULL

dim(mice.df)
# z-value for 95% CI
Z = 1.96
  


## linear model
m.lm <- lmer(y ~ type * factor(period) +  (1|mouse.name) , data=mice.df)
waldInterval(model=m.lm,Z,FUN=identity)
r <- residuals(m.lm)
timesResiduals(residuals(m.lm),df=mice.df,main='lm')



y.sim <- simulate(m.lm,nsim=10)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=10)
for (i in 1:10) { 
    r <- residuals(lmer( unlist(y.sim[i]) ~ type * factor(period) +  (1|mouse.name) , data=mice.df))  
    rsim[,i] <- r
}


policePlot(residuals(m.lm),rsim,xl= c(-0.1,0.1))
plots(df=mice.df,fit=predict(m.lm),r=r <- residuals(m.lm),qq=T)


# show fitted vs actual
f <- fitted(m.lm)
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.4),main='fitted vs observed (red) Normal model')
grid()
points(1:38,mice.df$y,col='red')


###
### Beta regression
###
library(glmmTMB)

m.beta <- glmmTMB(y~type * factor(period)  + (1|mouse.name), 
                  data=mice.df, family=list(family="beta",link="logit"))

s <- summary(m.beta)
# exp to get ORs
waldIntervalBeta(model=m.beta,Z=Z,FUN=exp)

timesResiduals(residuals(m.beta),df=mice.df,main='beta')


y.sim <- simulate(m.beta,nsim=10)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=10)
for (i in 1:10) { 
  r <- residuals(glmmTMB( unlist(y.sim[i]) ~ type * factor(period)  + (1|mouse.name), 
                         data=mice.df, family=list(family="beta",link="logit")))
  rsim[,i] <- r
}

policePlot(residuals(m.beta),rsim,xl= c(-0.1,0.1))
plots(df=mice.df,fit=predict(m.beta),r=r <- residuals(m.beta),qq=F)

 
# show fitted vs actual
f.beta <- fitted(m.beta)
par(mfrow=c(1,1))
plot(1:38,f.beta,ylim=c(0,.4),main='fitted vs observed(red), beta-regression')
grid()
points(1:38,mice.df$y,col='red')
plots(df=mice.df,fit=predict(m.beta),r=residuals(m.beta))

par(mfrow=c(1,1))
plot(1:38,f.beta,ylim=c(0,.4),main='fitted lm vs beta regression(red)')
points(1:38,fitted(m.lm),col='red')

##
## linear model on the logit scale
##
##
m.lm <- lmer(logity ~ type * factor(period) +  (1|mouse.name) , data=mice.df)
waldInterval(model=m.lm,Z=Z,FUN=exp)
r <- residuals(m.lm)

#boxplot(r~mice.df$type)
#boxplot(r~mice.df$period)



y.sim <- simulate(m.lm,nsim=10)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=10)
for (i in 1:10) { 
  r <- residuals(lmer( unlist(y.sim[i]) ~ type * factor(period) +  (1|mouse.name) , data=mice.df))  
  rsim[,i] <- r
}


policePlot(residuals(m.lm),rsim,xl= c(-1,1))
plots(df=mice.df,fit=invlogit(predict(m.lm)),r=r <- residuals(m.lm),qq=T)


# show fitted for all trhe models
f <- invlogit(fitted(m.lm))
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.4),main='fitted vs observed')
grid()
points(1:38,mice.df$y,col='red')
points(1:38,fitted(m.beta),col='blue')
legend(20,.4, 
       c('logit lm,','beta','lm'),
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=c("black","blue","red"))


timesResiduals(residuals(m.beta),df=mice.df,main='logit lm')



new.mouse <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5),mouse.name='mouse8')

new.mouse.wt <- data.frame(type=rep("aWT",5),period=c(1,2,3,4,5))
new.mouse.t2 <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5))


predict(m.lm,newdata=new.mouse.wt,re.form=~0)


#predict(m.lm,newdata=new.mouse.wt,re.form=~0)
par(mfrow=c(1,1))
plot('',xlim=c(1,5),ylim=c(0,0.06))
lines(1:5,invlogit(predict(m.lm,newdata=new.mouse.wt,re.form=~0)))
lines(1:5,invlogit(predict(m.lm,newdata=new.mouse.t2,re.form=~0)),col='red')
