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
C = 7
celltypes[C]
# create dataset in the long format for the given cell type
mice.df <- createMiceDF(celltype = celltypes[C], mice )
head(mice.df)

mice=NULL

dim(mice.df)
# z-value for 95% CI
Z = 1.96
  

### nonlinear in x, weeks numeric
mice.df$t <- as.numeric(as.character(mice.df$period))
mice.df[which(mice.df$t==5),'t'] <- 24
mice.df[which(mice.df$t==4),'t'] <- 17
mice.df[which(mice.df$t==3),'t'] <- 10
mice.df[which(mice.df$t==2),'t'] <- 3
mice.df[which(mice.df$t==1),'t'] <- 0



mice.df$tsq <- mice.df$t ^ 2
table(mice.df$t,mice.df$period)
library(nlme)
m.nl2 <- lme(logity ~ type * (tsq + t)  ,random=~1|mouse.name, 
                data=mice.df)
summary(m.nl2)
Anova(m.nl2)

##
## linear model on the logit scale
##
##
m.logit <- lmer(logity ~ type * factor(period) +  (1|mouse.name) , 
                data=mice.df)

library(nlme) # will get p-values
mm <- lme(logity ~ type * factor(period) ,random=~1|mouse.name, 
          data=mice.df)
  
anova(mm)

summary(m.logit)
anova(m.logit)
m.logit2 <- lmer(logity ~ type + factor(period) +  (1|mouse.name) , 
                data=mice.df)
#AIC(m.logit)
anova(m.logit,m.logit2)

m1 <- lmer(logity ~ type + factor(period) +  (1|mouse.name) , 
                 data=mice.df)
m2 <- lmer(logity ~ type + (1|mouse.name) , 
                 data=mice.df)
anova(m1,m2)

m1 <- lmer(logity ~ factor(period) +  (1|mouse.name) , 
           data=mice.df)
m2 <- lmer(logity ~ factor(period) + type + (1|mouse.name) , 
           data=mice.df)
anova(m1,m2)



anova(m.logit)


# library(car)
#Anova(m.logit)


w <- waldInterval(model=m.logit,Z=Z,FUN=exp)
w
r <- residuals(m.logit)



y.sim <- simulate(m.logit,nsim=10)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=10)
for (i in 1:10) { 
  r <- residuals(lmer( unlist(y.sim[i]) ~ type * factor(period) +  (1|mouse.name) , data=mice.df))  
  rsim[,i] <- r
}


policePlot(residuals(m.logit),rsim,xl= c(-1,1))
plots(df=mice.df,fit=invlogit(predict(m.logit)),r=r <- residuals(m.logit),qq=T)


new.mouse <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5),mouse.name='mouse8')

new.mouse.wt <- data.frame(type=rep("aWT",5),period=c(1,2,3,4,5))
new.mouse.t2 <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5))
predict(m.logit,newdata=new.mouse.wt,re.form=~0) %>% invlogit
predict(m.logit,newdata=new.mouse.t2,re.form=~0) %>% invlogit

# WT period 0
V <- vcov(m.logit)
x1 <- matrix(c(1,0,0,0,0,0,0,0,0,0,
               1,0,1,0,0,0,0,0,0,0,
               1,0,0,1,0,0,0,0,0,0,
               1,0,0,0,1,0,0,0,0,0,
               1,0,0,0,0,1,0,0,0,0,
               1,1,0,0,0,0,0,0,0,0,
               1,1,1,0,0,0,1,0,0,0,
               1,1,0,1,0,0,0,1,0,0,
               1,1,0,0,1,0,0,0,1,0,
               1,1,0,0,0,1,0,0,0,1
               ),
             
             nrow=10,ncol=10,byrow=T)

se <- sqrt(diag(x1 %*% V %*% t(x1)))
 
# x1 <- matrix(c(1,0,1,0,0,0,0,0,0,0),nrow=1)
# se <- sqrt(diag(x1 %*% V %*% t(x1)))

#
p.wt <- predict(m.logit,newdata=new.mouse.wt,re.form=~0)
p.wt.low <- invlogit( predict(m.logit,newdata=new.mouse.wt,re.form=~0) - Z * se[1:5]  )
p.wt.upp <- invlogit( predict(m.logit,newdata=new.mouse.wt,re.form=~0) + Z * se[1:5] )
p.t2.low <- invlogit( predict(m.logit,newdata=new.mouse.t2,re.form=~0) - Z * se[6:10])
p.t2.upp <- invlogit( predict(m.logit,newdata=new.mouse.t2,re.form=~0) + Z * se[6:10])


t.obs <- c(0,3,10,17,24)
par(mfrow=c(1,1))
plot('',xlim=c(0,24),ylim=c(0,0.65), 
     main=paste(celltypes[C], ' mean proportion'),
     xaxt='n')
axis(1,t.obs)
lines(t.obs,invlogit(predict(m.logit,newdata=new.mouse.wt,re.form=~0)))
lines(t.obs,invlogit(predict(m.logit,newdata=new.mouse.t2,re.form=~0)),col='red')
segments(t.obs,p.wt.low,t.obs,p.wt.upp)
segments(t.obs+0.05,p.t2.low,t.obs+0.05,p.t2.upp,col='red')


timesResiduals(residuals(m.logit),df=mice.df,main='lm')
# show fitted vs actual
f <- invlogit(fitted(m.logit))
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.06),main='fitted vs observed (red) Normal model')
grid()
points(1:38,mice.df$y,col='red')



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


# show fitted for all trhe models
f <- invlogit(fitted(m.lm))
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.06),main='fitted vs observed')
grid()
points(1:38,mice.df$y,col='red')
points(1:38,fitted(m.beta),col='blue')
points(1:38,invlogit(fitted(m.logit)),col='black')

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

fw<-data.frame(fertilizer="added",week=6)
p1<-predict(a1,fw,re.form=NA)
head(model.matrix(a1));
x1<-matrix(c(1,0,6,0),1,4)
V1<-vcov(a1);
se1<-sqrt(diag(x1%*%V1%*%t(x1)))
#[and likewise for a2 and a3...]
p<-c(p1,p2,p3);se<-c(se1,se2,se3)
z<-qnorm(0.975);low<-p-z*se;upp<-p+z*se
res<-rbind(p,se,low,upp);print(round(res,digits=3))


p1 <- predict(m.logit,newdata=new.mouse.wt,re.form=~0)
V1 <- vcov(m.logit)


for (i in 1:8) {
  print(celltypes[i])
  mice.df <- createMiceDF(celltype = celltypes[i], mice )
  m.logit <- lmer(logity ~ type * factor(period) +  (1|mouse.name) , 
                  data=mice.df)
  #AIC(m.logit)
  waldInterval(model=m.logit,Z=Z,FUN=exp)
}

walds <- celltypes[1:8] %>%
    map(createMiceDF,mice=mice) %>%
    map((function(x) {lmer(logity ~ type * factor(period) +  (1|mouse.name) , 
                     data=x)}) )%>%
    map(waldInterval,Z-Z,FUN=exp)
walds$celltype <- celltypes[1:8]
walds


mice.df$typetime <- factor(paste(mice.df$type,mice.df$period))
m.logit1 <- lmer(logity~type * period + (1|mouse.name),data=mice.df)
m.logit2 <- lmer(logity~type + period + (1|mouse.name),data=mice.df)


m.logit1 <- lmer(logity~type * period + (1|mouse.name) ,data=mice.df)
m.logit2 <- lmer(logity~type + period + (1|mouse.name) ,data=mice.df)

anova(m.logit2,m.logit1)
library(lmtest)
lrtest(m.logit1,m.logit2)


w <- waldInterval(m.logit)

mice.df$typetime <- relevel(mice.df$typetime,
                           ref='aWT 2')

