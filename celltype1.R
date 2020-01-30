mice = readFile(filename='data/mice.csv') 
celltypes <- getcelltypes(mice)
Z =  1.96 # 2.58 #1.96

###
### celltype 1
###
cellType = celltypes[1]
mice.df <- createMiceDF(celltype = cellType, mice )
  
m1 <- lme(logity ~ type  ,random=~1|mouse.name,  
          data=mice.df)
m2 <- lme(logity ~ factor(period) + type,random=~1|mouse.name, 
           data=mice.df)
m2b <- lme(logity ~ type + factor(period) ,random=~1|mouse.name,
          data=mice.df)
anova(m2)
anova(m2b)
anova(m1)

library(lme4)


m1 <- lmer(logity ~ type  + (1|mouse.name), REML=F ,
          data=mice.df)
AIC(m1)
BIC(m2)
m2 <- lmer(logity ~ type  + factor(period) + (1|mouse.name),  REML=F,
           data=mice.df)
AIC(m2)
BIC(m2)

summary(m1)

## policeplot

m.lm <- lmer(y ~ type +  (1|mouse.name) , data=mice.df)
waldInterval(model=m.lm,Z,FUN=identity)
r <- residuals(m.lm)



y.sim <- simulate(m.lm,nsim=10)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=10)
for (i in 1:10) { 
  r <- residuals(lmer( unlist(y.sim[i]) ~ type  +  (1|mouse.name) , data=mice.df))  
  rsim[,i] <- r
}

policePlot(residuals(m.lm),rsim,xl= c(-0.1,0.1))

y.sim <- simulate(m1,nsim=100)
rsim <- matrix(nrow = dim(y.sim)[1],ncol=100)
aics <- vector()
for (i in 1:100) { 
  aics[i] <- AIC(lmer( unlist(y.sim[i]) ~ type  + (1|mouse.name) ,REML=F, data=mice.df))  
  
}
par(mfrow=c(1,1))
hist(aics)
abline(v=AIC(m2),col='red')
abline(v=AIC(m1),col='blue')

###
### celltype 2
###
cellType = celltypes[2]
mice.df <- createMiceDF(celltype = cellType, mice )

m1 <- lme(logity ~ type  ,random=~1|mouse.name,  
          data=mice.df)
m2a <- lme(logity ~ type + factor(period) ,random=~1|mouse.name,
           data=mice.df)

m2b <- lme(logity ~ factor(period) + type,random=~1|mouse.name, 
          data=mice.df)
anova(m1)
anova(m2a)
anova(m2b)
### seems to show that both time and type are relevant.

m <- lmer(logity ~ factor(period) + type + (1|mouse.name), 
           data=mice.df)



m1 <- lmer(logity ~ type  + factor(period) + (1|mouse.name), 
           data=mice.df,REML=F)

new.mouse.wt <- data.frame(type=rep("aWT",5),period=c(1,2,3,4,5))
new.mouse.t2 <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5))

#predict(m.logit,newdata=new.mouse.wt,re.form=~0) %>% invlogit
#predict(m.logit,newdata=new.mouse.t2,re.form=~0) %>% invlogit

p <- doCelltype(cellType,mice=mice,Z=Z,returnType='plot',interaction=F)
p
d <- doCelltype(cellType,mice=mice,Z=Z,returnType='data',interaction=F)
round(d$wt.low,2)
round(d$wt.upp,2)
round(d$wt.est,2)
round(d$t2.est,2)

summary(m1)
w <- waldInterval(m1,Z,exp)
w


###
### celltype 3
###
cellType = celltypes[3]
mice.df <- createMiceDF(celltype = cellType, mice )

m1 <- lme(logity ~ type  ,random=~1|mouse.name,  
          data=mice.df)
m2 <- lme(logity ~ factor(period) + type,random=~1|mouse.name, 
          data=mice.df)
m2b <- lme(logity ~ type + factor(period) ,random=~1|mouse.name,
           data=mice.df)
anova(m2)
anova(m2b)
anova(m1)

library(lme4)


m1 <- lmer(logity ~ type  + (1|mouse.name), REML=F ,
           data=mice.df)
AIC(m1)
BIC(m1)
m2 <- lmer(logity ~ type  + factor(period) + (1|mouse.name),  REML=F,
           data=mice.df)
AIC(m2)
BIC(m2)
summary(m1)


m <- lmer(logity ~ type   + (1|mouse.name),  REML=F,
           data=mice.df)

##
## celltype 4
##

cellType = celltypes[4]
mice.df <- createMiceDF(celltype = cellType, mice )

m1 <- lme(logity ~ type  ,random=~1|mouse.name,  
          data=mice.df)
m2 <- lme(logity ~ factor(period) + type,random=~1|mouse.name, 
          data=mice.df)
m2b <- lme(logity ~ type + factor(period) ,random=~1|mouse.name,
           data=mice.df)
m3 <- lme(logity ~ type * factor(period) ,random=~1|mouse.name,
           data=mice.df)

anova(m3)
anova(m2b)
anova(m1)

library(lme4)


m1 <- lmer(logity ~ type  + (1|mouse.name), REML=F ,
           data=mice.df)
AIC(m1)
BIC(m2)
m2 <- lmer(logity ~ type  + factor(period) + (1|mouse.name),  REML=F,
           data=mice.df)
AIC(m2)
BIC(m2)

p <- doCelltype(cellType,mice=mice,Z=Z,returnType='plot',interaction=T)
p
d <- doCelltype(cellType,mice=mice,Z=Z,returnType='data',interaction=T)
round(d$wt.low,2)
round(d$wt.upp,2)
round(d$wt.est,2)
round(d$t2.est,2)



##
## celltype 5
##

p <- doCelltype(cellType,mice=mice,Z=Z,returnType='plot',interaction=T)

cellType = celltypes[5]
mice.df <- createMiceDF(celltype = cellType, mice )

m1 <- lme(logity ~ type  ,random=~1|mouse.name,  
          data=mice.df)
m2 <- lme(logity ~ factor(period) + type,random=~1|mouse.name, 
          data=mice.df)
m2b <- lme(logity ~ type + factor(period) ,random=~1|mouse.name,
           data=mice.df)
m3 <- lme(logity ~ type * factor(period) ,random=~1|mouse.name,
          data=mice.df)

anova(m3)
anova(m2)
anova(m2b)
anova(m1)

library(lme4)

m2 <- lmer(logity ~ factor(period) + (1|mouse.name),  REML=F,
           data=mice.df)


AIC(m2)
AIC(m3)

BIC(m2)

p <- doCelltype(cellType,mice=mice,Z=Z,returnType='plot',interaction=T)
p
d <- doCelltype(cellType,mice=mice,Z=Z,returnType='data',interaction=T)
round(d$wt.low,2)
round(d$wt.upp,2)
round(d$wt.est,2)
round(d$t2.est,2)





