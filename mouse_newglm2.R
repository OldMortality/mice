#
#
# 2019. Took this from mice.df.rmd.
#      
#
rm(list = ls())
library(stringr)
library(lme4)
library(nlme)
library(gridExtra)
library(ggplot2)
library(gridExtra)

setwd("~/Documents/mice")
source('createMiceDF.R')
source('waldInterval.R')

mice = readFile(filename='data/mice.csv') 
celltypes <- getcelltypes(mice)

# get data in long format
C = 1
celltypes[C]
# z-value for 95% CI
Z =  1.96 # 2.58 #1.96



doCelltype <- function(cellType,mice,Z, returnType='plot') {
  mice.df <- createMiceDF(celltype = cellType, mice )
   
  ## linear model on the logit scale
  
  library(nlme) # will get p-values for anova
  m1 <- lme(logity ~ type * factor(period) ,random=~1|mouse.name, 
            data=mice.df)
  f <- anova(m1)
  f
  summary(m1)
  
  # without interaction
  m0a <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
            data=mice.df)
  f0a <- anova(m0a)
  summary(m0a)
  exp(coefficients(m0a))
  # the other way around
  m0b <- lme(logity ~  factor(period) + type ,random=~1|mouse.name, 
            data=mice.df)
  f0b <- anova(m0b)
  
  # use lme4 for CI's
  m.logit <- lmer(logity ~ type * factor(period) + (1|mouse.name), 
                  data=mice.df)
  w1 <- waldInterval(model=m.logit,Z=Z,FUN=exp)
  
  
  
  
  new.mouse.wt <- data.frame(type=rep("aWT",5),period=c(1,2,3,4,5))
  new.mouse.t2 <- data.frame(type=rep("T2",5),period=c(1,2,3,4,5))
  
  #predict(m.logit,newdata=new.mouse.wt,re.form=~0) %>% invlogit
  #predict(m.logit,newdata=new.mouse.t2,re.form=~0) %>% invlogit
  
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
  
  
  p.wt <- predict(m.logit,newdata=new.mouse.wt,re.form=~0)
  p.wt.low <- invlogit( predict(m.logit,newdata=new.mouse.wt,re.form=~0) - Z * se[1:5]  )
  p.wt.upp <- invlogit( predict(m.logit,newdata=new.mouse.wt,re.form=~0) + Z * se[1:5] )
  p.t2.low <- invlogit( predict(m.logit,newdata=new.mouse.t2,re.form=~0) - Z * se[6:10])
  p.t2.upp <- invlogit( predict(m.logit,newdata=new.mouse.t2,re.form=~0) + Z * se[6:10])
  
  
  t.obs <- c(0,3,10,17,24)
  par(mfrow=c(1,1))
  plot('',xlim=c(0,24),ylim=c(0,0.65), 
       main=paste(celltypes[C], ' mean proportion'),
       xlab='day',
       xaxt='n')
  axis(1,t.obs)
  lines(t.obs,invlogit(predict(m.logit,newdata=new.mouse.wt,re.form=~0)))
  lines(t.obs,invlogit(predict(m.logit,newdata=new.mouse.t2,re.form=~0)),col='red')
  segments(t.obs,p.wt.low,t.obs,p.wt.upp)
  segments(t.obs+0.05,p.t2.low,t.obs+0.05,p.t2.upp,col='red')
  
  df.plot <- data.frame(t.obs = t.obs,
                        wt = invlogit(predict(m.logit,newdata=new.mouse.wt,re.form=~0)),
                        t2 = invlogit(predict(m.logit,newdata=new.mouse.t2,re.form=~0))         
                        )
  
  #breaks=seq(0,0.5,0.1)
  p <- 
    ggplot(data=df.plot) + geom_line(aes(x=t.obs,y=wt)) +
         geom_line(aes(x=t.obs,y=t2,colour='red')) + 
         geom_errorbar(aes(x=t.obs,     ymin=p.wt.low,ymax=p.wt.upp,alpha=0.5)) + 
         geom_errorbar(aes(x=t.obs,ymin=p.t2.low,ymax=p.t2.upp),alpha=0.5,colour='red') +
         scale_x_continuous(breaks=t.obs) +
         scale_y_continuous(limits=c(0,max(p.wt.upp,p.t2.upp)*1.1)) +
         ggtitle(cellType) + 
         ylab('proportion') +
         xlab('day') + 
         theme(legend.position="none") + 
         theme(plot.title = element_text(size=12))
  #c(0,max(df.plot$wt,df.plot$t2)+0.1) +
  
  mice.df$typetime <- factor(paste(mice.df$type,mice.df$period,sep='-'))
  
  reflevs <- levels(mice.df$typetime)[1:5]
  for (i in 1:5) {
    mice.df$typetime <- relevel(mice.df$typetime,
                                ref=reflevs[i])
    m.logit <- lme(logity ~ typetime ,random = ~1|mouse.name, 
              data=mice.df)
    
    w <- waldInterval.nlme(model=m.logit,Z=Z,FUN=exp)
    if (i==1) {
      results <- data.frame(w[4+i,] )
    } else {
      results <- rbind(results,w[4+i,])
    }
    #print(w[4+i,])
  }
  
  stars <- data.frame(
    x <- t.obs[which(results$signif=="*")],
    y <- rep(0,length(x))
  )
  p <- p + geom_text(data=stars,aes(x=x,y=y),label='*')
  p
  # for (i in 1:5) {
  #   if (results$signif[i] == "*") {
  #     p <- p + geom_point(aes(x=x[t.obs[i],y=0),shape='*',col='red')
  #     
  #   }
  # }
  if (returnType=='plot') {
    return(p)
  } else {
    anova.sum <- round(c(f[4,4],f0a[3,4],f0b[3,4]),4)
    return(list(w=w1,f=f,f0a=f0a,f0b=f0b,r=results,a=anova.sum))
  }
  
}

C = 1
par(mfrow=c(3,2))
pvals <- vector()
cells <- vector()
for (C in 1:length(celltypes)) {
  print(C)
  l <- doCelltype(cellType=celltypes[C],mice=mice,Z=Z,returnType='data')
  pvals <- c(pvals,l$r$p.value)
  cells <- c(cells,rep(C,length(l$r$p.value)))
}
cells
df.p <- data.frame(pvals=pvals,cells=cells)
# p-values inm ascending order
df.p <- df.p[with(df.p,order(pvals)), ]
df.p$rank <- seq(1:length(df.p$pvals))

FDR = 0.10
df.p$th <- FDR * df.p$rank / 55
df.p$pvals < df.p$th
length(which(df.p$pvals < df.p$th))
length(which(df.p$pvals<0.05))
df.p$reject <- FALSE 
# the largest p-value which is smaller than the threshold
cutoff <- max(which(df.p$pvals < df.p$th))
# reject H0 for all p-values not larger than cutoff
df.p[1:cutoff,"reject"] <- TRUE
df.p


plot('',xlim=c(0,55),ylim=c(0,1))
for (i in 1:55) {
  if (df.p$reject[i]) {
    col <- 'red'
  } else {
    col <- 'black'
  }
  points(i,df.p$pvals[i],col=col)
}




C = 1
l <- doCelltype(cellType=celltypes[C],mice=mice,Z=1.96,returnType='data')



par(mfrow=c(3,2))
library(purrr)
r <- celltypes[1:2] %>% map(doCelltype,mice=mice,Z=Z) 

returnType='plot'
{
p1 <-  doCelltype(celltypes[1],  mice = mice, Z = Z, returnType = returnType) 
p2 <-  doCelltype(celltypes[2],  mice = mice, Z = Z, returnType = returnType) 
p3 <-  doCelltype(celltypes[3],  mice = mice, Z = Z, returnType = returnType) 
p4 <-  doCelltype(celltypes[4],  mice = mice, Z = Z, returnType = returnType) 
p5 <-  doCelltype(celltypes[5],  mice = mice, Z = Z, returnType = returnType) 
p6 <-  doCelltype(celltypes[6],  mice = mice, Z = Z, returnType = returnType) 
p7 <-  doCelltype(celltypes[7],  mice = mice, Z = Z, returnType = returnType) 
p8 <-  doCelltype(celltypes[8],  mice = mice, Z = Z, returnType = returnType) 
p9 <-  doCelltype(celltypes[9],  mice = mice, Z = Z, returnType = returnType) 
p10 <- doCelltype(celltypes[10], mice = mice,Z = Z,returnType = returnType) 
p11 <- doCelltype(celltypes[11], mice = mice,Z = Z, returnType = returnType) 
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
grid.arrange(p7,p8,p9,p10,p11,ncol=2)
 }


celltypes[1]
p1$r$p.value <- round(p1$r$p.value,3)
p1

l1 <- doCelltype(celltypes[1],mice=mice,Z=Z,returnType='plot') 
l1 <- doCelltype(celltypes[1],mice=mice,Z=Z,returnType='data') 

l1$a




p1$a
p2$a
p3$a
p4$a
p5$a
p6$a
p7$a
p8$a
p9$a
p10$a
p11$a


# create dataset in the long format for the given cell type





# test
# mice.df$typetime <- relevel(mice.df$typetime,
#                             ref='aWT-1')
# m.logit <- lmer(logity ~ typetime + (1|mouse.name), 
#                 data=mice.df)
# w <- waldInterval(model=m.logit,Z=Z,FUN=exp)
# w

 

## linear model on the logit scale

library(nlme) # will get p-values for anova
m1 <- lme(logity ~ type * factor(period) ,random=~1|mouse.name, 
          data=mice.df)
f <- anova(m1)
f
summary(m1)

# without interaction
m0a <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
           data=mice.df)
f0a <- anova(m0a)
exp(coefficients(m0a))
# the other way around
m0b <- lme(logity ~  factor(period) + type ,random=~1|mouse.name, 
           data=mice.df)
f0b <- anova(m0b)

# use lme4 for CI's
m.logit <- lmer(logity ~ type * factor(period) + (1|mouse.name), 
                data=mice.df)
w1 <- waldInterval(model=m.logit,Z=Z,FUN=exp)

m.logit2 <- lme(logity ~ type + factor(period) ,random= ~1|mouse.name, 
                data=mice.df)
summary(m.logit2)


### models for each cell type
library(nlme) # will get p-values for anova

celltype = celltypes[1]
celltype
mice.df <- createMiceDF(celltype = cellType, mice )

m1 <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)
summary(m1)$tTable
waldInterval.nlme(m1,Z = Z,FUN = exp)


celltype = celltypes[2]
celltype
mice.df <- createMiceDF(celltype = celltype, mice )

m2 <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)
summary(m2)$tTable
waldInterval.nlme(m2,Z = Z,FUN = exp)

celltype = celltypes[3]
celltype
mice.df <- createMiceDF(celltype = celltype, mice )

m3 <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)
summary(m3)$tTable
waldInterval.nlme(m3,Z = Z,FUN = exp)

celltype = celltypes[5]
celltype
mice.df <- createMiceDF(celltype = celltype, mice )

m5 <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)

summary(m5)$tTable


celltype = celltypes[7]
celltype
mice.df <- createMiceDF(celltype = celltype, mice )

m7 <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)
anova(m7)
summary(m7)$tTable
waldInterval.nlme(m7,Z = Z,FUN = exp)

celltype = celltypes[11]
celltype

celltype <- celltypes[1]
mice.df <- createMiceDF(celltype = celltype, mice )

m8.i <- lme(logity ~ type * factor(period) ,random=~1|mouse.name, 
          data=mice.df)
anova(m8.i)
m8.p <- lme(logity ~ type + factor(period) ,random=~1|mouse.name, 
          data=mice.df)
anova(m8.p)
m8.type <- lme(logity ~ factor(period) + type   ,random=~1|mouse.name, 
          data=mice.df)
anova(m8.type)


mice.df$typetime <- factor(paste(mice.df$type,mice.df$period,sep='-'))
reflevs <- levels(mice.df$typetime)
mice.df$typetime <- relevel(mice.df$typetime,
                            ref=reflevs[9])
m.logit <- lme(logity ~ typetime ,random = ~1|mouse.name, 
               data=mice.df)

w <- waldInterval.nlme(model=m.logit,Z=Z,FUN=exp)
w


## check T with T-1 as reference level, to see whether
##   there is change vs the previous period, for each type.
##   expecting 11 * 2 * 4 = 88 rows
results.t <- NULL
results.t <- data.frame(celltype = character(88),
                      mousetype = character(88),
                      t0 = character(88),
                      t1 = character(88),
                      pval =  numeric(88),
                      stringsAsFactors = F)
counter <- 0
for (celltype in celltypes) {
  mice.df <- createMiceDF(celltype = celltype, mice )
  mice.df$typetime <- factor(paste(mice.df$type,mice.df$period,sep='-'))
  
  for (mousetype in c('aWT','T2')) {
    for (i in 1:4) {
      counter <- counter + 1
      colref <- paste(type,i,sep='-')
      nextref <- paste('typetime',mousetype,'-',i+1,sep='')
      mice.df$typetime <- relevel(mice.df$typetime,
                                  ref=colref)
      m.logit <- lme(logity ~ typetime ,random = ~1|mouse.name, 
                     data=mice.df)
      w <- waldInterval.nlme(model=m.logit,Z=Z,FUN=exp)
      
      results.t[counter,] <- c(toString(celltype),
                         toString(mousetype),
                         toString(colref),
                         toString(nextref),
                         w[nextref,]$p.value)
    }
  }
}
dim(results.t)
results.t$pval <- round(as.numeric(as.character(results.t$pval)),3)
head(results.t)






anova(m8)
summary(m8)$tTable
waldInterval.nlme(m8,Z = Z,FUN = exp)


mice = readFile(filename='data/mice.csv') 
celltypes <- getcelltypes(mice)

# get data in long format
C = 1
celltypes[C]
# z-value for 95% CI
Z =  1.96 # 2.58 #1.96



cellType = celltypes[1]
ts <- vector()
mice.df <- createMiceDF(celltype = cellType, mice )
for (i in 1:5) {
  t.df <- mice.df[which(mice.df$period==i),c('y','type')]
  t <- t.test(y~type,data=t.df)
  ts[i] <- t$p.value
  }
round(ts,3)
