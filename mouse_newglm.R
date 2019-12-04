#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#


rm(list = ls())
library(stringr)


library(lme4)
# one or the other
#detach(lme4)
#library(glmmTMB)

        
setwd("~/Documents/mice")
source('createMatrix.R')

mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]
mice.s <- read.csv("~/Documents/mice/data/mice2.csv",header=TRUE)[1:40,]
stroke_size <- scale(mice.s$stroke_size)


{
mice$female <- grepl('female',mice$Animal_type)
mice$type_1 <- grepl('PLT',mice$Animal_type)
mice$type_2 <- grepl('Pound',
                mice$Animal_type)
mice$type <- 'WT'
mice[which(mice$type_1),"type"] <- 'PLT'
mice[which(mice$type_2),"type"] <- 'T2'

celltype <- "CD19pos_B220"
# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24)
# these are labeled 1 through to 5
periods <- 5

mice <- mice[-which(mice$type=='PLT'),]
table(mice$type)


# N_0 <- mice$Live_cells_0
# R_0 <- mice$CD19pos_B220pos_0
# N_3 <- mice$Live_cells_3
# R_3 <- mice$CD19pos_B220pos_3
# N_10 <- mice$Live_cells_10
# R_10 <- mice$CD19pos_B220pos_10
# N_17 <- mice$Live_cells_17
# R_17 <- mice$CD19pos_B220pos_17
# N_24 <- mice$Live_cells_24
# R_24 <- mice$CD19pos_B220pos_24
  
female <- grepl('female',mice$Animal_type)
type_1 <- grepl('PLT',mice$Animal_type)
type_2 <- grepl('Pound',mice$Animal_type)

                


celltype <- "CD11cpos"

N_matrix <- getNMatrix(mice,celltype)
R_matrix <- getRMatrix(mice,celltype)



# long notation (i.e. observations as rows)
per_long <- vector()
mouse_long <- vector()
female_long <- vector()
type_1_long <- vector()
type_2_long <- vector()
 
N_long <- vector()
R_long <- vector()
stroke_size_long <- vector()
  
counter = 0
skipped = 0
for (i in 1: N) {
  for (j in 1:5) {
    if (is.na(N_matrix[i,j])) {
      skipped <- skipped + 1
      next
    } else {
      counter <- counter + 1
    }
    # periods run from [1,5] ~ [0,3,10,17,24]
    per_long[counter]    <- j
    mouse_long[counter]  <- str_c('mouse',i,sep='')
    female_long[counter] <- female[i]
    type_1_long[counter] <- type_1[i]
    type_2_long[counter] <- type_2[i]
    stroke_size_long[counter] <- stroke_size[i]
    N_long[counter]      <- N_matrix[i,j]
    R_long[counter]      <- R_matrix[i,j]
  }
}
  
  
  
} 
  mice2 <- data.frame(
    y = R_long / N_long,
    period = factor(per_long),
    mouse = mouse_long,
    female = female_long,
    type1 = type_1_long,
    type2 = type_2_long,
    stroke_size = stroke_size_long
  )
  
  head(mice2)
  mice2$type <- "aWT"
  mice2[which(mice2$type2),"type"] <- "T2"
  
  mice2$type <- factor(mice2$type)
  relevel(mice2$type,ref='WT')
  
  

  head(mice2)
      
  table(mice2$type)

  m <- lmer(y~ type + factor(period) * female +  (1|mouse) , data=mice2)
  s <- summary(m)
  CI_lower <- s$coefficients[,2] - 1.96*s$coefficients[,2]
  CI_upper <- s$coefficients[,2] + 1.96*s$coefficients[,2]
  cbind(CI_lower,CI_upper)
  
  
  
  
  
mice2 <- mice2[-which(mice2$female),]
head(mice2)
dim(mice2)

m <- lmer(y~ type * factor(period) +  (1|mouse) , data=mice2)

s <- summary(m)
CI_lower <- s$coefficients[,1] - 1*s$coefficients[,2]
CI_upper <- s$coefficients[,1] + 1*s$coefficients[,2]
round(cbind(CI_lower,CI_upper),2)

round(confint(m,level=0.9),2)

plot(m)
r<- residuals(m)

boxplot(r~mice2$type)
boxplot(r~mice2$period)

# police method
par(mfrow=c(2,5))
real.pos <- floor(runif(1,min=1,max=11))
for (i in 1:10) {
  y.sim <- unlist(simulate(m))
  m.sim <- lmer(y.sim~type * period + (1|mouse) , data=mice2)
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
plots(d=mice2,p=p,r=r)

qqnorm(residuals(m,scale=T))
abline(0,1,col='red')


# show fitted vs actual
f <- fitted(m)
par(mfrow=c(1,1))
plot(1:38,f,xlim=c(1,38),ylim=c(0,0.4),main='fitted vs observed (red) Normal model')
grid()
points(1:38,mice2$y,col='red')





