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
#library(glmmTMB)

setwd("~/Documents/mice")

mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]

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
  #mice <- mice[-which(mice$female),]
  #table(mice$type,mice$female)
  
  
  
  
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
  
  female <- grepl('female',mice$Animal_type)
  type_1 <- grepl('PLT',mice$Animal_type)
  type_2 <- grepl('Pound',mice$Animal_type)
  
  
  N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24))
  R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24))
  
  
  # long notation (i.e. observations as rows)
  per_long <- vector()
  mouse_long <- vector()
  female_long <- vector()
  type_1_long <- vector()
  type_2_long <- vector()
  
  N_long <- vector()
  R_long <- vector()
  
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
      N_long[counter]      <- N_matrix[i,j]
      R_long[counter]      <- R_matrix[i,j]
    }
  }
  
  
  
  #
} 
mice2 <- data.frame(
  y = R_long / N_long,
  period = factor(per_long),
  mouse = mouse_long,
  female = female_long,
  type1 = type_1_long,
  type2 = type_2_long
)

head(mice2)
mice2$type <- "WT"
mice2[which(mice2$type2),"type"] <- "T2"

mice2$type1 <- NULL
mice2$type2 <- NULL

head(mice2)


mice2 <- mice2[-which(mice2$female),]
head(mice2)

m <- lmer(y~type * factor(period) +  (1|mouse) , data=mice2,REML=F)
AIC(m)
m <- lmer(y~type + factor(period) +  (1|mouse) , data=mice2,REML=F)
AIC(m)
s <- summary(m)
CI_lower <- s$coefficients[,1] - 1.96*s$coefficients[,2]
CI_upper <- s$coefficients[,1] + 1.96*s$coefficients[,2]
cbind(CI_lower,CI_upper)

confint(m)



plot(m)
r<- residuals(m)
boxplot(r~mice2$type)
boxplot(r~mice2$period)

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

library(betareg)
# does not support mixed models
m <-betareg(y~type * factor(period)  , data=mice2)

library(glmmTMB)
m <- glmmTMB(y~type * factor(period)  + (1|mouse), 
             data=mice2, family=list(family="beta",link="logit"))

s <- summary(m)
CI_lower <- s$coefficients$cond[,2] - 1.96*s$coefficients$cond[,2]
CI_upper <- s$coefficients$cond[,2] + 1.96*s$coefficients$cond[,2]
cbind(exp(CI_lower),exp(CI_upper))
