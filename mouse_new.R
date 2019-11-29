#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#
#library(rjags)
library(R2jags)

# macFile <- "/Users/micheldelange/Documents/mice/data/mice.csv"
# ubuntuFile <- "~/Documents/mice/data/mice.csv"
# 
# if (file.exists(macFile)) {
#   mice <- read.csv(macFile,header=TRUE)[1:40,]
# } else {
#   mice <- read.csv(ubuntuFile,header=TRUE)[1:40,]
# }
# 

setwd("~/Documents/mice")

mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]


source('setinits.R')

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
periods <- 5

createList <- function(mice) {

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
    # periods run from [0,4] ~ [0,3,10,17,24]
    per_long[counter]    <- j-1
    mouse_long[counter]  <- i
    female_long[counter] <- female[i]
    type_1_long[counter] <- type_1[i]
    type_2_long[counter] <- type_2[i]
    N_long[counter]      <- N_matrix[i,j]
    R_long[counter]      <- R_matrix[i,j]
    }
  }
  
  # 
  theList <- list('R' = R_long,
                  'Ncells' = N_long,
                  'period' = per_long,
                  'mouse'  = mouse_long,
                  'female' = female_long,
                  # type as dummy variable
                  'type_1' = type_1_long,
                  'type_2' = type_2_long,
                  # time as dummy variable
                  'time3' =  per_long == 1,
                  'time10' = per_long == 2,
                  'time17' = per_long == 3,
                  'time24' = per_long == 4,
                  'N' = length(R_long))
  return(theList)
}

theList <- createList(mice)

mouseModel <- function()
{
  for ( i in 1 : N) {
    R[i] ~ dbin(p[mouse[i], period[i]],Ncells[i])
    b[mouse[i], period[i]] ~ dnorm(0.0, tau)
    logit(p[mouse[i], period[i]]) <- a0 + 
      a1 * female[i] +
      a2 * type_1[i] +
      a3 * type_2[i] +
      a4 * time3[i] +
      a5 * time10[i] +
      a6 * time17[i] +
      a7 * time24[i] +
      b1 * time3[i] * type_2[i] +
      b2 * time10[i] * type_2[i] +
      b3 * time17[i] * type_2[i] +
      b4 * time24[i] * type_2[i] +
      g1 * time3[i] * type_1[i] +
      g2 * time10[i] * type_1[i] +
      g3 * time17[i] * type_1[i] +
      g4 * time24[i] * type_1[i] +
      b[mouse[i],period[i]]
  }
  
  # priors
  a0 ~ dnorm(0.0,1.0E-6)
  a1 ~ dnorm(0.0,1.0E-6)
  a2 ~ dnorm(0.0,1.0E-6)
  a3 ~ dnorm(0.0,1.0E-6)
  a4 ~ dnorm(0.0,1.0E-6)
  a5 ~ dnorm(0.0,1.0E-6)
  a6 ~ dnorm(0.0,1.0E-6)
  a7 ~ dnorm(0.0,1.0E-6)
  b1 ~ dnorm(0.0,1.0E-6) 
  b2 ~ dnorm(0.0,1.0E-6) 
  b3 ~ dnorm(0.0,1.0E-6) 
  b4 ~ dnorm(0.0,1.0E-6) 
  g1 ~ dnorm(0.0,1.0E-6) 
  g2 ~ dnorm(0.0,1.0E-6) 
  g3 ~ dnorm(0.0,1.0E-6) 
  g4 ~ dnorm(0.0,1.0E-6) 
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}


 

##
## Rename bugs output, so it is a bit more user friendly
## 
renameSamples <- function(samples){
  return(data.frame(samples$BUGSoutput$sims.list))
    # result <- data.frame(
    #   a0 = samples$BUGSoutput$sims.list$a0,
    #   a1 = samples$BUGSoutput$sims.list$a1,
    #   a2 =  samples$BUGSoutput$sims.list$a2,
    #   a3 =  samples$BUGSoutput$sims.list$a3,
    #   a4 =  samples$BUGSoutput$sims.list$a4,
    #   a5 =  samples$BUGSoutput$sims.list$a5,
    #   a6 =  samples$BUGSoutput$sims.list$a6,
    #   a7 =  samples$BUGSoutput$sims.list$a7,
    #   a8 =  samples$BUGSoutput$sims.list$a8,
    #   b1 =  samples$BUGSoutput$sims.list$b1,
    #   b2 =  samples$BUGSoutput$sims.list$b2,
    #   b3 =  samples$BUGSoutput$sims.list$b3,
    #   b4 =  samples$BUGSoutput$sims.list$b4,
    #   g1 =  samples$BUGSoutput$sims.list$g1,
    #   g2 =  samples$BUGSoutput$sims.list$g2,
    #   g3 =  samples$BUGSoutput$sims.list$g3,
    #   g4 =  samples$BUGSoutput$sims.list$g3,
    #   sigma =  samples$BUGSoutput$sims.list$sigma)
      #p =  samples$BUGSoutput$sims.list$p,
      #b = samples$BUGSoutput$sims.list$b)
  return(result)
}

#res <- renameSamples(samples)


parameters <- c("a0", "a1","a2","a3","a4","a5",
                "a6","a7","b1","b2","b3","b4",
                "sigma","p","b")

getSamples <- function(model,data,parameters) {

  samples <- jags(model = mouseModel ,
                data = theList, 
                parameters ,
                n.chains = 4,
                n.iter = 5000,
                n.burnin = 1000,
                inits = setInits(i1 = 0,i2 = 0.1,i3 = -0.1,i4 = 0.245),
                n.thin = 1,
                DIC = T)
   return(renameSamples(samples))
}

samples <- getSamples(model=mouseModel,data=theList,parameters=parameters)



#samples$BUGS$summ
# not sure how this samples from the 4 MC's. I think I am 
#   probably just sampling from the first one this way.


expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}



# get all male WT mice

getMaleWT <- function() {
  return(which(!mice$female))
}




# get linear predictor
#   s is the mcmc samples df
#   period should be [0..4]; type ["WT","T2"]
#   we don't do period 5, ie day 31

getLPFactory <- function(model) {
  if (model=="mouseModel") {
    getLP <- function(s,female, type,period) {
      lp <- s$a0 +  
            s$a1 * female + 
            s$a2 * (type == "T1") +
            s$a3 * (type == "T2") +
        (period == 1) * (s$a4 + s$b1 * (type == "T2") ) +
        (period == 2) * (s$a5 + s$b2 * (type == "T2") ) +
        (period == 3) * (s$a6 + s$b3 * (type == "T2") ) +
        (period == 4) * (s$a7 + s$b4 * (type == "T2") )
      return(lp)
    }
  }
  if (model == "mouseModel2") {
    print('returning LP for mousemodel2')
    ## no time:type interactions, ie. no beta
    getLP <- function(s,type,period) {
      lp <- s$a0 +  s$a1 +  s$a3 * (type == "T2") +
        (period == 1) * (s$a4  ) +
        (period == 2) * (s$a5  ) +
        (period == 3) * (s$a6  ) +
        (period == 4) * (s$a7  )
      return(lp)
    }
  }
  if (model == "mouseModel3") {
    ## no time main effect
    getLP <- function(s,type,period) {
      lp <- s$a0 +  s$a1 +  s$a3 * (type == "T2") 
      return(lp)
    }
  }
  return(getLP)
  
  
}

#getLP <- getLPFactory("mouseModel3")
 


# get Odds Ratio posterior distribution for male T2 vs male WT at
#   the given period. Works for period on [0,4] ~ [day0,day24]
getORs <- function(samples,period,model) {
  getLP <- getLP <- getLPFactory(model) 
  return(exp(getLP(samples,type = "T2", period = period)) /
           exp(getLP(samples,type = "WT", period = period)))
}

# Odds ratio for given type (T2 or WT), period [0.4]
#      compared to the same time, period 0
getORsTime <- function(samples,type,period,model) {
  getLP <- getLP <- getLPFactory(model) 
  return(exp(getLP(samples,type = type, period = period)) /
           exp(getLP(samples,type = type, period = 0)))
}


par(mfrow=c(1,1))

#period should be [0..4] [0:0 ; 1:3 ; 2:10 ; 3:17; 4:24]

{
par(mfrow=c(2,3))

    
  OR <- getORs(samples,period = 0,"mouseModel")
  p <- hist(OR,main=paste("OR T2/WT ", celltype, "day=0"))
  abline(v=quantile(OR,0.025),col='blue')
  abline(v=quantile(OR,0.975),col='blue')
  abline(v=1,col='red')
  
   
  OR <- getORs(samples,period = 1,"mouseModel")
  p <- hist(OR,main=paste("OR T2/WT " , celltype, "day=3"))
  abline(v=quantile(OR,0.025),col='blue')
  abline(v=quantile(OR,0.975),col='blue')
  abline(v=1,col='red')
  
   
  OR <- getORs(samples,period = 2,"mouseModel")
  p <- hist(OR,main=paste("OR T2/WT ", celltype, "day=10"))
  abline(v=quantile(OR,0.025),col='blue')
  abline(v=quantile(OR,0.975),col='blue')
  abline(v=1,col='red')
  
  
  OR <- getORs(samples,period = 3,"mouseModel")
  p <- hist(OR,main=paste("OR T2/WT ", celltype, "day=17"))
  abline(v=quantile(OR,0.025),col='blue')
  abline(v=quantile(OR,0.975),col='blue')
  abline(v=1,col='red')
  
 
  OR <- getORs(samples,period = 4,"mouseModel")
  p <- hist(OR,main=paste("OR T2/WT ", celltype, "day=24"))
  abline(v=quantile(OR,0.025),col='blue')
  abline(v=quantile(OR,0.975),col='blue')
  abline(v=1,col='red')
  
}  
  
{
par(mfrow=c(2,2))  
hist(b1_sample,main="b1")
abline(v=quantile(b1_sample,0.025),col='blue')
abline(v=quantile(b1_sample,0.975),col='blue')
abline(v=0,col='red')

hist(b2_sample,main="b2")
abline(v=0,col='red')
abline(v=quantile(b2_sample,0.025),col='blue')
abline(v=quantile(b2_sample,0.975),col='blue')


hist(b3_sample,main="b3")
abline(v=0,col='red')
abline(v=quantile(b3_sample,0.025),col='blue')
abline(v=quantile(b3_sample,0.975),col='blue')

hist(b4_sample,main="b4")
abline(v=quantile(b4_sample,0.025),col='blue')
abline(v=quantile(b4_sample,0.975),col='blue')

abline(v=0,col='red')
} 


par(mfrow=c(2,2))
# WT mice
for (period in c(1,2,3,4)) {
  main = paste("WT, t=", period ,' vs period 0')
  ORs <- getORsTime(samples,type="WT",period=period,"mouseModel") 
  hist(ORs,main=main )
  abline(v=quantile(ORs,0.025),col='blue')
  abline(v=quantile(ORs,0.975),col='blue')
  abline(v=1,col='red')
}

# pound mice
for (period in c(1,2,3,4)) {
  main = paste("Pound, t=", period ,' vs period 0')
  ORs <- getORsTime(type="T2",period=period,model="mouseModel") 
  hist(ORs,main=main )
  abline(v=quantile(ORs,0.025),col='blue')
  abline(v=quantile(ORs,0.975),col='blue')
  abline(v=1,col='red')
}
  



  




