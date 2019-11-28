#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#
#library(rjags)
library(R2jags)

macFile <- "/Users/micheldelange/Documents/mice/data/mice.csv"
ubuntuFile <- "~/mice/data/mice.csv"

if (file.exists(macFile)) {
  mice <- read.csv(macFile,header=TRUE)[1:40,]
} else {
  mice <- read.csv(ubuntuFile,header=TRUE)[1:40,]
}


mice$gender <- grepl('female',mice$Animal_type)
mice$type_1 <- grepl('PLT',mice$Animal_type)
mice$type_2 <- grepl('Pound',
                mice$Animal_type)
mice$type <- 'WT'
mice[which(mice$type_1),"type"] <- 'PLT'
mice[which(mice$type_2),"type"] <- 'T2'

celltype <- "CD19pos_B220"
# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

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
                

  N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
  R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))

  
  # long notation (i.e. observations as rows)
  per_long <- vector()
  mouse_long <- vector()
  gender_long <- vector()
  type_1_long <- vector()
  type_2_long <- vector()
 
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
    # periods run from [0,4] ~ [0,3,10,17,24]
    per_long    <- c(per_long,j-1)
    mouse_long  <- c(mouse_long,i)
    female_long <- c(female_long,female[i])
    type_1_long <- c(type_1_long,type_1[i])
    type_2_long <- c(type_2_long, type_2[i])
    N_long      <- c(N_long,N_matrix[i,j])
    R_long      <- c(R_long,R_matrix[i,j])
    
    }
  }
  t3 <- per_long == 3
  t10 <- per_long == 10
  t17 <- per_long = 17
  t24 <- per_long = 24
  
  theList <- list('R' = R_long,
                  'Ncells' = N_long,
                  'period' = per_long,
                  'mouse'  = mouse_long,
                  'female' = female_long,
                  'type_1' = type_1_long,
                  'type_2' = type_2_long,
                  'time3' =  t3,
                  'time10' = t10,
                  'time17' = t17,
                  'time24' = t24,
                  'N' = length(R_long))
  return(theList)
}

theList <- createList(mice)

mouseModel <- function()
{
  for ( i in 1 : N) {
    R[i] ~ dbin(p[mouse[i], period[i]],Ncells[i])
    b[mouse[i], period[i]] ~ dnorm(0.0, tau); 
    logit(p[mouse[i], period[i]]) <- a0 + a1 * gender[i] +
      a2 * type_1[i] +
      a3 * type_2[i] +
      a4 * time3[i] +
      a5 * time10[i] +
      a6 * time17[i] +
      a7 * time24[i] +
      a8 * time31[i] +
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
  a8 ~ dnorm(0.0,1.0E-6) 
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




# initialise the chains


setInits <- function(i1,i2,i3,i4) { 
return( list(
  list(a0 = i1, 
       a1 = i1,
       a2 = i1 ,
       a3 = i1,
       a4 = i1,
       a5 = i1,
       a6 = i1,
       a7 = i1,
       a8 = i1,
       b1 = i1,
       b2 = i1,
       b3 = i1,
       b4 = i1,
       g1 = i1,
       g2 = i1,
       g3 = i1,
       g4 = i1,
       sigma = 0.9),
  list(a0 = i2, 
       a1 = i2,
       a2 = i2 ,
       a3 = i2,
       a4 = i2,
       a5 = i2,
       a6 = i2,
       a7 = i2,
       a8 = i2,
       b1 = i2,
       b2 = i2,
       b3 = i2,
       b4 = i2,
       g1 = i2,
       g2 = i2,
       g3 = i2,
       g4 = i2,
       sigma = 1.1),
  list(a0 = i3, 
       a1 = i3,
       a2 = i3,
       a3 = i3,
       a4 = i3,
       a5 = i3,
       a6 = i3,
       a7 = i3,
       a8 = i3,
       b1 = i3,
       b2 = i3,
       b3 = i3,
       b4 = i3,
       g1 = i3,
       g2 = i3,
       g3 = i3,
       g4 = i3,
       sigma = 1),
  list(a0 = i4, 
       a1 = i4,
       a2 = i4,
       a3 = i4,
       a4 = i4,
       a5 = i4,
       a6 = i4,
       a7 = i4,
       a8 = i4,
       b1 = i4,
       b2 = i4,
       b3 = i4,
       b4 = i4,
       g1 = i4,
       g2 = i4,
       g3 = i4,
       g4 = i4,
       sigma = 0.7)
  
)
)}

##
## Rename bugs output, so it is a bit more user friendly
## 
renameSamples <- function(samples){ 
    result <- data.frame(
      a0 = samples$BUGSoutput$sims.list$a0,
      a1 = samples$BUGSoutput$sims.list$a1,
      a2 =  samples$BUGSoutput$sims.list$a2,
      a3 =  samples$BUGSoutput$sims.list$a3,
      a4 =  samples$BUGSoutput$sims.list$a4,
      a5 =  samples$BUGSoutput$sims.list$a5,
      a6 =  samples$BUGSoutput$sims.list$a6,
      a7 =  samples$BUGSoutput$sims.list$a7,
      a8 =  samples$BUGSoutput$sims.list$a8,
      b1 =  samples$BUGSoutput$sims.list$b1,
      b2 =  samples$BUGSoutput$sims.list$b2,
      b3 =  samples$BUGSoutput$sims.list$b3,
      b4 =  samples$BUGSoutput$sims.list$b4,
      sigma =  samples$BUGSoutput$sims.list$sigma)
      #p =  samples$BUGSoutput$sims.list$p,
      #b = samples$BUGSoutput$sims.list$b)
  return(result)
}

#res <- renameSamples(samples)


parameters <- c("a0", "a1","a2","a3","a4","a5",
                "a6","a7","a8","b1","b2","b3","b4",
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
  return(which(!mice$gender))
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
  



  




