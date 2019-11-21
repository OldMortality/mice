#
#
# 2019. Took this from mice2.rmd.
#
#
library(R2jags)
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
type_2 <- grepl('Pound',
                mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))

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
      beta1 * time3[i] * type_2[i] +
      beta2 * time10[i] * type_2[i] +
      beta3 * time17[i] * type_2[i] +
      beta4 * time24[i] * type_2[i] +
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
  beta1 ~ dnorm(0.0,1.0E-6) 
  beta2 ~ dnorm(0.0,1.0E-6) 
  beta3 ~ dnorm(0.0,1.0E-6) 
  beta4 ~ dnorm(0.0,1.0E-6) 
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}



time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)

# initialise the chains


setInits <- function(i1,i2,i3,i4) { 
return( list(
  list(alpha0 = i1, 
       alpha1 = i1,
       alpha2 = i1 ,
       alpha3 = i1,
       alpha4 = i1,
       alpha5 = i1,
       alpha6 = i1,
       alpha7 = i1,
       alpha8 = i1,
       beta1 = i1,
       beta2 = i1,
       beta3 = i1,
       beta4 = i1,
       sigma = 0.9),
  list(alpha0 = i2, 
       alpha1 = i2,
       alpha2 = i2 ,
       alpha3 = i2,
       alpha4 = i2,
       alpha5 = i2,
       alpha6 = i2,
       alpha7 = i2,
       alpha8 = i2,
       beta1 = i2,
       beta2 = i2,
       beta3 = i2,
       beta4 = i2,
       sigma = 1.1),
  list(alpha0 = i3, 
       alpha1 = i3,
       alpha2 = i3,
       alpha3 = i3,
       alpha4 = i3,
       alpha5 = i3,
       alpha6 = i3,
       alpha7 = i3,
       alpha8 = i3,
       beta1 = i3,
       beta2 = i3,
       beta3 = i3,
       beta4 = i3,
       sigma = 1),
  list(alpha0 = i4, 
       alpha1 = i4,
       alpha2 = i4,
       alpha3 = i4,
       alpha4 = i4,
       alpha5 = i4,
       alpha6 = i4,
       alpha7 = i4,
       alpha8 = i4,
       beta1 = i4,
       beta2 = i4,
       beta3 = i4,
       beta4 = i4,
       sigma = 0.7)
  
)
)}



parameters <- c("alpha0", "alpha1","alpha2","alpha3","alpha4","alpha5",
                "alpha6","alpha7","alpha8","beta1","beta2","beta3","beta4",
                "sigma","p","b")
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
                n.chains = 4,
                n.iter = 5000,
                n.burnin = 1000,
                inits = setInits(i1 = 0,i2 = 0.1,i3 = -0.1,i4 = 0.245),
                n.thin = 1,
                DIC = T)



#samples$BUGS$summ
# not sure how this samples from the 4 MC's. I think I am 
#   probably just sampling from the first one this way.
alpha0_sample <- samples$BUGSoutput$sims.list$alpha0
alpha1_sample <- samples$BUGSoutput$sims.list$alpha1
alpha2_sample <- samples$BUGSoutput$sims.list$alpha2
alpha3_sample <- samples$BUGSoutput$sims.list$alpha3
alpha4_sample <- samples$BUGSoutput$sims.list$alpha4
alpha5_sample <- samples$BUGSoutput$sims.list$alpha4
alpha6_sample <- samples$BUGSoutput$sims.list$alpha4
alpha7_sample <- samples$BUGSoutput$sims.list$alpha4
alpha8_sample <- samples$BUGSoutput$sims.list$alpha4
beta1_sample <- samples$BUGSoutput$sims.list$beta1
beta2_sample <- samples$BUGSoutput$sims.list$beta2
beta3_sample <- samples$BUGSoutput$sims.list$beta3
beta4_sample <- samples$BUGSoutput$sims.list$beta4
sigma_sample <- samples$BUGSoutput$sims.list$sigma
p_sample <- samples$BUGSoutput$sims.list$p
b_sample <-samples$BUGSoutput$sims.list$b

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}


###
# Compare type_2 (Pound) with wild type mice (type_1=type_2 == 0) for
# male mice at time 1 (all time_1 ==0)
type2_t0 <- which(type_2_long & !gender_long & !
        (time_3_long | time_10_long | time_17_long | time_24_long | time_31_long))

R_long[type2_t0]
N_long[type2_t0]


ytype2 <- R_long[type2_t0]/N_long[type2_t0]
xtype2 <- rep(0,length(R_long[type2_t0]))
plot(ytype2~xtype2,ylim=c(0,.4),xlim=c(0,5),
     col='red')
yWT <- R_long[WT_t0]/N_long[WT_t0]
xWT <- rep(0.2,length(R_long[WT_t0]))
points(yWT~xWT,col='black')

# get all male WT mice

getMaleWT <- function() {
  return(which(!gender_long ))
}

# returns all records for a given period
#   period 0  0 days
#          1  3
#          2  10
#          3  17
#          4  24
#          5  31 days
#
getByPeriod <- function(period) {
  result <- vector()
  if (period == 0) 
    result <- which(!(time_3_long | time_10_long | time_17_long | time_24_long | time_31_long))
  if (period == 1) 
    result <- which(time_3_long)
  if (period == 2) 
    result <- which(time_10_long)
  if (period == 3) 
    result <- which(time_17_long)
  if (period == 4) 
    result <- which(time_24_long)
  if (period == 5) 
    result <- which(time_31_long)
  
  return(result)
}
getByPeriod(1)

getMaleWTbyPeriod(period) {
  WT_t0 <- which(! (type_2_long | type_1_long) & !gender_long & !
                   (time_3_long | time_10_long | time_17_long | time_24_long | time_31_long))
  R_long[WT_t0]
  N_long[WT_t0]
  
}

getMaleWTByPeriod <- function(period) {
  which(getMaleWT() ) %in% getByPeriod(period))
}

mice[getMaleWT(),]$Animal_type

mice[1,]$Animal_type

colnames(mice)

par(mfrow=c(1,2))

getORs <- function(period) {
  Odds_WT <- exp(alpha0_sample)
  etaT2 <- alpha0_sample +  alpha3_sample
  if (period == 1) {
    etaT2 <- etaT2 + beta1_sample
  }
  if (period == 2) {
    etaT2 <- etaT2 + beta2_sample
  }
  if (period == 3) {
    etaT2 <- etaT2 + beta3_sample
  }
  if (period == 4) {
    etaT2 <- etaT2 + beta4_sample
  }
  Odds_type2 <- exp(etaT2)
  OR <- Odds_type2/Odds_WT
  return(OR)
}

par(mfrow=c(1,1))
OR <- getORs(period=5)
hist(OR)

abline(v=quantile(OR,0.025),col='blue')
abline(v=quantile(OR,0.975),col='blue')
abline(v=1,col='red')


