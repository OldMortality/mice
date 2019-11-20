# uses R2jags, as opposed to rjags.
## uses wide data. Does not work because of NA's. I had replaced them with zero's

library(R2jags)
mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]
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
type_2 <- grepl('Pound',mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))

 






mouseModel <- function()
{
  for ( j in 1 : N) {
    for ( k in 1 : 6) {
      R_matrix[j, k] ~ dbin(p[j, k],N_matrix[j,k])
      b[j, k] ~ dnorm(0.0, tau); 
      logit(p[j, k]) <- alpha0 + alpha1 * gender[j] +
      alpha2 * type_1[j] +
      alpha3 * type_2[j] +
      alpha4 * time3[k] +
      alpha5 * time10[k] +
      alpha6 * time17[k] +
      alpha7 * time24[k] +
      alpha8 * time31[k] +
      b[j,k]
    }
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
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}
 


time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)

i1 <- 0
i2 <- 0.1
i3 <- -0.1
i4 <- 0.245

my_inits <- list(
  list(alpha0 = i1, 
       alpha1 = i1,
       alpha2 = i1 ,
       alpha3 = i1,
       alpha4 = i1,
       alpha5 = i1,
       alpha6 = i1,
       alpha7 = i1,
       alpha8 = i1,
       sigma =0.9),
  list(alpha0 = i2, 
       alpha1 = i2,
       alpha2 = i2 ,
       alpha3 = i2,
       alpha4 = i2,
       alpha5 = i2,
       alpha6 = i2,
       alpha7 = i2,
       alpha8 = i2,
       sigma =1.1),
  list(alpha0 = i3, 
       alpha1 = i3,
       alpha2 = i3,
       alpha3 = i3,
       alpha4 = i3,
       alpha5 = i3,
       alpha6 = i3,
       alpha7 = i3,
       alpha8 = i3,
       sigma =1),
  list(alpha0 = i4, 
       alpha1 = i4,
       alpha2 = i4,
       alpha3 = i4,
       alpha4 = i4,
       alpha5 = i4,
       alpha6 = i4,
       alpha7 = i4,
       alpha8 = i4,
       sigma =0.7)
  
)


parameters <- c("alpha1","alpha2","alpha3","alpha4")
samples <- jags(model = mouseModel ,
                   data = list('R_matrix' = R_matrix,
                               'N_matrix' = N_matrix,
                               'gender' = gender,
                               'type_1' = type_1,
                               'type_2' = type_2,
                               'time3' =  time_3,
                               'time10' = time_10,
                               'time17' = time_17,
                               'time24' = time_24,
                               'time31' = time_31,
                               'N' = 40 ),
                   parameters,
                   n.chains=4,
                   n.iter=10000,
                   n.burnin=1000,
                   inits=my_inits,
                   n.thin=1,
                   DIC=T)


 
samples$BUGS$summ
alpha1_sample <- samples$BUGSoutput$sims.list$alpha1
alpha2_sample <- samples$BUGSoutput$sims.list$alpha2
alpha3_sample <- samples$BUGSoutput$sims.list$alpha3
alpha4_sample <- samples$BUGSoutput$sims.list$alpha4

hist(alpha1_sample,30,probability=T)
hist(alpha2_sample,30,probability=T)
hist(alpha3_sample,30,probability=T)
hist(alpha4_sample,30,probability=T)

#plot(alpha1_sample,pch='.')
#traceplot(samples$BUGSoutput,ask=F)

