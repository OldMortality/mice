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
type_2 <- grepl('Pound',mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))

for (i in 1:dim(N_matrix)[1]) {
  for (j in 1:dim(N_matrix)[2]) {
    if (is.na(N_matrix[i,j])) {
      N_matrix[i,j] <- 1
    }
  } 
}
for (i in 1:dim(R_matrix)[1]) {
  for (j in 1:dim(R_matrix)[2]) {
    if (is.na(R_matrix[i,j])) {
      R_matrix[i,j] <- 1
    }
  } 
}






model1.string<-"
model {
  for ( j in 1 : N) {
    for ( k in 1 : 6) {
      
      R_matrix[j, k] ~ dbin(mu[j, k],N_matrix[j,k])
      b[j, k] ~ dnorm(0.0, tau); 
      logit(mu[j, k]) <- alpha0 + alpha1 * gender[j] +
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
"

model1.spec <- textConnection(model1.string)

time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)

my_inits <-
  list(alpha0 = 0, 
       alpha1 = 0,
       alpha2 = 0 ,
       alpha3 = 0,
       alpha4 = 0,
       alpha5 = 0,
       alpha6 = 0,
       alpha7 = 0,
       alpha8 = 0,
       sigma =1)



jags <- jags.model(model1.spec,
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
                   n.chains=4,
                   n.adapt=100,
                   inits=my_inits)


#traceplot(jags$)

#jags$


update(jags, 10000)

#jags.samples(jags,
#             c('alpha1'),
#             1000)
  
alpha1_post <- coda.samples(jags,c('alpha1'),1000)
alpha2_post <- coda.samples(jags,c('alpha2'),1000)

#alpha1_post[[1]]
 
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
} 

traceplot(alpha1_post)
traceplot(alpha2_post)

hist(expit(alpha1_post[[2]]),30,main='gender')
hist(expit(alpha1_post[[2]]),30,main='gender')



  