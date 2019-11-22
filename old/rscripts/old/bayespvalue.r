library(R2jags)

# size in the binomial draw
N <- 100
# number of observations
M <- 100
# true value of p
p_true = 0.3

#### the observations
#### comment out one of the following:

# this is how we model it. The suspect should
#    not stand out in the line-up.
y <- rbinom(M, size=N,p=p_true)
# this is something else. You should be able to
#    spot the suspect.
y <- round(runif(M,min=0,max=N))



mymodel <- function()
{
  for ( i in 1 : M) {
    y[i] ~ dbin(p[i],N)
    logit(p[i]) <- alpha0 
  }
  
  # priors
  alpha0 ~ dunif(-1,1)
   
}


my_inits <- list(
  list(alpha0 = 0.2 ),
  list(alpha0 = 0.25),
  list(alpha0 = 0.65),
  list(alpha0 = 0.8))


parameters <- c("alpha0","p")
samples <- jags(model = mymodel ,
                data = list('y' = y,
                            'N' = N,
                            'M' = M),
                parameters,
                n.chains=4,
                n.iter=10000,
                n.burnin=5000,
                inits=my_inits,
                n.thin=1,
                DIC=T)




 
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

sims <- samples$BUGSoutput$sims.list

a <- samples$BUGSoutput$sims.list$alpha0
p <- samples$BUGSoutput$sims.list$p

#hist(expit(a),30)

#traceplot(samples$BUGSoutput,ask=T)

# compute D
ds <- numeric()
counter <- 0
for (t in 1000:4000) {
  result <- 0
  counter <- counter + 1
  for (i in 1:M) {
     result <- result + ((y[i]-N*p[t,i])^2)/(N*p[t,i])
  }
  ds[counter] <- result
}
head(ds)  
hist(ds,30)  

# compute D*
# data simulated from the model
ds2 <- numeric()
counter <- 0
for (t in 1000:4000) {
  result <- 0
  counter <- counter + 1
  for (i in 1:M) {
    y2 <- rbinom(1,size=N,p=p[t,i])
    result <- result + ((y2-N*p[t,i])^2)/(N*p[t,i])
  }
  ds2[counter] <- result
}
#head(ds2)  
#hist(ds2,30)  


plot(ds,ds2,
     xlim=c(min(ds,ds2),max(ds,ds2)),
     ylim=c(min(ds,ds2),max(ds,ds2)))
     
abline(0,1,col='red')
sum(ds2<ds)
p.value <- sum(ds2<ds)/length(ds)
print("==============")
p.value


# police line-up
max_lineup = 16 #(must match the next line)
par(mfrow=c(4,4))
y2 <- matrix(nrow=max_lineup,ncol=M)
# leave one position for the observations
for (lineup in 1:max_lineup - 1) {
  t <- round(runif(1,min=1000,max=4000))
  for (i in 1:M) {
    y2[lineup,i] <- rbinom(1,size=N,p[t,i])  
  }
}


suspect <- round(runif(1,min=1,max=max_lineup))
counter <- 1
for (pos in 1:max_lineup) {
  if (pos == suspect) {
    # show observation
    hist(y[],30,main=toString(pos),xlab="y")    
  } else {
    
    # show data simulated from the model
    hist(y2[counter,],30.,main=toString(pos),xlab="y")  
    counter <- counter + 1
  }
}

 print(suspect)



