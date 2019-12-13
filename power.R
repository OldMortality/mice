library(MASS)


head(df)
powers <- vector()
# t.test on the difference
N.SIMS <- 10000
for (N in 35:60) {
  #
  base <- 75
  sd <- 15  #18.6
  effect.size <-  -7 #10.93
  ids <- c(seq(1:N),seq(1:N))
  errs <- 0
  for (sim in 1:N.SIMS) {
    if (sim %% 1000 == 0) { print(c(N,sim,1-errs/sim))}  
    # control group
    Means <- c(base,base)
    d.cov <- 0.7 * sd ^ 2 
    Sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 
    simulation <- mvrnorm(n = N, Means, Sigma)
    t0 <- simulation[,1]
    t1 <- simulation[,2]
    df.control <- data.frame(response=c(t0,t1),
                             id = seq(1,N),
                             time = c(rep('0',N),
                                      rep('1',N)),
                             type='control')
    # treatment group
    Means <- c(base,base+effect.size)
    d.cov <- 0.7 * sd ^ 2 
    Sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 
    simulation <- mvrnorm(n = N, Means, Sigma)
    t0 <- simulation[,1]
    t1 <- simulation[,2]
    df.treat <- data.frame(response=c(t0,t1),
                           id = seq(N+1,2*N),
                           time = c(rep('0',N),
                                    rep('1',N)),
                           type='treat')
    
    df <- rbind(df.control,df.treat)
    df$id <- factor(df$id)
    df$time <- factor(df$time)
    #sqrt(var(x))
    #sqrt(var(y))
    #cor(x,y)
    
    diff <- vector()
    for (id in seq(1,2*N)) {
      diff[id] <-  df[which(df$id==id),c("response","type")][2,1]-df[which(df$id==id),c("response","type")][1,1]
    }
    #df[which(df$id==33),]
    
    type <- c(rep('control',N),rep('treat',N))
    df2 <- data.frame(diff=diff,type=type)
    pval <- t.test(diff~type,data=df2)$p.val
    
    
    
    
    if ( pval > 0.05 ) {
      # ci contains zero. No difference detected
      errs <- errs + 1
    }
    
    
  }
  print(c(N,1-errs/N.SIMS))
  powers[N] <- 1 - errs/N.SIMS
}


sd <- 1.4
d.cov <- 0.8 * sd ^ 2 
sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 


diff <- -0.57
x <- mvrnorm(n = 10000, mu= c(8.9,8.9 + diff) ,Sigma = sigma)
x1 <- x[,1]
x2 <- x[,2]
mean(x1)
mean(x2)
cor(x1,x2)
sqrt(var(x1))
sqrt(var(x2))

plot(density(x1))
lines(density(x2),col='red')
abline(v=7,col='blue')
plot(density(x1-x2))
mean(x1-x2)
var(x1-x2)



rm(list = ls())
sd <- 1.4
d.cov <- 0.8 * sd ^ 2 
sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 

diff <- -0.57
pow <- 0
N.SIMS <- 10000
x.m <- vector()
s.m <- vector()
for (sim in 1:N.SIMS) {
  x <- mvrnorm(n = 41, mu= c(8.9,8.9 + diff) ,Sigma = sigma)
  x1 <- x[,1]
  x2 <- x[,2]
  x.m[sim] <- mean(x1-x2)
  s.m[sim] <- sqrt(var(x1-x2))
  if (t.test(x1-x2)$p.val < 0.05) {
    pow <- pow + 1
  }
}
pow/N.SIMS

hist(s.m,prob=T)
abline(v=0.8,col='red')
mean(s.m)
hist(x.m)
mean(x.m)
sqrt(var(x.m))

N <- 41
# sampling distribution of x1-x2
par(mfrow=c(1,1))
x.bar.h0 <- rnorm(10000,mean=0   ,sd=sqrt(0.8/N))
x.bar.h1 <- rnorm(10000,mean=0.57,sd=sqrt(0.8/N))
plot(density(x.bar.h0),xlim=c(-1,1))
lines(density(x.bar.h1),col='red',xlim=c(-1,1))



#acceptance region for H0
upp <- 1.96 * sqrt(0.8/N)
abline(v= upp  ,col='blue')
# power
sum(x.bar.h1 > upp) /length(x.bar.h1)
