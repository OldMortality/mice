rm(list = ls())
library(MASS)



getpow2 <- function(n, diff, sd, rho,  N.SIMS = 1000) {
  
  d.cov <- rho * sd ^ 2 
  sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 
  pow <- 0
  base <- 8.6
  for (sim in 1:N.SIMS) {
    xt <- mvrnorm(n , mu= c(base,base + diff) ,Sigma = sigma)
    # control group
    xc <- mvrnorm(n , mu= c(base,base ) ,Sigma = sigma)
    if (t.test(xt[,1]-xt[,2],
               xc[,1]-xc[,2])$p.val < 0.05) {
      pow <- pow + 1
    }
  }
  pow/N.SIMS
  
}
par(mfrow=c(2,1))
plot(density(difft),xlim=c(-1,1),col='red')
lines(density(diffc),col='black')
grid()

sqrt(var(diffc))
1.96 * 0.885/sqrt(41)

sum(x.bar.h1<1.96 * 0.885/sqrt(41))/length(x.bar.h1)

library(pipeR)

par(mfrow=c(1,1))

lowest <- 30
highest <- 50
ns <- seq(lowest,highest,1)
ns %>% map(getpow2,diff = 0.57,sd=1.4, rho = 0.8, N.SIMS=5000) %>% 
  unlist %>% 
  ( function(x) { c(rep(NA,lowest),x) } ) %>>% 
  (~ powers) %>%
  plot
lines(powers)  
grid()
abline(h=0.8,col='blue')
min(which(powers>0.8))

#pows <- map(ns,getpow)

#ns %>>%  (~ fred)



#unlist(pows2)


N <- 41
# sampling distribution of x1-x2
par(mfrow=c(2,1))
x.bar.h0 <- rnorm(10000,mean=0   ,sd=sqrt(0.784/N))
x.bar.h1 <- rnorm(10000,mean=0.57,sd=sqrt(0.784/N))
plot(density(x.bar.h0),xlim=c(-1,1))
lines(density(x.bar.h1),col='red',xlim=c(-1,1))


#acceptance region for H0
upp <- 1.96 * sqrt(0.784) /sqrt(N)
low <- -1.96 * sqrt(0.784/N)
abline(v= c(low,upp)  ,col='blue')
# area under red curve, to theleft of upp
sum(x.bar.h1 < upp) /length(x.bar.h1)
# power
1-sum(x.bar.h1 < upp) /length(x.bar.h1)

# power
1-pnorm(upp,mean=0.57,sd=sqrt(0.8/N))


#effect of rho

rho = 0.2
d.cov <- rho * sd ^ 2 
sigma <- matrix(c(sd^2,d.cov,d.cov,sd^2),nrow = 2) 
base <- 8.6
x <- mvrnorm(n , mu= c(base,base + diff) ,Sigma = sigma)
xc <- mvrnorm(n , mu= c(base,base ) ,Sigma = sigma)
par(mfrow=c(2,1))
plot(density(x2-x1),xlim=c(-10,10),ylim=c(0,1),main=rho)
lines(density(x2c-x1c),col='red')


s <- vector()
for (i in 1:1000) {
  s[i] <- mean(rnorm(30,mean=0,sd=0.784))
}
sqrt(var(s))
0.784/sqrt(30)
