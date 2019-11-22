# install.packages("rjags")
library(rjags)
set.seed(432104)
n <- 1000
x <- rnorm(n, 0, 5)

model1.string <-"
  model {
    for (i in 1:N){
    x[i] ~ dnorm(mu, tau)
    }
  mu ~ dnorm(0,.0001)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}
"
model1.spec<-textConnection(model1.string)

jags <- jags.model(model1.spec,
                   data = list('x' = x,
                               'N' = n),
                   n.chains=4,
                   n.adapt=100)

update(jags, 1000)

jags.samples(jags,
             c('mu', 'tau'),
             1000)


