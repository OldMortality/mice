
p1 <- 0.8

# probability of g2 | g1


jdistance <- function(p2) {

  J <- vector()
  for (sim in 1:10000) {
    g1 <- rbinom(5,size=1,prob=p1)
    g2 <- g1 * rbinom(5,size=1,prob=p2)
    # Jaccard distance
    J[sim] <- 1- (sum(g1 & g2))/sum(g1 | g2)
  }
  J <- round(J,1)
  table(J)
  return(mean(J,na.rm=T))
}

JS <- vector()
counter <- 0

for (p2 in seq(0.01,1,length.out = 100)) {
  counter <- counter+ 1
  JS[counter] <- jdistance(p2)
}
par(mfrow=c(1,1))
plot(JS~seq(0.01,1,length.out = 100),
     ylab='mean Jdistance',
     xlab='P(g2|g1)')
