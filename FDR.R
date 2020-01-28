
getDR <- function(M = 100) {
N <- 20

# 50 true different, 950 zeros
d <- c(rnorm(10*N,mean=1),
       rnorm((M-10)* N,mean=0))
e <- rnorm(M* N,mean=0)

y1 <- matrix(data=d,nrow=M,byrow=T)
y2 <- matrix(data=e,nrow=M,byrow=T)
ps <- vector()
for (i in 1:M) {
  t <- t.test(y1[i,],y2[i,])
  ps[i] <- t$p.value
}

p2 <- data.frame(ps = ps)
p2$ranks <- NA
head(p2)
p2 <- p2[with(p2,order(ps)),]
p2$ranks <- seq(1,M)

head(p2)
 

#plot(p2$ps)
FDR = 0.05
p2$th <- FDR * p2$ranks/ M
r1 <- max(which(p2$ps < p2$th))
r2 <- length(which(p2$ps < 0.05))
return(r1)
}


r <- vector()
for (i in 1:1000) {
  r[i] <- getDR()
}
hist(r)
hist((r-10)/r,50)
mean((r-10)/r)

mean(r)

p <- getDR(M=100000)
hist(p$p)
