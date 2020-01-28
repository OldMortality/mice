N <- 50
# number of experiments
M <- 100

pvals <- vector()
s.fwer <- vector()
s.raw <- vector()

for (sim in 1:1000) {
  yref <- rnorm(N,mean=0)

  for (i in 1:M) {
    if (i <= 20) {
      y <- rnorm(N,mean=1)
    } else {
      y <- rnorm(N)
    }
    t <- t.test(y,yref)
    pvals[i] <- t$p.value 
  }

  df.p <- data.frame(pvals=pvals,number=1:length(pvals))
  df.p$pvals
  colnames(df.p)
  f <- FWER(df.p,plot=F)
  s.fwer[sim] <- sum(f$reject)
  s.raw[sim]  <-sum(df.p$pvals < 0.05)

}

plot(s.fwer~s.raw)
abline(a=0,b=1,col='red')
mean(s.fwer)
mean(s.raw)
hist(s.fwer)
hist(s.raw)



errs <- vector()
N <- 5
for (sim in 1:10000) {
  y1 <- rnorm(N,mean=0)
  y2 <- rnorm(N,mean=0)
  y3 <- rnorm(N,mean=0)
  y4 <- rnorm(N,mean=0)
  y5 <- rnorm(N,mean=0)
  
  p2 <- t.test(y2,y1)$p.value
  p3 <- t.test(y3,y2)$p.value
  p4 <- t.test(y4,y3)$p.value
  p5 <- t.test(y5,y4)$p.value

  err2 <- p2 < 0.05
  err3 <- p3 < 0.05
  err4 <- p4 < 0.05
  err5 <- p5 < 0.05
  errs[sim] <- err2 + err3 + err4 + err5
}
plot(errs,pch='.')
sum(errs==0)  
sum(errs==1)
sum(errs==2)
sum(errs==3)
sum(errs==4)  





  

plot(s.fwer~s.raw)
abline(a=0,b=1,col='red')
mean(s.fwer)
mean(s.raw)
hist(s.fwer)
hist(s.raw)
