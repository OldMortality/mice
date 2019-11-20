N <- 100000
x1 <- rnorm(N,10,1)
x2 <- rnorm(N,100,5)
x <- x1 + x2
hist(x1,30)
hist(x2,30)

hist(x,30)

x_std <- (x - mean(x))/(var(x))^0.5
hist(x_std,30)

