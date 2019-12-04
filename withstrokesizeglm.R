mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]
mice.s <- read.csv("~/Documents/mice/data/mice2.csv",header=TRUE)[1:40,]

mice$stroke_size <- scale(mice.s$stroke_size)

x <- rnorm(1000,10,5)
var(scale(x))


m <- lmer(y~ stroke_size + type * factor(period) +  (1|mouse) , data=mice2)

s <- summary(m)
CI_lower <- s$coefficients[,2] - 1.96*s$coefficients[,2]
CI_upper <- s$coefficients[,2] + 1.96*s$coefficients[,2]
cbind(CI_lower,CI_upper)
