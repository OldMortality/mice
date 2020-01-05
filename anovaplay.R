
N <- 30
x <- rnorm(N)
y <- rnorm(N)


obs <- 2 * x + 3 * y + rnorm(N)
plot(obs ~ y)
m0 <- lm(obs ~ y * x)
m1 <- lm(obs ~ y + x)
summary(m1)
errs <- 0
for (i in 1:1000) {
  z <- rnorm(N)
  m2 <- lm(obs ~ y + z)  
  a <- anova(m2)
  if (a$`Pr(>F)`[2] < 0.05) {
    errs <- errs + 1
  }
}
errs
anova(m0,m1)
