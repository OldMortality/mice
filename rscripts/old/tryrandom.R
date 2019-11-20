
set.seed(123)
type <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
mouse <- factor(c(1,2,3,4,5,6,7,8,9,10,11,12))

# pab = probability for type a, at time b
p11 <- 0.4
p12 <- 0.55
p21 <- 0.5
p22 <- 0.67

# random effects
r <- numeric()
for (i in 1:12) {
  r[i] <- runif(1,min=-0.2,max=0.2)
}

# time 1 
y1 <- numeric()
for (i in 1:12) {
  if (i<7) {
    # type 1
    y1[i] <- rbinom(1,size=100,p=p11 + r[i]  )/100 
  } else {
    # type 2
    y1[i] <- rbinom(1,size=100,p=p21 + r[i])/100
  }
}

time <- factor(c('1','1','1','1','1','1','1','1','1','1','1','1'))
d1 <- data.frame(mouse,type,time,y1)

# time 2 
y2 <- numeric()
for (i in 1:12) {
  if (i<7) {
    # type 1
    y2[i] <- rbinom(1,size=100,p=p12 + r[i]  )  / 100
  } else {
    # type 2
    y2[i] <- rbinom(1,size=100,p=p22 + r[i])  / 100
  }
}

time2 <- c('2','2','2','2','2','2','2','2','2','2','2','2')
d2 <- data.frame(mouse,type,time2,y2)
colnames(d2) <- c('mouse','type','time','y1')

d <- data.frame(rbind(d1,d2))
#d$time
#d$y1
#d$y1 <- as.numeric(d$y1)
## just fixed effects

boxplot(d$y1[1:12]~d$type[1:12],main='time 1')
boxplot(d$y1[13:24]~d$type[13:24],main='time 2')

ml1 <- glm(d$y1~1 + d$type,family="quasibinomial")
summary(ml1)


ml2 <- glm(d$y1~1 + d$type + d$time,family="quasibinomial")
summary(ml2)

ml3 <- glm(d$y1~1 + d$type * d$time,family="quasibinomial")
summary(ml3)

## random effects model
library(lme4)

# use numbers, not proportions
#d$y1 = d$y1 / 100 
d$successes <- d$y1 * 100
d$failures <- 100- d$successes
response <- cbind(d$successes,d$failures)

ml3 <- glmer(response ~ type + time + type:time +  (1|mouse),family="binomial",data=d)
summary(ml3)

ml4 <- glmer(response ~ d$type + d$time +  (1|mouse),family="binomial",data=d)
summary(ml4)


ml5 <- glmer(response ~ d$type  +  (1|mouse),family="binomial",data=d)
summary(ml5)

ml6 <- glmer(response ~ 1  +  (1 + time |mouse),family="binomial",data=d)
summary(ml6)

