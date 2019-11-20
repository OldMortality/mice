##
##  Away with all that Bayesian stuff
##  Long data
##  quasibinomial model
##
##
library(lme4)
library(knitr)
library(stringr)

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice2.csv",header=TRUE)[1:35,]

# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24)
periods <- 5

N_0 <- mice$Live_cells_0
R_0 <- mice$CD19pos_B220pos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD19pos_B220pos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD19pos_B220pos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD19pos_B220pos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD19pos_B220pos_24
N_31 <- mice$Live_cells_31
R_31 <- mice$CD19pos_B220pos_31

female <- grepl('female',mice$Animal_type)
type_1 <- grepl('PLT',mice$Animal_type)
type_2 <- grepl('Pound',mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24))
P_matrix <- R_matrix/N_matrix


time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)



 
# long notation (i.e. observations as rows)
per_long <- vector()
mouse_long <- vector()
female_long <- vector()
type_long <- vector()
N_long <- vector()
R_long <- vector()
stroke_size_long <- vector()
 

counter = 0
skipped = 0
for (i in 1: N) {
  for (j in 1:5) {
    
    counter <- counter + 1
    if (is.na(N_matrix[i,j])) {
      skipped <- skipped + 1
      next
    }
    per_long    <- c(per_long,str_c("PERIOD",j))
    mouse_long  <- c(mouse_long,i)
    female_long <- c(female_long,female[i])
    stroke_size_long <- c(stroke_size_long,as.numeric(toString(mice$Stroke.volume..mm3.[i])))
    
    
    type = "aWT"
    if (type_1[i] == 1) {
      type <- "PLT"
    }
    if (type_2[i] == 1) {
      type <- "Pound"
    }
    type_long <- c(type_long,type)
    N_long      <- c(N_long,N_matrix[i,j])
    R_long      <- c(R_long,R_matrix[i,j])
    
  }
  
}

long_data <- data.frame(mouse_long,
                   female_long,
                   type_long,
                   per_long,
                   N_long,
                   R_long,
                   stroke_size_long)


dim(long_data)
#write.csv(long_data,"/Users/micheldelange/Documents/mice/data/long_data.csv")

d <- long_data

colnames(d) <- c("mouse","female","type","period","total","cells","stroke_size")

d <- d[which(!is.na(d$stroke_size)),]
dim(d)

d$stroke_size
d$stroke_size <- scale(d$stroke_size)
d$p <- d$cells / d$total
d$logitp <- logit(d$cells / d$total)

#write.csv(d,"/Users/micheldelange/Documents/mice/data/d.csv")



logit <- function(x) {
  return(log(x/(1-x)))
}

############################ linear models ##########
d$p <- d$cells / d$total
d$logitp <- logit(d$cells / d$total)

par(mfrow=c(2,1))
hist(d$p,30)
hist(d$logitp,30)


m0 <- lmer(p ~ 1   + (1 | mouse) , data=d, REML=F  )
summary(m0)
AIC(m0)
r <- residuals(m0)
p <- predict(m0)
plots()

m0 <- lmer(p ~ 1   + (1 | mouse) , data=d, REML=T )
ds <- numeric()
for (i in 1:50) {
  sim <- simulate(m0,nsim=1)
  m3 <- lmer(p ~ 1   + (1 | mouse) , data=d, REML=T  )
  ds[i] <- REMLcrit(m3)
}
plot(ds)
abline(h=REMLcrit(m0),col='red')



plots <- function() {
  par(mfrow=c(3,2))
  plot(r~p,main="residuals ~ fitted")
  plot(r~d$stroke_size)
  boxplot(r~d$female)
  boxplot(r~d$type)
  boxplot(r~d$period)
}

plots()


##==============I am using this model ==============
## AIC says no evidence for interactions

m1 <- lmer(logitp ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d, REML=F  )
summary(m1)
AIC(m1)
r <- residuals(m1)
p <- predict(m1)
plots()
qqnorm(residuals(m1,scale=T))
abline(0,1,col='red')

##==================without sex ================================

m1 <- lmer(logit(p) ~ 1 + stroke_size +  type +  period  + (1 | mouse) , data=d, REML=F )
summary(m1)
AIC(m1)







m2 <- lmer(p ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d, REML=F  )
summary(m2)
AIC(m2)
r <- residuals(m2)
p <- predict(m2)
plots()
anova(m2,m1)

library(lmerTest)
m2 <- lmer(p ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d)
summary(m2)
ds <- numeric()
dev <- REMLcrit(m2)

for (i in 1:50) {
  sim <- simulate(m2,nsim=1)
  m3 <- lmer(unlist(sim) ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d  )
  ds[i] <- REMLcrit(m3)
}
plot(ds)
abline(h=dev,col='red')


##
## arcsin transformation
##

d$p2 <- asin(sqrt(d$p))
m3 <- lmer(p2 ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d , REML=F )
summary(m3)
AIC(m3)
r <- residuals(m3)
p <- predict(m3)
plots()

m3 <- lmer(prop2 ~ 1 + stroke_size  + female * type +  period  + (1 | mouse) , data=d, REML=F  )
summary(m3)
AIC(m3)
r3 <- residuals(m3)
p3 <- predict(m3)
plot(r3~p3,main="residuals ~ fitted")
boxplot(r3~d$female)
boxplot(r3~d$type)
boxplot(r3~d$period)


m4a<- lmer(logitp ~ 1 +  stroke_size          + type *  period  + (1 | mouse) , data=d, REML=F)
m4b<- lmer(logitp ~ 1 +  stroke_size + female + type *  period  + (1 | mouse) , data=d, REML=F  )
anova(m4a,m4b)

AIC(m4a)
AIC(m4b)

summary(m4)


r4 <- residuals(m4)
p4 <- predict(m4)
plot(r4~p4,main="residuals ~ fitted")
boxplot(r4~d$female)
boxplot(r4~d$type)
boxplot(r4~d$period)

step(m3)







############ GLM ###########




m1 <- glm(cells/total ~ 1 + female + type +  period   , family="quasipoisson", data=d  )
summary(m1)

m2 <- glm(cbind(cells,total) ~ 1 + female + type + period, family="quasibinomial", data=d)
summary(m2)
exp(confint(m2))

qqnorm(rstudent(m1))
abline(0,1,col='red')

qqnorm(rstudent(m2))
abline(0,1,col='red')

est <- exp(summary(m1)$coefficients[, 1])
low <- exp(summary(m1)$coefficients[, 1]-2*summary(m1)$coefficients[, 2])
upp <- exp(summary(m1)$coefficients[, 1]+2*summary(m1)$coefficients[, 2])

fem  <- c(est[2],low[2],upp[2])
plt  <- c(est[3],low[3],upp[3])
pou  <- c(est[4],low[4],upp[4])
pe3  <- c(est[5],low[5],upp[5])
pe10 <- c(est[6],low[6],upp[6])
pe17 <- c(est[7],low[7],upp[7])
pe24 <- c(est[8],low[8],upp[8])
pe31 <- c(est[9],low[9],upp[9])




## mixed model
library(MASS)
m3 <- glmmPQL(cells/total ~ 1 + female * type + period, 
        random=~1|mouse,family=quasibinomial,data=d)
pred <- predict(m3,se=TRUE)
summary(m3)
?predict.glmmPQL

m3 <- glmmPQL(cbind(cells,total-cells) ~ 1 + female * type + period, 
              random=~1|mouse,family=quasibinomial,data=d)
pred <- predict(m3,se=TRUE)
summary(m3)


m4 <- glmer(cbind(cells,total-cells) ~ 1 + female * type + period +(1|mouse),
           family=binomial,data=d)
pred <- predict(m3,se=TRUE)
summary(m4)
anova(m4)
m4





par(mfrow=c(1,1))
max_plot = 167
plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  low <- exp(pred$fit[i] - 2 * pred$se.fit[i])
  upp <- exp(pred$fit[i] + 2 * pred$se.fit[i])
  segments(x0=i,y0=low,x1=i,y1=upp,col='black')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
}





library(lme4)
myglm3 <- glmer(cbind(R_long,N_long) ~ 1 + female_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long + (1 | mouse_long), family="binomial")
summary(myglm3)
pre3 <- predict(myglm3)
max_plot = 167
plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  points(i,expit(pre3[i]),col='black',pch='x')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
  segments(i,R_long[i]/N_long[i],i,expit(pre3[i]))
}




library(lme4)
myglm3c <- glm(cbind(R_long,N_long) ~ 1 + female_long + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")
myglm3d <- glm(cbind(R_long,N_long) ~ 1 + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")

summary(myglm3c)
summary(myglm3d)
#(763946-8*2)/5106 + 2 * (8+1)
#1] 167.6142
#> (754808-9*2)/5106 + 2 * (9+1)
#[1] 167.8241

myglm3a <- glm(cbind(R_long,N_long) ~ 1 + female_long + type_1_long + type_2_long + time_3_long + time_10_long
                + time_17_long + time_24_long + time_31_long, family="quasibinomial")
myglm3b <- glm(cbind(R_long,N_long) ~ 1 + type_1_long + type_2_long + time_3_long + time_10_long
              + time_17_long + time_24_long + time_31_long, family="quasibinomial")


anova(myglm3a,myglm3b,test='F')


s <- rstudent(myglm3)
qqnorm(s) 
abline(0,5601^0.5,col='red')
summary(myglm3)
hist(s,probability = T,30)
plot(density(s))
### david

fits<-fitted(myglm)

plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  low <- qbinom(0.025,N_long[i],fits[i])/N_long[i]
  upp <- qbinom(0.975,N_long[i],fits[i])/N_long[i]
  segments(x0=i,y0=low,x1=i,y1=upp,col='black')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
}





