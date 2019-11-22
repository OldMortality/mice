##
##  Away with all that Bayesian stuff
##  Long data
##  quasibinomial model
##
##

mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]

# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

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

gender <- grepl('female',mice$Animal_type)
type_1 <- grepl('PLT',mice$Animal_type)
type_2 <- grepl('Pound',mice$Animal_type)

N_matrix <- data.frame(cbind(N_0,N_3,N_10,N_17,N_24,N_31))
R_matrix <- data.frame(cbind(R_0,R_3,R_10,R_17,R_24,R_31))
P_matrix <- R_matrix/N_matrix


time_3 <-  c(0,1,0,0,0,0)
time_10 <- c(0,0,1,0,0,0)
time_17 <- c(0,0,0,1,0,0)
time_24 <- c(0,0,0,0,1,0)
time_31 <- c(0,0,0,0,0,1)


 
# long notation (i.e. observations as rows)
per_long <- vector()
mouse_long <- vector()
gender_long <- vector()
type_1_long <- vector()
type_2_long <- vector()
time_3_long <- vector()
time_10_long <- vector()
time_17_long <- vector()
time_24_long <- vector()
time_31_long <- vector()

N_long <- vector()
R_long <- vector()
 

counter = 0
skipped = 0
for (i in 1: N) {
  for (j in 1:6) {
    
    counter <- counter + 1
    if (is.na(N_matrix[i,j])) {
      skipped <- skipped + 1
      next
    }
    per_long    <- c(per_long,j)
    mouse_long  <- c(mouse_long,i)
    gender_long <- c(gender_long,gender[i])
    type_1_long <- c(type_1_long,type_1[i])
    type_2_long <- c(type_2_long, type_2[i])
    N_long      <- c(N_long,N_matrix[i,j])
    R_long      <- c(R_long,R_matrix[i,j])
    time_3_long <- c(time_3_long,time_3[j])
    time_10_long <-c(time_10_long,time_10[j])
    time_17_long <- c(time_17_long,time_17[j])
    time_24_long <- c(time_24_long,time_24[j])
    time_31_long <- c(time_31_long,time_31[j])
  }
  
}

long_data <- data.frame(mouse_long,
                   per_long,
                   gender_long,
                   type_1_long,
                   type_2_long,
                   N_long,
                   R_long,
                   time_3_long,
                   time_10_long,
                   time_17_long,
                   time_24_long,
                   time_31_long)

write.csv(long_data,"/Users/micheldelange/Documents/mice/data/long_data.csv")


############ GLM ###########

myglm <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long, family="quasibinomial" )

myglm2 <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
             + time_17_long + time_24_long + time_31_long
             + gender_long * type_1_long + gender_long * type_2_long 
             , 
             family="quasibinomial" )


summary(myglm)

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

pred <- predict(myglm,se=TRUE)

max_plot = 167
plot(0,0,xlim=c(0,max_plot),ylim=c(0,0.5),ylab='prop alive')
for (i in 1:max_plot) {
  low <- expit(pred$fit[i] - 2 * pred$se.fit[i])
  upp <- expit(pred$fit[i] + 2 * pred$se.fit[i])
  segments(x0=i,y0=low,x1=i,y1=upp,col='black')
}
for (i in 1:max_plot) {
  points(i,R_long[i]/N_long[i],col='red',pch='x')
}

library(lme4)
myglm3 <- glmer(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
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
myglm3c <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")
myglm3d <- glm(cbind(R_long,N_long) ~ 1 + type_1_long + type_2_long + time_3_long + time_10_long
               + time_17_long + time_24_long + time_31_long, family="binomial")

summary(myglm3c)
summary(myglm3d)
#(763946-8*2)/5106 + 2 * (8+1)
#1] 167.6142
#> (754808-9*2)/5106 + 2 * (9+1)
#[1] 167.8241

myglm3a <- glm(cbind(R_long,N_long) ~ 1 + gender_long + type_1_long + type_2_long + time_3_long + time_10_long
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


### beta binomial Bayes


