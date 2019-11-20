##
##  Away with all that Bayesian stuff
##  Long data
##  linear models against logit p
##
##
library(lme4)
library(knitr)
library(stringr)

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit <- function(x) {
  return(log(x/(1-x)))
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
 
d <- long_data

colnames(d) <- c("mouse","female","type","period","total","cells","stroke_size")

d <- d[which(!is.na(d$stroke_size)),]
dim(d)

d$stroke_size
d$stroke_size <- scale(d$stroke_size)
d$p <- d$cells / d$total
 
 




plots <- function() {
  par(mfrow=c(3,2))
  plot(r~p,main="residuals ~ fitted")
  plot(r~d$stroke_size)
  boxplot(r~d$female)
  boxplot(r~d$type)
  boxplot(r~d$period)
}

 

############################ linear models ##########
d$p <- d$cells / d$total
 

m0 <- lmer(p ~ 1   + (1 | mouse) , data=d, REML=T )
ds <- numeric()
for (i in 1:50) {
  sim <- simulate(m0,nsim=1)
  m3 <- lmer(p ~ 1   + (1 | mouse) , data=d, REML=T  )
}
r <- residuals(m0)
p <- predict(m0)
plots()



##==============I am using this model ==============
## AIC says no evidence for interactions

m1 <- lmer(logit(p) ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d, REML=F  )
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
r <- residuals(m1)
p <- predict(m1)
plots()
qqnorm(residuals(m1,scale=T))
abline(0,1,col='red')


## police line-up of residuals
par(mfrow=c(5,4))
for (i in 1:19) {
  sim <- unlist(simulate(m1,nsim=1))
  m <- lmer(sim ~ 1 + stroke_size + female + type +  period  + (1 | mouse) , data=d  )
  hist(residuals(m3,scale=T),30)
  
}
hist(residuals(m1,scale=T),30)
abline(h=dev,col='red')


 