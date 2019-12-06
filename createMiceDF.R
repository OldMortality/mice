createMiceDF <- function(celltype) {

  setwd("~/Documents/mice")
  source("createMatrix.R")
  
  mice <- read.csv("data/mice.csv",header=TRUE)[1:40,]
  # file with stroke size
  mice.s <- read.csv("data/mice2.csv",header=TRUE)[1:40,]
  stroke_size <- scale(mice.s$stroke_size)
  
  
  
  mice$female <- grepl('female',mice$Animal_type)
  mice$type_1 <- grepl('PLT',mice$Animal_type)
  mice$type_2 <- grepl('Pound',
                       mice$Animal_type)
  mice$type <- 'aWT'
  mice[which(mice$type_1),"type"] <- 'PLT'
  mice[which(mice$type_2),"type"] <- 'T2'
  
  # number of mice
  N <- dim(mice)[1]
  # number of time periods (0,3,10,17,24)
  # these are labeled 1 through to 5
  periods <- 5
  
  mice <- mice[-which(mice$type=='PLT'),]
  table(mice$type)
  female <- grepl('female',mice$Animal_type)
  N_matrix <- getNMatrix(mice,celltype)
  R_matrix <- getRMatrix(mice,celltype)
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
      if (is.na(N_matrix[i,j])) {
        skipped <- skipped + 1
        next
      } else {
        counter <- counter + 1
      }
      # periods run from [1,5] ~ [0,3,10,17,24]
      per_long[counter]    <- j
      mouse_long[counter]  <- str_c('mouse',i,sep='')
      female_long[counter] <- female[i]
      type_long[counter] <- mice$type[i]
      stroke_size_long[counter] <- stroke_size[i]
      N_long[counter]      <- N_matrix[i,j]
      R_long[counter]      <- R_matrix[i,j]
    }
  }
  
  
  
  
  mice2 <- data.frame(
    y = R_long / N_long,
    period = factor(per_long),
    mouse = mouse_long,
    female = female_long,
    type = type_long,
    stroke_size = stroke_size_long
  )
  
  head(mice2)
  
  mice2$type <- factor(mice2$type)
  table(mice2$type)
  return(mice2)
}

#
# gets the Wald CI at the linear predictor scale.
#
#
getCI <- function(model,z) {
  s <- summary(model)
  CI_lower <- s$coefficients[,1] - z*s$coefficients[,2]
  CI_upper <- s$coefficients[,1] + z*s$coefficients[,2]
  return(round(cbind(CI_lower,CI_upper),2))
}

# r: a vector of residuals (the real ones)
# sim.r   : matrix, each column is a set of simulated residuals
# xl are the limits on the x-axis
policePlot <- function(r,rsim,xl) {
  
  par(mfrow=c(2,5))
  real.pos <- floor(runif(1,min=1,max=11))
  
  for (i in 1:dim(rsim)[2]) {
    if (i == real.pos) {
      hist(r,main = "" ,15,xlim = xl,xlab='*')
    } else {
      hist(rsim[,i],main = "",15,xlim = xl,xlab="")
    }
  }
  return(real.pos)
}


##  
## plots residuals~fitted, ~ type, ~ time
##
##   input: df=dataframe, fit=fitted, r=residuals,
##      qq should qqplot be drawn
##
plots <- function(df,fit,r,qq=NULL) {
  par(mfrow=c(2,2))
  col <- rep('black',dim(df)[1])
  col[which(df$type=='T2')] <- 'red'
  plot(r~fit,main="residuals ~ fitted",col=col)
  abline(h=0,col='red')
  
  boxplot(r~df$type,main='residuals~type')
  boxplot(r~df$period,main='residuals~time')
  
  if (qq) {
    qqnorm(scale(r))
    abline(0,1,col='red')
  }
  
}

logit <- function(x) {
  stopifnot(x<1,x>0)
  return(log(x/(1-x)))
}

invlogit <- function(x) {
  return(1/(1+exp(-x)))
}


##
## Plots to show autocorrelation of residuals
##
timesResiduals <- function(r,df,main) {
  df <- mice.df
  mice <- unique(df$mouse)
  par(mfrow=c(2,3))
  for (t in 2:5) {
    plot('',xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),
         main=paste(main,'t=',t))
    abline(h=0,col='red')
    abline(v=0,col='red')
    for (mouse in mice) {
      x = which(df$mouse==mouse & df$period == t)
      y = which(df$mouse==mouse & df$period == t-1)
      if (length(x)==1 & length(y)==1)
        points(r[x],r[y],pch='x')
    }
  }
  
}

