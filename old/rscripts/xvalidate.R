getPredictions(m,i) {}
getErrs() {}


fitModel <- function(model,df,parameters) {
  getSamples(model,df,parameters)
}

samples <- fitModel(mouseModel,mice,parameters) 

allMice <- createList(mice)

withinSampleSSE <- function(model,df.mice,allMice,modelName) {
  # fit the model, with one mouse left out.
  samples <- fitModel(model,createList(df.mice),parameters) 
  # mean of posterior predictive for this mouse type, and time
  sse <- 0
  getLPred <- getLPFactory(modelName)
  for (mouseId in 1:dim(df.mice)[1]) {
    for (period in seq(0:4)) {
      expected <- mean(getLPred(s = samples,type = mice[mouseId,"type"],period = period))
      # get actual 
      observed <- getObserved(allMice,mouseId,period)
      observed.lp <- expit(observed)
      sse <- sse + (expected-observed.lp)^2
    }
  }
  return(sse)
}

Sumofsquares  <- withinSampleSSE(mouseModel,  mice,allMice,"mouseModel")



# get SSE by leaving one mouse out, fitting the model
#   and comparing mean posterior for the mouse left out
#   with the actual proportion of cells, 
#   at the linear predictor scale.
#
#   input: df.mice   data.frame with all mice
#          allMice   list version of same
#          mouseId   the mouse to be left out.
#   returns: vector of (SSE), 1 number for each period. could be NA
#
getSSE <- function(mouseId,model,df.mice,allMice,modelName) {

  # mice minus one mouse
  mice_1 <- df.mice[-mouseId,]
  # fit the model, with one mouse left out.
  samples <- fitModel(model,createList(mice_1),parameters) 
  # mean of posterior predictive for this mouse type, and time
  sse <- vector()
  getLPred <- getLPFactory(modelName)
  for (period in seq(0:4)) {
    expected <- mean(getLPred(s = samples,type = mice[mouseId,"type"],period = period))
    # get actual 
    observed <- getObserved(allMice,mouseId,period)
    observed.lp <- expit(observed)
    sse <- c(sse ,(expected-observed.lp)^2)
  }
  return(sse)
}

getObserved <- function(allMice,mouse.id,period) {
  df <- data.frame(allMice)
  allMice$R[which(allMice$mouse == mouse.id & allMice$period == period) ] /
    allMice$Ncells[which(allMice$mouse == mouse.id & allMice$period == period)] 
}

#getSSE(mouseId = 8,mice,allMice.model='mouseModel')

interestingMice <- which(!mice$gender & !mice$type=="PLT")
mm <- data.frame(interestingMice)
dim(mm)

{
Sumofsquares  <- apply(mm,1, FUN=getSSE, mouseModel,  mice,allMice,"mouseModel")
# no interactions
Sumofsquares2 <- apply(mm,1, FUN=getSSE, mouseModel2, mice,allMice,"mouseModel2")
# only depends on type
Sumofsquares3 <- apply(mm,1, FUN=getSSE, mouseModel3, mice,allMice,"mouseModel3")


ss <- mean(unlist(Sumofsquares))
ss2 <- mean(unlist(Sumofsquares2))
ss3 <- mean(unlist(Sumofsquares3))
}
ss
ss2
ss3

# no time:type interactions

mouseModel2 <- function()
{
  for ( i in 1 : N) {
    R[i] ~ dbin(p[mouse[i], period[i]],Ncells[i])
    b[mouse[i], period[i]] ~ dnorm(0.0, tau); 
    logit(p[mouse[i], period[i]]) <- alpha0 + alpha1 * gender[i] +
      alpha2 * type_1[i] +
      alpha3 * type_2[i] +
      alpha4 * time3[i] +
      alpha5 * time10[i] +
      alpha6 * time17[i] +
      alpha7 * time24[i] +
      alpha8 * time31[i] +
      b[mouse[i],period[i]]
  }
  
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  alpha3 ~ dnorm(0.0,1.0E-6)
  alpha4 ~ dnorm(0.0,1.0E-6)
  alpha5 ~ dnorm(0.0,1.0E-6)
  alpha6 ~ dnorm(0.0,1.0E-6)
  alpha7 ~ dnorm(0.0,1.0E-6)
  alpha8 ~ dnorm(0.0,1.0E-6) 
  
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}

mouseModel3 <- function()
{
  for ( i in 1 : N) {
    R[i] ~ dbin(p[mouse[i], period[i]],Ncells[i])
    b[mouse[i], period[i]] ~ dnorm(0.0, tau); 
    logit(p[mouse[i], period[i]]) <- alpha0 + alpha1 * gender[i] +
      alpha2 * type_1[i] +
      alpha3 * type_2[i] +
      b[mouse[i],period[i]]
  }
  
  # priors
  alpha0 ~ dnorm(0.0,1.0E-6)
  alpha1 ~ dnorm(0.0,1.0E-6)
  alpha2 ~ dnorm(0.0,1.0E-6)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100)
}






# periods <- seq(0,4)
# expand.grid(maleMice,periods)


getSSE(interestingMice,mouseModel, allMice,8)


# xvalidate <- function(df) {
#   for (i in 1:dim(df)[1]) {
#     dropm <- df[which(df$)]
#     m <- fitModel(df[-i,],model='?')
#     getPredictions(m,i)
#     getErrs(i)
#   }