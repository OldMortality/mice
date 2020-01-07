identity <- function(x) {
  return(x)
}

addStar <- function(df) {
  df$signif <- '*'
  df[which(df$lower < 1 & df$upper > 1),"signif"] <- ""
  df <- df[-1,]
  df[,c(1,2,3)] <- round(df[,c(1,2,3)],2)
  return(df)
}

waldInterval <- function(model,Z,FUN=exp) {
  s <- summary(model)
  df <- data.frame(estimate = s$coefficients[,1],
                   lower    = s$coefficients[,1] - Z * s$coefficients[,2],
                   upper    = s$coefficients[,1] + Z * s$coefficients[,2]) %>% 
    FUN %>% 
    addStar %>%
    return
  #return(df[-1,c(1,2,3)] %>% round(2))   # remove intercept
  
}

waldInterval.nlme <- function(model,Z,FUN=exp) {
  s <- summary(model)$tTable
  df <- data.frame(estimate = s[,1],
                   lower    = s[,1] - Z * s[,2],
                   upper    = s[,1] + Z * s[,2]
                   ) %>% 
    FUN %>% 
    addStar %>%
    cbind(s[-1,5])
  
    colnames(df)[5] <- 'p.value'
    return(df)
  #return(df[-1,c(1,2,3)] %>% round(2))   # remove intercept
  
}

waldIntervalBeta <- function(model,Z,FUN=exp) {
  s <- summary(model)$coefficients$cond
  df <- data.frame(estimate = s[,1],
                   lower    = s[,1] - Z * s[,2],
                   upper    = s[,1] + Z * s[,2]) %>% FUN 
  df$pvalue <- s[,4]
  return(df[-1,] %>% round(2))   # remove intercept
}
