identity <- function(x) {
  return(x)
}

waldInterval <- function(model,Z,FUN=exp) {
  s <- summary(model)
  df <- data.frame(estimate = s$coefficients[,1],
                   lower    = s$coefficients[,1] - 1.96 * s$coefficients[,2],
                   upper    = s$coefficients[,1] + 1.96 * s$coefficients[,2]) %>% FUN 

  return(df[-1,] %>% round(2))   # remove intercept
}

waldIntervalBeta <- function(model,Z,FUN=exp) {
  s <- summary(model)$coefficients$cond
  df <- data.frame(estimate = s[,1],
                   lower    = s[,1] - 1.96 * s[,2],
                   upper    = s[,1] + 1.96 * s[,2]) %>% FUN 
  df$pvalue <- s[,4]
  return(df[-1,] %>% round(2))   # remove intercept
}
