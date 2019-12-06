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

timesResiduals(residuals(m),df.mice,"logit lm")
