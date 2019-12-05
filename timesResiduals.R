timesResiduals <- function(r,df) {
  df <- mice.df
  r <- residuals(m)
  mice <- unique(df$mouse)
  par(mfrow=c(1,1))
  for (t in 2:2) {
    plot('',xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),main=paste('t='t))
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
