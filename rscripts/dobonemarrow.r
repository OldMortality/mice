 
  
  d <- data.frame(N,R,female,type) 
  m <- glm(cbind(R,N) ~ female + type,family="quasibinomial",data=mice_bm)
  summary(m)
  pred <- predict(m,se.fit=T)
  r <- residuals(m)
  rs <- rstandard(m)
  #qqnorm(rs)
  #abline(0,1,col='red')
  #plot(r ~ pred$fit)
  #plot(rs~pred$fit)
  #plot(r~female)
  #plot(r~type)
  
  
  m.grid <- ref.grid(m)
  #m.grid
  #summary(m.grid)
  
  
  lsm <- lsmeans(m.grid,'female',by="type")
  #str(lsm)
  plot(lsm,'female',by="type",type="response",main=celltype)
 