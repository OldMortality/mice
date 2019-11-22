 
library("lsmeans")
library("effects")


invlogit <- function(x) {
  return(1/(1+exp(-x)))
}

mice_bm <- read.csv("~/Documents/mice/data/BM2.csv",na.strings=c("","NA"),header=TRUE)[1:28,]
mice_bm$Animal_ID[c(1,28)]
dim(mice_bm)
female <- grepl('female',mice_bm$Animal_type)
length(female)

type <- rep("aWT",28)
type[grepl('PLT',mice_bm$Animal_type)] <- "PLT"
type[grepl('Pound',mice_bm$Animal_type)] <- "Pound"
female <- factor(female)
type <- factor(type)




doit <- function(celltype, N, R) {
  
  
  p <- R/N
  logitp <- log(p/(1-p))
  
  
  m2 <- glm(logitp ~ female * type)
  summary(m2)
  
  new.data <- data.frame(female=rep(c("TRUE","FALSE"),3),
                         type=c("aWT","aWT","PLT","PLT","Pound","Pound"))
  
  pre <- predict(m2,newdata <- new.data,se.fit=T)
   
  
  plot(0,0,pch='',xlim=c(0,7),ylim=c(0,0.5),xaxt = 'n',ylab="Prob")
  axis(1,at=1:6,labels=c("WT f","WT m","PLT f","PLT m","Pound f","Pound m"))
  grid()
  for (i in 1:6) {
    segments(i,invlogit(pre$fit[i] - 2 * pre$se.fit[i]),i,invlogit(pre$fit[i] + 2 * pre$se.fit[i]))
  }
  
   
  
  
   
  
  #pre <- predict(m2,se.fit=T)
  ##r <- residuals(m2)
  #qqnorm(rstandard(m2))
  ##abline(0,1,col='red')
  #plot(r)
  #plot(r~female)
  #plot(r~type)
  #plot(r~pre$fit)
  
  
}







##CD19pos_B220pos_0

 

N <- mice_bm$live_cells_0
R <- mice_bm$CD19pos_B220pos_0

doit("Bonemarrow: CD19pos_B220pos_0",N,R)

R <- mice_bm$CD19pos_B220pos_MHCIIpos_0
doit("Bonemarrow: CD19pos_B220pos_0",N,R)
