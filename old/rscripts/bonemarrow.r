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
  
  d <- data.frame(N,R,female,type) 
  m <- glm(cbind(R,N) ~ female + type,family="quasibinomial",data=mice_bm)
  summary(m)
  pred <- predict(m,se.fit=T)
  r <- residuals(m)
  rs <- rstandard(m)
  qqnorm(rs)
  abline(0,1,col='red')
  plot(r ~ pred$fit)
  plot(rs~pred$fit)
  plot(r~female)
  plot(r~type)
  
  
  m.grid <- ref.grid(m)
  m.grid
  summary(m.grid)
  
  
  lsm <- lsmeans(m.grid,'female',by="type")
  str(lsm)
  plot(lsm,'female',by="type",type="response",main=celltype)
}


N <- mice_bm$live_cells_0
R <- mice_bm$CD19pos_B220pos_0
doit("Bonemarrow: CD19pos_B220pos_0",N,R)

R <- mice_bm$CD19pos_B220pos_MHCIIpos_0
doit(celltype="Bonemarrow: CD19pos_B220pos_MHCIIpos_0",N=N,R=R)


# spleen

sp <- read.csv("~/Documents/mice/data/spleen.csv",na.strings=c("","NA"),header=TRUE)[1:27,]
sp <- sp[-7,]
dim(sp)
female <- grepl('female',sp$Animal_type)
length(female)

type <- rep("aWT",dim(sp)[1])
type[grepl('PLT',sp$Animal_type)] <- "PLT"
type[grepl('Pound',sp$Animal_type)] <- "Pound"
female <- factor(female)
type <- factor(type)

N <- sp$live_cells_0
R <- sp$CD19pos_B220pos_0
doit("Spleen: CD19pos_B220pos_0",N,R)
