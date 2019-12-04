#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#
rm(list = ls())
library(stringr)


library(lme4)
# one or the other
#detach(lme4)
#library(glmmTMB)


setwd("~/Documents/mice")
source('createMatrix.R')

mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]


# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

celltype <- "CD19pos_B220"

# getColName <- function(celltype,day) {
#   return(paste(celltype,"pos_",day,sep=''))
# }
# 
# getColName(celltype,0)
# 
# 
# 
# create.df <- function(df,celltype) {
#    df.mice <- data.frame(
#     name = mice$Data.Set,
#     animal_type = mice$Animal_type,
#     N_0 = mice$Live_cells_0,
#     R_0 = mice[,getColName(celltype,0)],
#     N_3 = mice$Live_cells_3,
#     R_3 = mice[,getColName(celltype,3)],
#     N_10 = mice$Live_cells_10,
#     R_10 = mice[,getColName(celltype,10)],
#     N_17 = mice$Live_cells_17,
#     R_17 = mice[,getColName(celltype,17)],
#     N_24 = mice$Live_cells_24,
#     R_24 = mice[,getColName(celltype,24)],
#     N_31 = mice$Live_cells_31,
#     R_31 = mice[,getColName(celltype,31)],
#     female = grepl('female',mice$Animal_type),
#     type_1 = grepl('PLT',mice$Animal_type),
#     type_2 = grepl('Pound',mice$Animal_type)
#   )
# 
#   df.mice$type <- 'WT'
#   df.mice[which(df.mice$type_1),"type"] <- 'PLT'
#   df.mice[which(df.mice$type_2),"type"] <- 'Pound'
#   df.mice <- df.mice[ , !(names(df.mice) %in% c("type_1","type_2"))]
#   df.mice$y0 <- df.mice$R_0/df.mice$N_0
#   df.mice$y3 <- df.mice$R_3/df.mice$N_3
#   df.mice$y10 <- df.mice$R_10/df.mice$N_10
#   df.mice$y17 <- df.mice$R_17/df.mice$N_17
#   df.mice$y24 <- df.mice$R_24/df.mice$N_24
#   return(df.mice)
# }
# 
 create.df <- function(df,celltype,NMatrix,RMatrix) {
    df.mice <- data.frame(
     name = mice$Data.Set,
     animal_type = mice$Animal_type,
     N_0 = NMatrix[,1],
     R_0 = RMatrix[,1],
     N_3 = NMatrix[,2],
     R_3 = RMatrix[,2],
     N_10 = NMatrix[,3],
     R_10 = RMatrix[,3],
     N_17 = NMatrix[,4],
     R_17 = RMatrix[,4],
     N_24 = NMatrix[,5],
     R_24 = RMatrix[,5],
     
     female = grepl('female',mice$Animal_type),
     type_1 = grepl('PLT',mice$Animal_type),
     type_2 = grepl('Pound',mice$Animal_type)
   )

   df.mice$type <- 'WT'
   df.mice[which(df.mice$type_1),"type"] <- 'PLT'
   df.mice[which(df.mice$type_2),"type"] <- 'Pound'
      df.mice <- df.mice[ , !(names(df.mice) %in% c("type_1","type_2"))]
   df.mice$y0 <- df.mice$R_0/df.mice$N_0
   df.mice$y3 <- df.mice$R_3/df.mice$N_3
   df.mice$y10 <- df.mice$R_10/df.mice$N_10
   df.mice$y17 <- df.mice$R_17/df.mice$N_17
   df.mice$y24 <- df.mice$R_24/df.mice$N_24
   return(df.mice)
 }


df.mice <- create.df(df = mice,
                     celltype = celltype,
                     NMatrix = getNMatrix(df=mice,celltype),
                     RMatrix = getRMatrix(df=mice,celltype))
dim(df.mice)
head(df.mice)


getPropsByMouse <- function(df,id) {
  return(df[id,c("y0","y3","y10","y17","y24")])
}



getMouseType <- function(df,mouse) {
  return(df$type[mouse])
}
#getMouseType(df.mice,6)

addMouseToPlot <- function(mouse,df) {
  col <- 'black'
  pch <- 20
  if (getMouseType(df.mice,mouse) == "Pound") {
    col <- 'red' 
  }
  pts <- getPropsByMouse(df.mice,mouse)
  points(0,pts[1,1],col=col,pch=pch)
  points(3,pts[1,2],col=col,pch=pch)
  points(10,pts[1,3],col=col,pch=pch)
  points(17,pts[1,4],col=col,pch=pch)
  points(24,pts[1,5],col=col,pch=pch)
  lines(x=c(0,3,10,17,24),y=pts,col=col,pch=pch)

}


# for (mouse in maleWT.mice) {
#   print(round(getPropsByMouse(df.mice,mouse),2))
# }
# for (mouse in malePound.mice) {
#   print(round(getPropsByMouse(df.mice,mouse),2))
# }


x <- c(0,3,10,17,24)

addMiceToPlot <- function(mice,df) {
  for (i in 1:length(mice)) {
    addMouseToPlot(mice[i],df)
  }
   
}

addMeansToPlot <- function(df,theGroup,col) {
  lines(x,apply(getPropsByMouse(df.mice,theGroup),2,mean,na.rm=T),col=col,lty=10,lwd=4)
}


male.mice <- df.mice[which(!df.mice$female),]
maleWT.mice <- 
  which(!df.mice$female & df.mice$type == "WT")
malePound.mice <- 
  which(!df.mice$female & df.mice$type == "Pound")


female.mice <- df.mice[which(df.mice$female),]
femaleWT.mice <- 
  which(df.mice$female & df.mice$type == "WT")
femalePound.mice <- 
  which(df.mice$female & df.mice$type == "Pound")


createDataplot <- function(g1,g2,main,df,ylim) {
  
  { 
    
    plot('',xlim=c(0,24),ylim=c(0,ylim),xaxt='n',yaxt='n',
         main=main,xlab='day')
    xbreaks <- c(0,3,10,17,
                 24)
    abline(v=xbreaks,lty=3,col='grey')
    abline(h=seq(0,1,0.1),lty=3,col='grey')
    axis(side = 1, at=xbreaks)
    # ? this does not work
    axis(side = 2, at=c(0.1,0.2,0.3,0.4),labels=c(0.1,0.2,0.3,0.4))
    
    addMiceToPlot(g1,df)
    addMiceToPlot(g2,df)
  }
  
  addMeansToPlot(df,g1,'black')
  addMeansToPlot(df,g2,'red')
  
}

par(mfrow=c(1,1))

celltype = "CD11cpos"
createDataplot(maleWT.mice,malePound.mice,main=str_c(celltype," male"),mice.df,ylim=0.1)


CD11cpos

#createDataplot(femaleWT.mice,femalePound.mice,main=str_c(celltype," female"))



### t.tests
# per here is the column, so
#  1 = time0
#  2 = time3
#  3 = time10
#  4 = time17
#  5 = time24
do.ttest <- function(per) {
  g1 <- getPropsByMouse(df.mice,malePound.mice)[,per]
  g2 <- getPropsByMouse(df.mice,maleWT.mice)[,per]
  t.test(g1,g2)
}


for (i in 1:5) {
  print(i)
  t <- do.ttest(i)
  print(t)
}

par(mfrow=c(1,1))
main = paste(celltype,' t-tests: CIs mn(Pound) - mn(WT')
plot('',xlim=c(0,24),ylim=c(-1,1),xaxt='n',
     main=main,xlab='day')
xbreaks <- c(0,3,10,17,
             24)
abline(v=xbreaks,lty=3,col='grey')
abline(h=seq(-1,1,0.1),lty=3,col='grey')
axis(side = 1, at=xbreaks)
for (i in 1:5) {
  t <- do.ttest(i)
  ci <- t$conf.int
  segments(xbreaks[i],ci[1],xbreaks[i],ci[2])
}
abline(h=0,col='red')


# paired t-test for one type of mouse, period compared to period 0
doPairedttest <- function(per,group) {
  g1 <- getPropsByMouse(df.mice,group)[,per]
  g2 <- getPropsByMouse(df.mice,group)[,1]
  t.test(g1,g2,paired=T)
}

group = maleWT.mice
doPairedttest(2,group)$p.value
doPairedttest(3,group)$p.value
doPairedttest(4,group)$p.value
doPairedttest(5,group)$p.value
group = malePound.mice
doPairedttest(2,group)$p.value
doPairedttest(3,group)$p.value
doPairedttest(4,group)$p.value
doPairedttest(5,group)$p.value

