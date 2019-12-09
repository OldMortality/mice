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

celltypes <- getcelltypes(mice)  # from createMiceDF.R
celltype <- "CD19pos_B220"
celltype <- celltypes[2]


#
# Creates a dataset with proportions of cells for the given type
#
#
 create.df <- function(df,celltype) {
   
   NMatrix = getNMatrix(df=df,celltype=celltype)
   RMatrix = getRMatrix(df=df,celltype=celltype)
   
    result <- data.frame(
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

   result$type <- 'WT'
   result[which(result$type_1),"type"] <- 'PLT'
   result[which(result$type_2),"type"] <- 'Pound'
      result <- result[ , !(names(result) %in% c("type_1","type_2"))]
   result$y0 <- result$R_0/result$N_0
   result$y3 <- result$R_3/result$N_3
   result$y10 <- result$R_10/result$N_10
   result$y17 <- result$R_17/result$N_17
   result$y24 <- result$R_24/result$N_24
   return(result)
 }


# subset of mice df, for given celltype. Contains proportions R/N
# df.props <- create.df(df = mice,celltype)
# dim(df.props)
# head(df.props)


# return proportions for a give mouse row
getPropsByMouse <- function(df,id) {
  return(df[id,c("y0","y3","y10","y17","y24")])
}


# returns type of given mouse row
getMouseType <- function(df,mouse) {
  return(df$type[mouse])
}

#getMouseType(df.props,6)


# adds mouse to the plot
#  mouse is the row in the df
#  df would be df.props
addMouseToPlot <- function(mouse,df) {
  col <- 'black'
  pch <- 20
  if (getMouseType(df,mouse) == "Pound") {
    col <- 'red' 
  }
  pts <- getPropsByMouse(df,mouse)
  points(0,pts[1,1],col=col,pch=pch)
  points(3,pts[1,2],col=col,pch=pch)
  points(10,pts[1,3],col=col,pch=pch)
  points(17,pts[1,4],col=col,pch=pch)
  points(24,pts[1,5],col=col,pch=pch)
  lines(x=c(0,3,10,17,24),y=pts,col=col,pch=pch)
  print(paste(mouse,col,round(pts,2)))

}


# for (mouse in maleWT.mice) {
#   print(round(getPropsByMouse(df.props,mouse),2))
# }
# for (mouse in malePound.mice) {
#   print(round(getPropsByMouse(df.props,mouse),2))
# }


x <- c(0,3,10,17,24)

addMiceToPlot <- function(mice,df) {
  for (i in 1:length(mice)) {
    addMouseToPlot(mice[i],df)
  }
   
}

addMeansToPlot <- function(df,theGroup,col) {
  lines(c(0,3,10,17,24),apply(getPropsByMouse(df,theGroup),2,mean,na.rm=T),col=col,lty=10,lwd=4)
}


#male.mice <- df.props[which(!df.props$female),]
maleWT.mice <- 
  which(!df.props$female & df.props$type == "WT")
malePound.mice <- 
  which(!df.props$female & df.props$type == "Pound")


# female.mice <- df.props[which(df.props$female),]
# femaleWT.mice <- 
#   which(df.props$female & df.props$type == "WT")
# femalePound.mice <- 
#   which(df.props$female & df.props$type == "Pound")
# 

createDataplot <- function(g1,g2,main,df) {
  
    ylim <- max(df[c(g1,g2),c("y0","y3","y10","y17","y24")] ,na.rm=T) + 0.05
    plot('',xlim=c(0,24),ylim=c(0,ylim),xaxt='n',yaxt='n',
         main=main,xlab='day')
    xbreaks <- c(0,3,10,17,24)
    abline(v=xbreaks,lty=3,col='grey')
    abline(h=seq(0,1,0.1),lty=3,col='grey')
    axis(side = 1, at=xbreaks)
    # ? this does not work
    axis(side = 2, at=c(0.1,0.2,0.3,0.4),labels=c(0.1,0.2,0.3,0.4))
    
    addMiceToPlot(g1,df)
    addMiceToPlot(g2,df)
  
  
  addMeansToPlot(df,g1,'black')
  addMeansToPlot(df,g2,'red')
  
}

par(mfrow=c(1,1))


celltype = celltypes[2] 
celltype

round(getPropsByMouse(df=create.df(df = mice,celltype=celltypes[2]),maleWT.mice),2)
round(getPropsByMouse(df=create.df(df = mice,celltype=celltypes[2]),malePound.mice),2)

par(mfrow=c(1,1))
createDataplot(g1 = maleWT.mice,
               g2 = malePound.mice,
               main = str_c(celltype),
               df = create.df(df = mice,celltype)
                              
               )

create.df(df = mice,
               NMatrix = getNMatrix(df=mice,celltype),
               RMatrix = getRMatrix(df=mice,celltype))

dim(dfx)

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
  g1 <- getPropsByMouse(df.props,malePound.mice)[,per]
  g2 <- getPropsByMouse(df.props,maleWT.mice)[,per]
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
  g1 <- getPropsByMouse(df.props,group)[,per]
  g2 <- getPropsByMouse(df.props,group)[,1]
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

