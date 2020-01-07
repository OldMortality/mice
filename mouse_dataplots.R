#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#
rm(list = ls())
library(stringr)
library(purrr)



setwd("~/Documents/mice")
source('createMatrix.R')
source('createMiceDF.R')

mice <- read.csv("~/Documents/mice/data/mice.csv",header=TRUE)[1:40,]
celltypes <- getcelltypes(mice)


# number of mice
N <- dim(mice)[1]

C <- 1
celltypes <- getcelltypes(mice)  # from createMiceDF.R
celltype <- celltypes[C]
celltype


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



# return proportions for a give mouse row
getPropsByMouse <- function(df,id) {
  return(df[id,c("y0","y3","y10","y17","y24")])
}

# returns type of given mouse row
getMouseType <- function(df,mouse) {
  return(df$type[mouse])
}


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


addMiceToPlot <- function(mice,df) {
  for (i in 1:length(mice)) {
    addMouseToPlot(mice[i],df)
  }
   
}

addMeansToPlot <- function(df,theGroup,col) {
  lines(c(0,3,10,17,24),apply(getPropsByMouse(df,theGroup),2,mean,na.rm=T),col=col,lty=10,lwd=4)
}


df.props <- create.df(df = mice,celltype)
maleWT.mice <-
  which(!df.props$female & df.props$type == "WT")
malePound.mice <-
  which(!df.props$female & df.props$type == "Pound")


getMaleWT <- function(df) {
  which(!df$female & df$type == "WT")
}

getMalePound <- function(df) {
  which(!df$female & df$type == "Pound")
}


createDataplot <- function(main,df) {
  
    maleWT.mice <- getMaleWT(df)
    malePound.mice <- getMalePound(df)
  
    ylim <- max(df[c(maleWT.mice,malePound.mice),c("y0","y3","y10","y17","y24")] ,na.rm=T) * 1.1
    plot('',xlim=c(0,24),ylim=c(0,ylim),xaxt='n',yaxt='n',
         main=main,xlab='day')
    xbreaks <- c(0,3,10,17,24)
    abline(v=xbreaks,lty=3,col='grey')
    abline(h=seq(0,1,0.1),lty=3,col='grey')
    axis(side = 1, at=xbreaks)
    # ? this does not work
    axis(side = 2, at=c(0.1,0.2,0.3,0.4),labels=c(0.1,0.2,0.3,0.4))
    
    addMiceToPlot(maleWT.mice,df)
    addMiceToPlot(malePound.mice,df)
  
  
  addMeansToPlot(df,maleWT.mice,'black')
  addMeansToPlot(df,malePound.mice,'red')
  
}

par(mfrow=c(1,1))
C = 8
celltype <- celltypes[C]
df <- create.df(df = mice, celltype = celltype)
celltype
round(getPropsByMouse(df=df,getMaleWT(df)),2)
round(getPropsByMouse(df=df,getMalePound(df)),2)

par(mfrow=c(1,1))
C = 11
createDataplot(main = str_c(celltypes[C]),
               df = create.df(df = mice,celltypes[C]))

