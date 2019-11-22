#
#
# 2019. Took this from mice2.rmd.
#       Fits Bayesian model. 
#
#
#library(rjags)
library(R2jags)

macFile <- "/Users/micheldelange/Documents/mice/data/mice.csv"
ubuntuFile <- "~/mice/data/mice.csv"

if (file.exists(macFile)) {
  mice <- read.csv(macFile,header=TRUE)[1:40,]
} else {
  mice <- read.csv(ubuntuFile,header=TRUE)[1:40,]
}

# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

celltype <- "CD19pos_B220"

getColName <- function(celltype,day) {
  return(paste(celltype,"pos_",day,sep=''))
}

getColName(celltype,0)



create.df <- function(df,celltype) {
   df.mice <- data.frame(
    name = mice$Data.Set,
    animal_type = mice$Animal_type,
    N_0 = mice$Live_cells_0,
    R_0 = mice[,getColName(celltype,0)],
    N_3 = mice$Live_cells_3,
    R_3 = mice[,getColName(celltype,3)],
    N_10 = mice$Live_cells_10,
    R_10 = mice[,getColName(celltype,10)],
    N_17 = mice$Live_cells_17,
    R_17 = mice[,getColName(celltype,17)],
    N_24 = mice$Live_cells_24,
    R_24 = mice[,getColName(celltype,24)],
    N_31 = mice$Live_cells_31,
    R_31 = mice[,getColName(celltype,31)],
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


df.mice <- create.df(df = mice,celltype = celltype)
dim(df.mice)
head(df.mice)


getPropsByMouse <- function(df,id) {
  return(df[id,c("y0","y3","y10","y17","y24")])
}




getMouseType <- function(df,mouse) {
  return(df$type[mouse])
}
#getMouseType(df.mice,6)

addMouseToPlot <- function(mouse) {
  print(mouse)
  col <- 'black'
  if (getMouseType(df.mice,mouse) == "Pound") {
    col <- 'red'
  }
  pts <- getPropsByMouse(df.mice,mouse)
  points(0,pts[1,1],col=col)
  points(3,pts[1,2],col=col)
  points(10,pts[1,3],col=col)
  points(17,pts[1,4],col=col)
  points(24,pts[1,5],col=col)
  lines(x=c(0,3,10,17,24),y=pts,col=col)
}

addMiceToPlot <- function(mice) {
  for (i in 1:length(mice)) {
    addMouseToPlot(mice[i])
  }
}

male.mice <- df.mice[which(!df.mice$female),]
maleWT.mice <- 
  which(!df.mice$female & df.mice$type == "WT")
malePound.mice <- 
  which(!df.mice$female & df.mice$type == "Pound")

plot('',xlim=c(0,31),ylim=c(0,0.4))
addMiceToPlot(maleWT.mice)
addMiceToPlot(malePound.mice)






