mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice2.csv",header=TRUE)[1:40,]

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}
logit <- function(x) {
  return(x/(1-x))
}

# number of mice
N <- dim(mice)[1]
mice <- mice[-which(mice$Stroke.volume..mm3. == "died"),]
N <- dim(mice)[1]

# number of time periods (0,3,10,17,24)
periods <- 5

######### FIRST LOT ######################
cell_type <- "CD19pos_B220pos"
N_0 <- mice$Live_cells_0
R_0 <- mice$CD19pos_B220pos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD19pos_B220pos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD19pos_B220pos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD19pos_B220pos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD19pos_B220pos_24

ss <- as.numeric(as.character(mice$Stroke.volume..mm3.))
stroke_size <- scale(as.numeric(as.character(mice$Stroke.volume..mm3.)))[1:35]

### some analysis of stroke size, unscaled ###

s <- data.frame
plot(ss~type)
length(ss)
length(type)
####


var(stroke_size)
### scale
{
  scaleIt <- function(R,N) {
    return(round(R * 1000 / N))
  }
  scaleN <- function(x) {
    for (i in 1:length(x)) {
      if (!is.na(x[i])) {
        x[i] <- 1000
      }
    }
    return(x)
  }
  
}

#R_0 <- scaleIt(R_0,N_0)
#R_3 <- scaleIt(R_3,N_3)
#R_10 <- scaleIt(R_10,N_10)
#R_17 <- scaleIt(R_17,N_17)
#R_24 <- scaleIt(R_24,N_24)
#N_0 <- scaleN(N_0)
#N_3 <-scaleN(N_3)
#N_10 <-scaleN(N_10)
#N_17 <-scaleN(N_17)
#N_24 <-scaleN(N_24)

source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")



######### SECOND LOT ####### N_0 <- mice$CD11cpos_0
cell_type = "CD11cpos"
R_0 <- mice$CD11cpos_0
R_3 <- mice$CD11cpos_3
R_7 <- mice$CD11cpos_7
R_10 <- mice$CD11cpos_10
R_17 <- mice$CD11cpos_17
R_24 <- mice$CD11cpos_24

# N's are the same as CD19B220

source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")

######## Third lot CD11bpos #########

cell_type = "CD11bpos"
R_0 <- mice$CD11bpos_0
R_3 <- mice$CD11bpos_3
R_7 <- mice$CD11bpos_7
R_10 <- mice$CD11bpos_10
R_17 <- mice$CD11bpos_17
R_24 <- mice$CD11bpos_24

# N's are the same as CD19B220
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")

######## 4th lot "mMDCSC"

cell_type = "mMDCSCs"
R_0 <- mice$mMDSCs_0
R_3 <- mice$mMDSCs_3
R_10 <- mice$mMDSCs_10
R_17 <- mice$mMDSCs_17
R_24 <- mice$mMDSCs_24
N_0 <- mice$CD11bpos_0
N_3 <- mice$CD11bpos_3
N_10 <- mice$CD11bpos_10
N_17 <- mice$CD11bpos_17
N_24 <- mice$CD11bpos_24
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")

### 4b mMDCSC as proportion of live cells.
N_0 <- mice$Live_cells_0
N_3 <- mice$Live_cells_3
N_10 <- mice$Live_cells_10
N_17 <- mice$Live_cells_17
N_24 <- mice$Live_cells_24

source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")
######## 5th lot "mMDCSC"

cell_type = "gmMDCSs"
R_0 <- mice$gMDSCs_0
R_3 <- mice$gMDSCs_3
R_10 <- mice$gMDSCs_10
R_17 <- mice$gMDSCs_17
R_24 <- mice$gMDSCs_24
N_0 <- mice$CD11bpos_0
N_3 <- mice$CD11bpos_3
N_10 <- mice$CD11bpos_10
N_17 <- mice$CD11bpos_17
N_24 <- mice$CD11bpos_24
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")


### 5b mMDCSC as proportion of live cells.
N_0 <- mice$Live_cells_0
N_3 <- mice$Live_cells_3
N_10 <- mice$Live_cells_10
N_17 <- mice$Live_cells_17
N_24 <- mice$Live_cells_24
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")

######## 6th lot "intMDCSC"

cell_type = "intMDCSs"
R_0 <- mice$intMDSCs_0
R_3 <- mice$intMDSCs_3
R_10 <- mice$intMDSCs_10
R_17 <- mice$intMDSCs_17
R_24 <- mice$intMDSCs_24
N_0 <- mice$CD11bpos_0
N_3 <- mice$CD11bpos_3
N_10 <- mice$CD11bpos_10
N_17 <- mice$CD11bpos_17
N_24 <- mice$CD11bpos_24
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")

### 6b intMDCSC as proportion of live cells.
N_0 <- mice$Live_cells_0
N_3 <- mice$Live_cells_3
N_10 <- mice$Live_cells_10
N_17 <- mice$Live_cells_17
N_24 <- mice$Live_cells_24
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")


### 7 MHC
mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice2.csv",header=TRUE)

### CD19pos_B220pos_MHCIIpos_0
###    We don't have the numbers N, only proportions, so we scale them to N=1000
mice2 <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[41:80,]
length(mice$CD19pos_B220pos_MHCIIpos_0[41:80])
x <- which(mice$Stroke.volume..mm3. == "died")
mice2 <- mice2[-x,]
dim(mice2)
R_0 <- 1000 * mice2$CD19pos_B220pos_MHCIIpos_0
R_3 <- 1000 * mice2$CD19pos_B220pos_MHCIIpos_3
R_10 <- 1000 * mice2$CD19pos_B220pos_MHCIIpos_10
R_17 <- 1000 * mice2$CD19pos_B220pos_MHCIIpos_17
R_24 <- 1000 * mice2$CD19pos_B220pos_MHCIIpos_24
N_0 <- 1000
N_3 <- 1000
N_10 <- 1000
N_17 <- 1000
N_24 <- 1000
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice_B220pos.r")



