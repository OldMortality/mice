mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]


# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

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
N_31 <- mice$Live_cells_31
R_31 <- mice$CD19pos_B220pos_31

source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice6.r")
dim(N_0)

mice$Animal_type <- c(rep("female PLT",7),
                      rep("female Pound",7),
                      rep("female WT",6),
                      rep("male Pound",7),
                      rep("male PLT",7),
                      rep("male WT",6))
N_0 <- rep(10000,40)
N_3 <- rep(10000,40)
N_10 <- rep(10000,40)
N_17<- rep(10000,40)
N_24 <- rep(10000,40)
N_31 <- rep(10000,40)
R_0 <- round( c(runif(20,4000,5000),runif(20,1000,2000)) )
R_3 <- round( c(runif(20,4000,5000),runif(20,1000,2000)) )
R_10 <- round( c(runif(20,4000,5000),runif(20,1000,2000)) )
R_17 <- round( c(runif(20,4000,5000),runif(20,1000,2000)) )
R_24 <- round( c(runif(20,4000,5000),runif(20,1000,2000)) )
R_31 <- c(rep(5000,20),rep(2000,20))

