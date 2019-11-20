mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]


# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6



######### THIRD LOT ######################

cell_type <- "CD11bpos" 
N_0 <- mice$Live_cells_0
R_0 <- mice$CD11bpos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD11bpos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD11bpos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD11bpos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD11bpos_24
N_31 <- mice$Live_cells_31
R_31 <- mice$CD11bpos_31
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice6.r")
