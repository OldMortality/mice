createMiceDF <- function(celltype) {

  setwd("~/Documents/mice")
  source("createMatrix.R")
  
  mice <- read.csv("data/mice.csv",header=TRUE)[1:40,]
  # file with stroke size
  mice.s <- read.csv("data/mice2.csv",header=TRUE)[1:40,]
  stroke_size <- scale(mice.s$stroke_size)
  
  
  
  mice$female <- grepl('female',mice$Animal_type)
  mice$type_1 <- grepl('PLT',mice$Animal_type)
  mice$type_2 <- grepl('Pound',
                       mice$Animal_type)
  mice$type <- 'aWT'
  mice[which(mice$type_1),"type"] <- 'PLT'
  mice[which(mice$type_2),"type"] <- 'T2'
  
  # number of mice
  N <- dim(mice)[1]
  # number of time periods (0,3,10,17,24)
  # these are labeled 1 through to 5
  periods <- 5
  
  mice <- mice[-which(mice$type=='PLT'),]
  table(mice$type)
  
  
  female <- grepl('female',mice$Animal_type)
  # type_1 <- grepl('PLT',mice$Animal_type)
  # type_2 <- grepl('Pound',mice$Animal_type)
  
  
  
  
  #celltype <- "CD11cpos"
  
  N_matrix <- getNMatrix(mice,celltype)
  R_matrix <- getRMatrix(mice,celltype)
  
  
  
  # long notation (i.e. observations as rows)
  per_long <- vector()
  mouse_long <- vector()
  female_long <- vector()
  type_long <- vector()
  
  N_long <- vector()
  R_long <- vector()
  stroke_size_long <- vector()
  
  counter = 0
  skipped = 0
  for (i in 1: N) {
    for (j in 1:5) {
      if (is.na(N_matrix[i,j])) {
        skipped <- skipped + 1
        next
      } else {
        counter <- counter + 1
      }
      # periods run from [1,5] ~ [0,3,10,17,24]
      per_long[counter]    <- j
      mouse_long[counter]  <- str_c('mouse',i,sep='')
      female_long[counter] <- female[i]
      type_long[counter] <- mice$type[i]
      stroke_size_long[counter] <- stroke_size[i]
      N_long[counter]      <- N_matrix[i,j]
      R_long[counter]      <- R_matrix[i,j]
    }
  }
  
  
  
  
  mice2 <- data.frame(
    y = R_long / N_long,
    period = factor(per_long),
    mouse = mouse_long,
    female = female_long,
    type = type_long,
    stroke_size = stroke_size_long
  )
  
  head(mice2)
  
  mice2$type <- factor(mice2$type)
  table(mice2$type)
  return(mice2)
}
