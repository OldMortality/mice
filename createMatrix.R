getNMatrix <- function(df,celltype) {
    N_0  <- df$Live_cells_0
    N_3  <- df$Live_cells_3
    N_10 <- df$Live_cells_10
    N_17 <- df$Live_cells_17
    N_24 <- df$Live_cells_24
    return(data.frame(cbind(N_0,N_3,N_10,N_17,N_24)))
  
}



getRMatrix <- function(df,celltype) {
  
  getColname <- function(num,celltype) {
    result <- paste(celltype,'pos_',num,sep='')
    if (celltype %in% c("mMDSCs","gMDSCs","intMDSCs")) {
      result <- result %>% str_remove('pos')
    }
    return(result)
  }
  
  if (celltype %in% 
      c("CD19pos_B220","CD11c","CD19","CD19pos_B220pos_MHCII","CD11b","mMDSCs","gMDSCs","intMDSCs",
        "mMDSCs_MHCII")) {
    days <- c(0,3,10,17,24)
    colnames <- days %>% map(getColname,celltype=celltype) %>% unlist
    result <- df[,colnames]
    
    colnames(result) <- paste("R_",days,sep="")
    return(result)
  }
  
}

#getRMatrix(mice,celltypes[7])
