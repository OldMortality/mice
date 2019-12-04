getNMatrix <- function(df,celltype) {
  if (celltype %in% c("CD19pos_B220","CD11cpos")) {
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
    return(data.frame(cbind(N_0,N_3,N_10,N_17,N_24)))
  }
  
  
}

getRMatrix <- function(df,celltype) {
  if (celltype == "CD19pos_B220") {
    R_0 <- mice$CD19pos_B220pos_0
    R_3 <- mice$CD19pos_B220pos_3
    R_10 <- mice$CD19pos_B220pos_10
    R_17 <- mice$CD19pos_B220pos_17
    R_24 <- mice$CD19pos_B220pos_24
    return(data.frame(cbind(R_0,R_3,R_10,R_17,R_24)))
  }
  if (celltype == "CD11cpos") {
    R_0 <- mice$CD11cpos_0
    R_3 <- mice$CD11cpos_3
    R_7 <- mice$CD11cpos_7
    R_10 <- mice$CD11cpos_10
    R_17 <- mice$CD11cpos_17
    R_24 <- mice$CD11cpos_24
    return(data.frame(cbind(R_0,R_3,R_10,R_17,R_24)))
  }
}
