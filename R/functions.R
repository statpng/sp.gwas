SetFolder <- function(FolderName){
    if( !file.exists( paste0("./",FolderName) )) system( paste0("mkdir ./", FolderName) )
    setwd(paste0("./", FolderName))
}



mat.order <- function(mat, row=TRUE, col=TRUE){
  if( row ){
    row.order <- mixedorder(rownames(mat))
  } else {
    row.order <- TRUE
  }
  
  if( col ){
    col.order <- mixedorder(colnames(mat))
  } else {
    col.order <- TRUE
  }
  return( mat[row.order, col.order] )
}