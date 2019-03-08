SetFolder <- function(FolderName){
    if( !file.exists( paste0("./",FolderName) )) system( paste0("mkdir ./", FolderName) )
    setwd(paste0("./", FolderName))
}
