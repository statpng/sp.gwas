#' @importFrom R.utils isDirectory

#' @export png.manhattan_from_dir
# png.manhattan_from_dir("./Kunyan/SNPs=2684 _RILs=157/lasso_perm_pop_mode3/", threshold = "Permuted", FD=c(1,5,10))
png.manhattan_from_dir <- function(path = "./sp_result", threshold = c("none", "Theoretical", "Permuted"), FD = NULL){
  # path : a path of directory containing the sp.gwas result
  # FD : a vector of the number of false discoveries to be controlled
  
  # path = "./Kunyan/SNPs=2684 _RILs=157/lasso_perm_pop_mode3/"
  # threshold <- "Theoretical"
  # FD = c(1,5,10)
  
  if( !isDirectory(path) ) stop("There doesn't exist the path provided")
  
  
  # Manhattan plot for the first phenotype (Y1)
  sp_results <- read.csv( paste0(path, "/[2]sp.results.csv") )
  sp_results$chr <- as.numeric( as.character(sp_results$chr) )
  sp_results$pos <- as.numeric( as.character(sp_results$pos) )
  
  
  if( file.exists( paste0(path, "/[2]sp.thresholds.csv") ) ){
    thresholds <- read.csv( paste0(path, "/[2]sp.thresholds.csv") )
  } else {
    thresholds <- read.csv( paste0(path, "/[2]sp.thresholds_before_perm.csv") )
  }
  
  colnames(thresholds)[2] <- "FalseDiscovery"
  q <- ncol(thresholds)-2
  
  threshold_final <- subset( thresholds, Method == threshold & FalseDiscovery %in% FD )
  
  if( nrow(threshold_final) == 0 ) {
    if( is.null(FD) ) stop("FD should be provided if you want to draw a cutoff line")
  }
  
  for( kk in 1:q ){
    df <- sp_results[,c(1:3, kk+3)]
    colnames(df)[ncol(df)] <- "Y"
    
    threshold_final_Y <- threshold_final[,c(1:2,kk+2)]
    colnames(threshold_final_Y)[3] <- "Y"
    
    for( fdfd in FD ){
      threshold_final_Y_FD <- subset(threshold_final_Y, FalseDiscovery == fdfd )
      rs.highlighted <- df$rs[df$Y >= threshold_final_Y_FD$Y]
      thr <- threshold_final_Y_FD$Y
      if( length(thr) == 0 ) thr <- FALSE
      
      jpeg(paste0(path, "/[3]manhattan_method=",threshold,"_Y=",kk,"_FD=",fdfd,".jpeg"), width=12, height=5, unit="in", res=600)
      qqman_manhattan(df,
                      chr="chr", bp="pos", snp="rs", p="Y", logp=FALSE,
                      col=c("green", "black", "red", "orange", "gold3", "deeppink", "green4", "lightblue4"),
                      suggestiveline = thr,
                      # highlight = list(highlight_theory, highlight_perm, highlight_intersect),
                      highlight = rs.highlighted,
                      # col.highlight = c("blue", "red", "purple"),
                      col.highlight = "blue",
                      ylab="Selection probabilities", ylim=c(0,1.01) )  
      dev.off()
    }
  }
  
  
}


