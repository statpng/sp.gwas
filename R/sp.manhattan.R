#' @import CMplot
#' @export sp.glmnet
sp.manhattan <- function(sp.df,
                         threshold,
                         save.path,
                         manhattan.type = "m",
                         plot.ylim = NULL,
                         plot.name = "",
                         plot.type = c("pdf", "jpg"),
                         dpi = 300){

    
    devAskNewPage(ask = FALSE)
    par(ask=F)
    
    PLOT.NAME <- paste0(save.path,"/[3]CircularManhattanPlot_",plot.name)
        if (plot.type == "jpg")
            jpeg(paste(PLOT.NAME, ".jpg", sep = ""),
                 width = 5.5 * 2 * dpi, height = 5.5 * 2 * dpi,
                 res = dpi, quality = 100)
        if (plot.type == "pdf")
            pdf(paste(PLOT.NAME, ".pdf", sep = ""),
                width = 5.5 * 2, height = 5.5 * 2)
        if (plot.type == "tiff")
            tiff(paste(PLOT.NAME, ".tiff",
                       sep = ""), width = 5.5 * 2 * dpi, height = 5.5 * 2 *
                     dpi, res = dpi)

    multracks <- TRUE
    file.output <- FALSE
    
    if( manhattan.type == "r" ){
        manhattan.type <- "m"
        multracks <- FALSE
        file.output <- TRUE
        gd <- getwd()
        setwd(save.path)
    }
    
    sp.df %>%
        CMplot(plot.type=manhattan.type,
               col=matrix(c("orange", "grey30",
                            "darkblue", "green",
                            "darkgreen", "darkmagenta",
                            "gold", "black"), nrow=4, byrow=T)[1:(ncol(.)-3),],
               chr.labels=paste("Chr",c(1:20),sep=""),
               r=0.4,
               outward=FALSE,
               cir.chr.h=0.5,
               cir.legend=TRUE,
               cir.legend.cex=0.5,
               cir.legend.col="black",
               threshold = threshold[nrow(threshold), ],
               threshold.col = 1:ncol(threshold),
               threshold.lwd = 2,
               ylim=plot.ylim,
               LOG10 = FALSE,
               ylab = "Selection Probabilities",
               amplify = TRUE,
               multracks = multracks,
               signal.col = "red",
               signal.cex=1.5,
               chr.den.col="black",
               file.output=file.output,
               file="jpg", dpi=dpi)
    dev.off()

    
    if( manhattan.type == "r" ){
        setwd("..")
    }
    
    
}

