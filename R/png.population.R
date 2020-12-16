#' @import ggplot2
#' @importFrom  gtools mixedsort
#' @import tsne
#' @import cluster
#' @importFrom stats kmeans
#' @export png.population
png.population <- function(snp, 
                           cluster.method = c("kmeans", "pam"), 
                           gap.method = c("firstSEmax", "globalSEmax", "globalmax", "Tibs2001SEmax"),
                           pca.scale = FALSE,
                           K.max = NULL,
                           nboot = 100,
                           tsne = FALSE,
                           save.path = "./"
){
  
  if( missing(cluster.method) ) cluster.method <- "kmeans"
  if( missing(gap.method) ) gap.method <- "firstSEmax"
  
  # colors <-
  #   c(`1`= "red",
  #     `2`= "green",
  #     `3`= "blue",
  #     `4`= "gold1",
  #     `5`= "pink",
  #     `6`= "skyblue",
  #     `7`= "orange",
  #     `0`= "black")
  
  colors <- c("red",
              "green",
              "blue",
              "gold1",
              "pink",
              "skyblue",
              "orange",
              "black",
              "antiquewhite4",
              "aquamarine4",
              "azure4",
              "bisque3",
              "black",
              "blue3",
              "blueviolet",
              "brown3",
              "burlywood4",
              "cadetblue4",
              "chartreuse4",
              "chocolate4",
              "cornflowerblue",
              "cyan4",
              "darkgoldenrod4",
              "darkgreen",
              "darkmagenta",
              "darkolivegreen4",
              "darkorange3",
              "darkorchid4",
              "darkred",
              "darksalmon",
              "darkslateblue",
              "darkseagreen4",
              "deepskyblue4")
  
  
  
  
  
  
  
  
  
  # PCA
  pca.fit <- prcomp(snp, scale=pca.scale)
  
  # Clustering --------------------------------------------------------------
  # We do not recommend using hierarchical clustering here
  # hclusCut <- function(x, k, d.meth = "euclidean", ...){
  #   list(cluster = cutree(hclust(dist(x, method=d.meth), ...), k=k))
  # }
  if( is.null(K.max) ){
    K.max <- floor( nrow(snp) * 0.1 )
  }
  
  if( cluster.method == "pam" ){
    pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
    fit.gap <- cluster::clusGap(snp, FUN=pam1, K.max=K.max, B=nboot)
  } else if(cluster.method == "kmeans"){
    fit.gap <- cluster::clusGap(snp, FUN=kmeans, K.max=K.max, B=nboot, iter.max = 100)
  }
  
  pdf(file=paste0(save.path, "/[3]GapPlot.pdf"), width=7, height=4)
  plot(fit.gap, xlab='The number of clusters', main=NA)
  dev.off()
  
  res.gap <- fit.gap$Tab
  nc <- maxSE(f = res.gap[,"gap"], 
              SE.f = res.gap[,"SE.sim"],
              method=gap.method, SE.factor=1)
  
  res.kmeans <- kmeans(pca.fit$x, centers = nc, iter.max = 100)
  
  group <- res.kmeans$cluster %>% as.character %>% as.factor %>% {factor(., levels = mixedsort(levels(.)) )}
  
  
  
  # PC plot ------------------------------------------------------------------
  
  ggplot( data=data.frame(pca.fit$x[,1:2], 
                          group=group ) ) +
    geom_point(aes(PC1, PC2, color=group), shape=18, size=3) +
    scale_color_manual(values=colors[unique(group)]) +
    xlab(paste0("PC1 (", round( {pca.fit$sdev[1]^2/sum(pca.fit$sdev^2)}, 3)*100 , "%)")) +
    ylab(paste0("PC2 (", round( {pca.fit$sdev[2]^2/sum(pca.fit$sdev^2)}, 3)*100 , "%)")) +
    theme_bw(base_size = 15) +
    # ggsave(filename = "./[Figure] PCplot_without_prop.pdf", width=8, height=5)
    ggsave(filename = paste0(save.path, "/[3]PCoA.pdf"), width=8, height=5)
  
  
  if( tsne ){
    tsne.list <- as.list(1:4)
    names(tsne.list) <- paste0("tsne", c(5, 10, 30, 50))
    
    count <- 0
    for( tsne.perplexity in c(5, 10, 30, 50) ){
      count <- count + 1
      # absKinship <- abs(Kinship)
      # Dissimilarity_from_Kinship <- diag(diag(sqrt(absKinship))) %*% absKinship %*% diag( diag(sqrt(absKinship)) )
      # res.tsne_kinship <- tsne::tsne( Dissimilarity_from_Kinship, perplexity=tsne.perplexity )
      # colnames(res.tsne_kinship) <- c("tSNE_1", "tSNE_2")
      # 
      # ggplot( data=data.frame(res.tsne_kinship, 
      #                         group=as.character(group) ) ) +
      #   geom_point(aes(tSNE_1, tSNE_2, color=group), shape=18, size=3) +
      #   scale_color_manual(values=colors[group],
      #                      labels=c(1:6)) +
      #   theme_bw(base_size = 15) +
      #   ggsave(filename = paste0(save.path, "[3]tSNE_Kinship_perplexity=", tsne.perplexity,".pdf"), width=7, height=5)
      
      
      res.tsne_euclidean <- tsne::tsne( dist(snp, "euclidean"), perplexity=tsne.perplexity )
      colnames(res.tsne_euclidean) <- c("tSNE_1", "tSNE_2")
      
      tsne.list[[count]] <- res.tsne_euclidean
      
      ggplot( data=data.frame(res.tsne_euclidean, 
                              group=group ) ) +
        geom_point(aes(tSNE_1, tSNE_2, color=group), shape=18, size=3) +
        scale_color_manual(values=colors[unique(group)]) +
        theme_bw(base_size = 15) +
        ggsave(filename = paste0(save.path, "/[3]tSNE_Euclidean_perplexity=", tsne.perplexity,".pdf"), width=7, height=5)
      
    }
  } else {
    tsne.list <- list( FALSE )
  }
  
  
  list( #Kinship = Kinship, 
        gap = fit.gap,
        group = group, 
        pc = pca.fit$x,
        tsne = tsne.list )
  
}




png.Kinship <- function(snp, method = c("VanRaden", "Zhang", "loiselle", "separation")){
  # snp : n individual rows and p snps columns
  
  if( kinship.method %in% c("VanRaden", "Zhang", "loiselle") ){
    Kinship <- switch( kinship.method,
                       "VanRaden" = GAPIT.kinship.VanRaden(snp),
                       "Zhang" = GAPIT.kinship.Zhang(snp),
                       "loiselle" = GAPIT.kinship.loiselle(t(snp))
    )
  }
  
  # PCA
  pca.fit <- prcomp(snp, scale=TRUE)
  
  if( kinship.method == "separation" ){
    Kinship <- GAPIT.kinship.separation(pca.fit$x)
  }
  
  return( Kinship )
}