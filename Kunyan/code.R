pkgs <- "ggplot2, compiler, glmnet, dplyr, CMplot, readxl, data.table, writexl, bestNormalize, gridExtra, tidyr, gtools, calibrate, genetics, pbapply, HardyWeinberg, tsne, cluster, stats, R.utils, RcppArmadillo, RcppProgress"

for( i in strsplit(pkgs, ", ")[[1]] ){
  library(i, character.only=TRUE)
}

install.packages("./SNPs=2684_RILs=157/sp.gwas_1.4.1.tar.gz", repos = NULL, type = "source")
# install.packages("./SNPs=2684_RILs=157/sp.gwas_1.4.1.zip", repos = NULL, type = "win.binary")


library(sp.gwas)

genotype <- read.csv("./SNPs=2684_RILs=157/RIL genotype.csv", stringsAsFactors = FALSE)
phenotype <- read.csv("./SNPs=2684_RILs=157/RIL phenotype.csv", stringsAsFactors = FALSE)

colnames(genotype)[3] <- "chr"
genotype <- rbind( colnames(genotype), genotype )

genotype$chr <- gsub( pattern = "chr([0-9])", "\\1", genotype$chr )

genotype[1:10,1:11]
str(phenotype)


png.venn2 = function(x,y,duplicated=FALSE){
  inner = intersect(x,y)
  xnyc = x[!x%in%inner]
  ynxc = y[!y%in%inner]
  if(!duplicated){
    inner = unique(inner)
    xnyc = unique(xnyc)
    ynxc = unique(ynxc)
  }
  list( x = xnyc , y = ynxc , inner = inner )
}

png.venn2(phenotype$Taxa, unlist(genotype[1,-(1:11)]))





# The lasso ---------------------------------------------------------------

sp.gwas(genotype = genotype, 
        phenotype = phenotype, 
        input.type = c("object", "path")[1], 
        QC = TRUE,
        imputation = TRUE,
        impute.type = "mode",
        population = TRUE, 
        remove.missingY = TRUE,
        save.path = "./SNPs=2684_RILs=157/lasso",
        y.id.col = 1, 
        y.col = 2:3, 
        normalization = FALSE, 
        method="lasso",
        family="multinomial",
        false.discovery = c(1:20),
        permutation = TRUE,
        plot.ylim = NULL,
        lambda.min.quantile = 0.5,
        alpha.seq = 1,
        n.lambda = 20,
        K = 100,
        psub = 0.8,
        manhattan.type = c("c", "r")[1],
        plot.name = "Test",
        plot.type = "jpg",
        plot.dpi = 300)





# Elastic-net -------------------------------------------------------------

sp.gwas(genotype = genotype, 
        phenotype = phenotype, 
        input.type = c("object", "path")[1], 
        QC = TRUE,
        imputation = TRUE,
        remove.missingY = TRUE,
        save.path = "./SNPs=2684_RILs=157/enet",
        y.id.col = 1, 
        y.col = 2:3, 
        normalization = TRUE, 
        method="enet",
        family="multinomial",
        false.discovery = c(1:20),
        permutation = TRUE,
        plot.ylim = NULL,
        lambda.min.quantile = 0.5,
        alpha.seq = seq(0.5,0.9,by=0.2),
        n.lambda = 20,
        K = 100,
        psub = 0.8,
        manhattan.type = c("c", "r")[1],
        plot.name = "Test",
        plot.type = "jpg",
        plot.dpi = 300)



sp.gwas::png.manhattan_from_dir("./SNPs=2684_RILs=157/lasso/lasso_seed=20201110__2020-11-10_153338/",threshold = "Theoretical", FD = c(1, 5, 10) )
sp.gwas::png.manhattan_from_dir("./SNPs=2684_RILs=157/lasso/lasso_seed=20201110__2020-11-10_153338/",threshold = "Permuted", FD = c(1, 5, 10) )

sp.gwas::png.manhattan_from_dir("./SNPs=2684_RILs=157/enet/enet_seed=20201110__2020-11-10_155308/",threshold = "Theoretical", FD = c(1, 5, 10) )
sp.gwas::png.manhattan_from_dir("./SNPs=2684_RILs=157/enet/enet_seed=20201110__2020-11-10_155308/",threshold = "Permuted", FD = c(1, 5, 10) )

