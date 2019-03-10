## sp.gwas
R package for an analysis of high-dimensional Single-Nucleotide-Polymorphism(SNP) data in the hapmap format.


## Installation
There is an install_github function to install R packages hosted on GitHub in the devtools package.
```{R}
install.packages("devtools")
libary(devtools)
```

Now then you can install **sp.gwas** package as follows.
```{R}
install_github("statpng/sp.gwas")
```

You may need to update your installed package which have dependencies with **sp.gwas**.
```{R}
update.packages()
```

If you cannot update some packages, try install **sp.gwas** after removing the packages that are not updated correctly from directory of R library, such as
```{R}
remove.packages("pkgname")
```
or removing package folder in R library("C:\Program Files\R\R-3.5.2\library").
