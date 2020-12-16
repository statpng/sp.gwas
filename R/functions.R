SetFolder <- function(FolderName){
    if( !file.exists( paste0("./",FolderName) )) system( paste0("mkdir ./", FolderName) )
    setwd(paste0("./", FolderName))
}


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




png.save.figure <- function(expr, file="./Figure/figure.jpeg", print=T, save=T, ...){
  
  extension <- strsplit( file, "\\." )[[1]] %>% .[length(.)]
  if(save){
    switch(extension,
           png = {
             png(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           jpeg = {
             jpeg(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           pdf = {
             pdf(file, ...)
             eval(parse(text=expr))
             dev.off()
           },
           eps = {
             setEPS()
             postscript(file, ...)
             # postscript("./Figure/Manhattan.eps", width=10, height=5)
             eval(parse(text=expr))
             dev.off()
           },
           
           {stop("filetype not recognized")}
    )
  }
  
  #print?
  if (print) eval(parse(text=expr))
  invisible(NULL)
  
}	




#' @export png.Numericalization
png.Numericalization <- function(x,bit=2,Major.allele.zero = FALSE, byRow=TRUE){
  #Object: To convert character SNP genotpe to numerical
  #Output: Coresponding numerical value
  #Authors: Feng Tian and Zhiwu Zhang
  # Last update: May 30, 2011 
  ##############################################################################################
  if(bit==1)  {
    x[x=="X"]="N"
    x[x=="-"]="N"
    x[x=="+"]="N"
    x[x=="/"]="N"
  }
  
  if(bit==2)  {
    x[x=="XX"]="N"
    x[x=="--"]="N"
    x[x=="++"]="N"
    x[x=="//"]="N"
    x[x=="NN"]="N"
  }
  
  n=length(x)
  lev=levels(as.factor(x))
  lev=setdiff(lev,"N")
  #print(lev)
  len=length(lev)
  #print(lev)
  
  
  
  #Genotype counts
  count=1:len
  for(i in 1:len){
    count[i]=length(x[(x==lev[i])])
  }
  
  
  
  if(Major.allele.zero){
    if(len>1 & len<=3){
      #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
      if(bit==1){ 
        count.temp = cbind(count, seq(1:len))
        if(len==3) count.temp = count.temp[-3,]
        count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
        if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
      }
      
      #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
      if(bit==2){
        count.temp = cbind(count, seq(1:len))
        if(len==3) count.temp = count.temp[-2,]
        count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
        if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
      }
      
      count = count[order]
      lev = lev[order]
      
    }   #End  if(len<=1 | len> 3)
  } #End  if(Major.allele.zero)
  
  
  
  #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
  if(bit==1 & len==3){
    temp=count[2]
    count[2]=count[3]
    count[3]=temp
  }
  position=order(count)
  
  
  #1status other than 2 or 3
  if(len<=1 | len> 3) x=0
  
  #2 status
  if(len==2) x=ifelse(x=="N",NA,ifelse(x==lev[1],0,1))
  
  #3 status
  if(bit==1){
    if(len==3) x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
  }else{
    if(len==3) x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],2,1)))
  }
  
  
  if(byRow) {
    result=matrix(x,n,1)
  }else{
    result=matrix(x,1,n)  
  }
  
  return(result)
}
