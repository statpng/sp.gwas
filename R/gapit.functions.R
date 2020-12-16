`GAPIT.Numericalization` <-
  function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
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
      x[x=="K"]="Z" #K (for GT genotype)is is replaced by Z to ensure heterozygose has the largest value
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
    if(len<=1 | len> 3)x=0
    
    #2 status
    if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,1))
    
    #3 status
    if(bit==1){
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
    }else{
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],2,1)))
    }
    
    #print(paste(lev,len,sep=" "))
    #print(position)
    
    #missing data imputation
    if(impute=="Middle") {x[is.na(x)]=1 }
    
    if(len==3){
      if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
      if(impute=="Major")  {x[is.na(x)]=position[len]-1}
      
    }else{
      if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
      if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
    }
    
    #alternative genetic models
    if(effect=="Dom") x=ifelse(x==1,1,0)
    if(effect=="Left") x[x==1]=0
    if(effect=="Right") x[x==1]=2
    
    if(byRow) {
      result=matrix(x,n,1)
    }else{
      result=matrix(x,1,n)  
    }
    
    return(result)
  }#end of GAPIT.Numericalization function
#=============================================================================================


















# # old ---------------------------------------------------------------------
# 
# `GAPIT.Numericalization` <-
#   function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
#     #Object: To convert character SNP genotpe to numerical
#     #Output: Coresponding numerical value
#     #Authors: Feng Tian and Zhiwu Zhang
#     # Last update: May 30, 2011
#     ##############################################################################################
#     if(bit==1)  {
#       x[x=="X"]="N"
#       x[x=="-"]="N"
#       x[x=="+"]="N"
#       x[x=="/"]="N"
#       x[x=="K"]="Z" #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
#     }
#     
#     if(bit==2)  {
#       x[x=="XX"]="N"
#       x[x=="--"]="N"
#       x[x=="++"]="N"
#       x[x=="//"]="N"
#       x[x=="NN"]="N"
#       x[x=="00"]="N"
#       
#     }
#     
#     n=length(x)
#     lev=levels(as.factor(x))
#     lev=setdiff(lev,"N")
#     #print(lev)
#     len=length(lev)
#     #print(len)
#     #Jiabo creat this code to convert AT TT to 1 and 2. 2018.5.29
#     if(bit==2)
#     {
#       inter_store=c("AT","AG","AC","TA","GA","CA","GT","TG","GC","CG","CT","TC")
#       inter=intersect(lev,inter_store)
#       if(length(inter)>1)
#       {
#         x[x==inter[2]]=inter[1]
#         n=length(x)
#         lev=levels(as.factor(x))
#         lev=setdiff(lev,"N")
#         #print(lev)
#         len=length(lev)
#       }
#       if(len==2&bit==2)
#       { #inter=intersect(lev,inter_store)
#         if(!is.na(inter[1]))
#         {
#           lev=union(lev,"UU")
#           len=len+1
#           
#         }
#       }
#       if(len==3&bit==2)
#       {
#         inter=intersect(lev,inter_store)
#       }
#       
#     }
#     #print(lev)
#     #print(len)
#     #Jiabo code is end here
#     
#     #Genotype counts
#     count=1:len
#     for(i in 1:len){
#       count[i]=length(x[(x==lev[i])])
#     }
#     
#     #print(count)
#     
#     if(Major.allele.zero){
#       if(len>1 & len<=3){
#         #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
#         if(bit==1){
#           count.temp = cbind(count, seq(1:len))
#           if(len==3) count.temp = count.temp[-3,]
#           count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
#           if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
#         }
#         
#         #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
#         if(bit==2){
#           count.temp = cbind(count, seq(1:len))
#           if(len==3) count.temp = count.temp[-2,]
#           count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
#           if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
#         }
#         
#         count = count[order]
#         lev = lev[order]
#         
#       }   #End  if(len<=1 | len> 3)
#     } #End  if(Major.allele.zero)
#     
#     #print(x)
#     
#     #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
#     if(bit==1 & len==3){
#       temp=count[2]
#       count[2]=count[3]
#       count[3]=temp
#     }
#     
#     #print(lev)
#     #print(count)
#     position=order(count)
#     
#     #Jiabo creat this code to convert AT TT to 1 and 2.2018.5.29
#     
#     lev1=lev
#     if(bit==2&len==3)
#     {
#       lev1[1]=lev[count==sort(count)[1]]
#       lev1[2]=lev[count==sort(count)[2]]
#       lev1[3]=lev[count==sort(count)[3]]
#       position=c(1:3)
#       lev=lev1
#     }
#     #print(lev)
#     #print(position)
#     #print(inter)
#     #Jiabo code is end here
#     
#     
#     #1status other than 2 or 3
#     if(len<=1 | len> 3)x=0
#     
#     #2 status
#     if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,2))
#     
#     #3 status
#     if(bit==1){
#       if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
#     }else{
#       if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[lev!=inter][1],0,ifelse(x==inter,1,2)))
#     }
#     
#     #print(paste(lev,len,sep=" "))
#     #print(position)
#     
#     #missing data imputation
#     if(impute=="Middle") {x[is.na(x)]=1 }
#     
#     if(len==3){
#       if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
#       if(impute=="Major")  {x[is.na(x)]=position[len]-1}
#       
#     }else{
#       if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
#       if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
#     }
#     
#     #alternative genetic models
#     if(effect=="Dom") x=ifelse(x==1,1,0)
#     if(effect=="Left") x[x==1]=0
#     if(effect=="Right") x[x==1]=2
#     
#     if(byRow) {
#       result=matrix(x,n,1)
#     }else{
#       result=matrix(x,1,n)
#     }
#     
#     return(result)
#   }#end of GAPIT.Numericalization function
# #=============================================================================================









`GAPIT.kinship.VanRaden` <-
  function(snps,hasInbred=TRUE) {
    # Object: To calculate the kinship matrix using the method of VanRaden (2009, J. Dairy Sci. 91:4414???C4423)
    # Input: snps is n individual rows by m snps columns
    # Output: n by n relationship matrix
    # Authors: Zhwiu Zhang
    # Last update: March 2, 2016 
    ############################################################################################## 
    print("Calculating kinship with VanRaden method...")
    #Remove invariants
    fa=colSums(snps)/(2*nrow(snps))
    index.non=fa>=1| fa<=0
    snps=snps[,!index.non]
    
    nSNP=ncol(snps)
    nInd=nrow(snps)
    n=nInd 
    
    ##allele frequency of second allele
    p=colSums(snps)/(2*nInd)
    P=2*(p-.5) #Difference from .5, multiple by 2
    snps=snps-1 #Change from 0/1/2 coding to -1/0/1 coding
    
    print("substracting P...")
    Z=t(snps)-P#operation on matrix and vector goes in direction of column
    print("Getting X'X...")
    #K=tcrossprod((snps), (snps))
    K=crossprod((Z), (Z)) #Thanks to Peng Zheng, Meng Huang and Jiafa Chen for finding the problem
    
    print("Adjusting...")
    adj=2*sum(p*(1-p))
    K=K/adj
    
    print("Calculating kinship with VanRaden method: done")
    
    return(K)
  }
#=============================================================================================

`GAPIT.kinship.Zhang` <-
  function(snps,hasInbred=TRUE) {
    # Object: To calculate ZHANG (Zones Harbored Adjustments of Negligent Genetic) relationship
    # Authors: Zhwiu Zhang
    # Last update: october 25, 2014 
    ############################################################################################## 
    print("Calculating ZHANG relationship defined by Zhiwu Zhang...")
    #Remove invariants
    fa=colSums(snps)/(2*nrow(snps))
    index.non=fa>=1| fa<=0
    snps=snps[,!index.non]
    
    het=1-abs(snps-1)
    ind.sum=rowSums(het)
    fi=ind.sum/(2*ncol(snps))
    inbreeding=1-min(fi)
    
    nSNP=ncol(snps)
    nInd=nrow(snps)
    n=nInd 
    snpMean= apply(snps,2,mean)   #get mean for each snp
    print("substracting mean...")
    snps=t(snps)-snpMean    #operation on matrix and vector goes in direction of column
    print("Getting X'X...")
    #K=tcrossprod((snps), (snps))
    K=crossprod((snps), (snps)) 
    if(is.na(K[1,1])) stop ("GAPIT says: Missing data is not allowed for numerical genotype data")
    
    print("Adjusting...")
    #Extract diagonals
    i =1:n
    j=(i-1)*n
    index=i+j
    d=K[index]
    DL=min(d)
    DU=max(d)
    floor=min(K)
    
    
    #Set range between 0 and 2
    top=1+inbreeding
    K=top*(K-floor)/(DU-floor)
    Dmin=top*(DL-floor)/(DU-floor)
    
    #Adjust based on expected minimum diagonal (1)
    if(Dmin<1) {
      print("Adjustment by the minimum diagonal")
      K[index]=(K[index]-Dmin+1)/((top+1-Dmin)*.5)
      K[-index]=K[-index]*(1/Dmin)
    }
    
    #Limiting the maximum offdiagonal to the top
    Omax=max(K[-index])
    if(Omax>top){
      print("Adjustment by the minimum off diagonal")
      K[-index]=K[-index]*(top/Omax)
    }
    
    print("Calculating kinship with Zhang method: done")
    return(K)
  }
#=============================================================================================

`GAPIT.kinship.loiselle` <-
  function(snps, method="additive", use="all") {
    # Object: To calculate the kinship matrix using the method of Loiselle et al. (1995)
    # Authors: Alex Lipka and Hyun Min Kang
    # Last update: May 31, 2011 
    ############################################################################################## 
    #Number of SNP types that are 0s
    n0 <- sum(snps==0,na.rm=TRUE)
    #Number of heterozygote SNP types
    nh <- sum(snps==0.5,na.rm=TRUE)
    #Number of SNP types that are 1s
    n1 <- sum(snps==1,na.rm=TRUE)
    #Number of SNP types that are missing
    nNA <- sum(is.na(snps))
    
    
    
    #Self explanatory
    dim(snps)[1]*dim(snps)[2]
    #stopifnot(n0+nh+n1+nNA == length(snps))
    
    
    #Note that the two lines in if(method == "dominant") and if(method == "recessive") are found in
    #if(method == "additive").  Worry about this only if you have heterozygotes, which you do not.
    if( method == "dominant" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
    }
    else if( method == "recessive" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
    }
    else if( ( method == "additive" ) && ( nh > 0 ) ) {
      dsnps <- snps
      rsnps <- snps
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      dsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      rsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
      snps <- rbind(dsnps,rsnps)
    }
    
    #mafs is a (# SNPs)x(# lines) matrix.  The columns of mafs are identical, and the ij^th element is the average
    #allele frequency for the SNP in the i^th row.
    
    #if(use == "all") imputes missing SNP type values with the expected (average) allele frequency.
    if( use == "all" ) {
      mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
      snps[is.na(snps)] <- mafs[is.na(snps)]
    }
    else if( use == "complete.obs" ) {
      mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
      snps <- snps[rowSums(is.na(snps))==0,]
    }
    mafs_comp <- 1-mafs
    snps_comp <- 1-snps
    
    
    n <- ncol(snps)
    K <- matrix(nrow=n,ncol=n)
    diag(K) <- 1
    #Create the k term on page 1422 of Loiselle et al. (1995)
    
    missing <- rep(NA, dim(snps)[1])  
    for(i in 1:dim(snps)[1]) {
      missing[i] <- sum(is.na(snps[i,]))
    }
    
    
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        Num_First_Term_1 <- (snps[,i]-mafs[,i])*(snps[,j]-mafs[,j])
        Num_First_Term_2 <- (snps_comp[,i]-mafs_comp[,i])*(snps_comp[,j]-mafs_comp[,j])
        First_Term <- sum(Num_First_Term_1)+sum(Num_First_Term_2)
        
        Num_Second_Term_1 <- mafs[,i]*(1-mafs[,i])
        Num_Second_Term_2 <- mafs_comp[,i]*(1-mafs_comp[,i])
        Num_Second_Term_Bias_Correction <- 1/((2*n)-missing - 1)
        Num_Second_Term <-  Num_Second_Term_1 + Num_Second_Term_2
        Second_Term <- sum(Num_Second_Term*Num_Second_Term_Bias_Correction)
        
        Third_Term <- sum(Num_Second_Term) 
        
        f <- (First_Term + Second_Term)/Third_Term
        
        K[i,j] <- f
        if(K[i,j]<0) K[i,j]=0
        
        K[j,i] <- K[i,j]
      }
    }
    return(K)
  }
#=============================================================================================

`GAPIT.kinship.separation` <-
  function(PCs=NULL,EV=NULL,nPCs=0 ){
    #Object: To calculate kinship from PCS
    #       PCs: the principal component as columns and individual as rows, the first column is taxa
    #       EV: Eigen values
    #       nPCs: the number of front PCs excluded to calculate kinship
    #Output: kinship
    #Authors: Huihui Li and Zhiwu Zhang
    #Last update: April 17, 2012
    ##############################################################################################
    print("Calling GAPIT.kinship.separation")  
    Total.number.PCs=ncol(PCs)
    n=nrow(PCs)
    print(Total.number.PCs)
    print(n)
    #Choose Total.number.PCs-nPCs PCs and EV to calculate K
    sep.PCs=PCs[, (nPCs+2):(Total.number.PCs)]  #first column is taxa
    sep.EV=EV[(nPCs+1):Total.number.PCs]
    
    Weighted.sep.EV=sep.EV/sum(sep.EV)
    
    #X=t(t(sep.PCs)*Weighted.sep.EV)  
    X=sep.PCs
    
    XMean= apply(X,2,mean)
    X=as.matrix(X-XMean)
    K=tcrossprod((X), (X))
    
    #Extract diagonals
    i =1:n
    j=(i-1)*n
    index=i+j
    d=K[index]
    DL=min(d)
    DU=max(d)
    floor=min(K)
    
    K=(K-floor)/(DL-floor)
    MD=(DU-floor)/(DL-floor)
    
    if(is.na(K[1,1])) stop ("GAPIT says: Missing data is not allowed for numerical genotype data")
    if(MD>2)K[index]=K[index]/(MD-1)+1
    print("GAPIT.kinship.separation called succesfuly")
    return (K)
  }
#=============================================================================================
