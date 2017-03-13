###############
## Functions to extend iBBiG
## EnsEMBL clustering
## 2017
## Aedin
##############

clusterWeights <-function(x) {
  log(apply(x@NumberxCol,1,sum)/ncol(x@NumberxCol)*(x@Clusterscores))
}

filter_iBBiG<-function(x,filterWeight=0, minClustSize=25) {
  require(iBBiG)
  which(clusterWeights(x) > filterWeight & rowSums(NumberxCol(x)*1) >= minClustSize)
}



ensemblClusters<-function(xx, savetmpfiles=FALSE, outpath="\tmp", outfile="ibbig_ensembl_merge", n=paste0("iter",1:length(xx)), nModules=10, alpha=0.3, ...) {
    # gsMat is  a list of iBBiG objects merge them
    # Need to define an iBBiG list S4 object, and method for it
    # in the meantime... will use this


    # EXAMPLE
    # binMat<-makeArtificial()
    # xx<-lapply(seq_along(1:10), function(x) iBBiG:::iBBiG(binMat@Seeddata, 10))
    # xx2<-ensemblClusters (xx)

    require(iBBiG)

    if(!is.list(xx)) stop("Expecting list")
    if (is.null(n)) n =names(xx)
    if (is.null(n)) n = paste0("iter", 1:length(xx))

    # Extract Seedata
    if(length(unique(sapply(xx, function(x) sum(x@Seeddata))))==1 | !identical(xx[[1]]@Seeddata, xx[[2]]@Seeddata)) print("list should contain multiple runs of iBBiG on the same data")

    #Extract RowScorexNumber from each iteration
    GS<- lapply(xx, function(x) x@RowScorexNumber)
    names(GS)<- n
    if(savetmpfiles) saveRDS(GS, file=file.path(outpath, paste0(outfile, "_GS.rds")))
    gsMat<-do.call(cbind, GS)
    nGS<-sapply(GS, ncol)
    colnames(gsMat)<-paste0(rep(names(nGS), nGS), "_", colnames(gsMat))

    #Extract NumberxCol from each iteration
    TT<-lapply(xx, function(x) t(x@NumberxCol))
    names(TT)<- n
    if(savetmpfiles) saveRDS(TT, file=file.path(outpath, paste0(outfile, "_Tumors.rds")))
    ttMat<-do.call(cbind, TT)
    nTT<-sapply(TT, ncol)
    colnames(ttMat)<-paste0(rep(names(nTT), nTT), "_", colnames(ttMat))

    # Check Consistency in Cluster (colnames) Names
    if(!identical(colnames(gsMat),colnames(ttMat))) print("RowScorexNumber and NumberxCol colnames differ")
    cluster= colnames(gsMat)

    # Check Consistency in Number of Clusters
    nc=unique(sapply(GS, ncol))
    if(!nc==unique(sapply(TT, ncol))) print("Col numbers of Tumors and GS differ")
    if(length(nc)>1) print("ncol in list of RDS files not differ")


    ##############################################################
    # Run Biclustering on all NumberxCol binary matrices,
    # to find consistent clusters
    ####################################################
     clus<-iBBiG(ttMat,nModules=nModules,...)

    # Check Consistency in Cluster (colnames) Names
    identical(colnames(clus@NumberxCol),colnames(gsMat))

    # Create a New NumberxCol Matrix
    nc= (clus@NumberxCol*1)%*%t(gsMat)

    # MAKE NEW IBBIG OBJECT OF THE ENSEMBL
    x<- new("iBBiG", Seeddata = xx[[1]]@Seeddata, RowxNumber = clus@RowxNumber, RowScorexNumber = clus@RowScorexNumber, NumberxCol = (nc>0) ,Number = clus@Number, Clusterscores=clus@Clusterscores)


    # Need to create some stats on freq of tumor in final clusters etc
    # For moment.. data dump
    x@info=list(IterRowScorexNumber= gsMat, IterNumberxCol=ttMat, NumberxScore=nc)
    if(!validObject(x)) print("Can't create valid iBBiG Object")

    return(x)
}



mergeClusters<- function(x, i, drop=FALSE) {
  # x is an ibbig object
  ind<-1:x@Number
  i= paste("M", i)
  i=which(names(x@Clusterscores)%in%i)
  nlab= paste0("M ",paste0(i, collapse = "_"))
  if (!all(i%in%ind)) stop(paste("Invalid select"))
  x@info<-list("originalOrder"=ind, "merged"=i)

  rn= rowMeans(x@RowxNumber[,i, drop=drop])
  x@RowxNumber<-x@RowxNumber[,-c(i) ,drop=drop]
  x@RowxNumber<-cbind(x@RowxNumber, rn)
  colnames(x@RowxNumber)[colnames(x@RowxNumber)=='rn'] = nlab


  rn= rowMeans(x@RowScorexNumber[,i, drop=drop])
  x@RowScorexNumber<-x@RowScorexNumber[,-c(i) ,drop=drop]
  x@RowScorexNumber<-cbind(x@RowScorexNumber, rn)
  colnames(x@RowScorexNumber)[colnames(x@RowScorexNumber)=='rn'] = nlab

  rc=colMeans(x@NumberxCol[c(i), ,drop=drop])
  x@NumberxCol<-x@NumberxCol[-c(i), , drop=drop]
  x@NumberxCol<-rbind(x@NumberxCol, rc)
  rownames(x@NumberxCol)[rownames(x@NumberxCol)=='rc'] = nlab

  if (!is.null(x@info$NumberxScore)){
    rc=colMeans(x@info$NumberxScore[c(i), ,drop=drop])
    x@info$NumberxScore<-x@info$NumberxScore[-c(i), , drop=drop]
    x@info$NumberxScore<-rbind(x@info$NumberxScore, rc)
    rownames(x@info$NumberxScore)[rownames(x@info$NumberxScore)=='rc'] = nlab
  }

}