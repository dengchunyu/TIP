#' The functions of compute infiltration.
#' 
#' There is a little difference from the source code from CIBERSORT,but no change in nature.
#' 
#' @param sig_matrix signature matrix data import
#' @param mixture_file heterogenous mixed expression,
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @param saveDir the storage path of results
#' @param sigMat signature matrix name for choose."LM22" or "LM14"
#' @param sample "multiple" or "single"
#' @export
#' 

CIBERSORT <- function(sig_matrix, mixture_file, perm=100, QN=TRUE,saveDir="",sigMat="LM22",sample="multiple"){
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  X <- sig_matrix
  Y <- mixture_file
  X <- data.matrix(X)
  X <- X[order(rownames(X)),]
  
 # Two processing path for parameter 'sample'
  if(sample=="multiple"){
    Y <- data.matrix(Y)
    Y <- Y[order(rownames(Y)),]
    
  }else if(sample=="single"){
    Y <- as.matrix(Y)
    rownames_y<-rownames(Y)[order(rownames(Y))]
    colnames_y<-colnames(Y)
    Y <- matrix(Y[order(rownames(Y)),],ncol = 1)
    rownames(Y)<-rownames_y
    colnames(Y)<-colnames_y
  }
  
  P <- perm #number of permutations
  
  
  if(max(Y) < 50) {Y <- 2^Y} 
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  if(sample=="multiple"){
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
  }else if(sample=="single"){
    
    Y <- matrix(Y[YintX,],ncol=1)
    rownames(Y)<- Ygns[YintX]
    colnames(Y)<- colnames_y
    XintY <- Xgns %in% row.names(Y)
    
  }
  
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE") 
  
  
  
  #print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    w <- round(result$w,3)
    mix_r <- round(result$mix_r,3)
    mix_rmse <- round(result$mix_rmse,3)
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,round(pval,5),mix_r,mix_rmse)
    
    if(itor == 1) {
      output <- out
    }else {
      output <- rbind(output, out)
    }
    
    itor <- itor + 1
    
  }
  
  
  
  #save results
  #Two processing path for the choice of parameter 'sigMat'
  if(sigMat=="LM14"){
    if(sample=="single"){
      #saveDir="C:/" filename="deng"
      filename<-paste0(saveDir,"/CIBER_lm14_User_SingleSample_Result.txt")
    }else if(sample=="multiple"){
      filename<-paste0(saveDir,"/CIBER_User_Allsample_Result.txt")
    }
    
    write.table(rbind(header,output), file=filename, sep="\t", row.names=F, col.names=F, quote=F)
    
  }else if(sigMat=="LM22"){
    
    if(sample=="single"){
      filename<-paste0(saveDir,"/CIBER_lm22_User_SingleSample_Result.txt")
    }else if(sample=="multiple"){
      filename<-paste0(saveDir,"/CIBER_User_Allsample_Result.txt")
    }
    
    write.table(rbind(header,output), file=filename, sep="\t", row.names=F, col.names=F, quote=F)
    
  }
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  if(sample=="single"){
    obj<-matrix(obj,nrow = 1)
  }
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}