#' @title 
#' Core algorithm
#' 
#' @description 
#' It's a subfunction of CIBERSORT. 
#' SVR algorithm, the core algorithm of CIBERSORT, referenced by CIBERSORT to compute infiltration of immune cells and it can't run alone.
#'
#' @details 
#' Download this file from https://cibersort.stanford.edu/
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' There is no change compared to the source code from CIBERSORT.
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' 
#' @param X 
#' cell-specific gene expression.
#' 
#' @param y 
#' mixed expression per sample.
#' 

CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}
