#' 
#' @title Do permutations
#' 
#' @description 
#' It's a subfunction of CIBERSORT. 
#' To do permutations for deconvolution algorithm. Tt is the key to produce p value.
#' It referenced by CIBERSORT to compute infiltration of immune cells and it can't run alone.
#' 
#' @details 
#' Download this file from https://cibersort.stanford.edu/
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' There is no change compared to the source code from CIBERSORT.
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#'
#' @param perm 
#' Number of permutations
#' 
#' @param X 
#' Cell-specific gene expression
#' 
#' @param y 
#' Mixed expression per sample

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}
