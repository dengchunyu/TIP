#' 
#' @title Do permutations
#' 
#' @description 
#' It's a subfunction of CIBERSORT. 
#' To do permutations for deconvolution algorithm. Tt is the key to produce p value.
#' It referenced by CIBERSORT to compute infiltration of immune cells and can't run alone.
#' 
#'
#' @param perm 
#' Number of permutations
#' 
#' @param X 
#' cell-specific gene expression
#' 
#' @param y 
#' mixed expression per sample

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
