#' @@ superpositionRank
#' @description Get rank of a matrix based on superposition by gradient. 
#' Required for getting rank of elements in every column of a matrix.
#' 
#' @param  raw.matrix A matrix consist of numeric values. 
#' @param  m Numeric value indicating number of row of raw.matrix.   
#' @param  n Numeric value indicating number of column of raw.matrix.   
#' @returnType matrix
#' @return rank of all values in the matrix

superpositionRank <- function(raw.matrix, m, n){
  max.element <- max(raw.matrix)
  over.matrix <- matrix(rep(seq(0, max.element*(n-1), by = max.element), each=m), nrow=m)
  superposition <- raw.matrix + over.matrix
  
  superposition.rank <- matrix(rank(superposition, ties.method = "first"), nrow = m)
  
  return(superposition.rank)
}



#' @@ rankMatrixByCol
#' @description With the rank of all values in the matrix(superposition.rank), restore the rank of elements in every column.
#' 
#' @param  superposition.rank A matrix indicating rank of all values in the raw matrix.
#' @param  m Numeric value indicating number of row of matrix 'superposition.rank'.  
#' @param  n Numeric value indicating number of column of matrix 'superposition.rank'.
#' @returnType matrix
#' @return rank of elements in every column

rankMatrixByCol <- function(superposition.rank, m, n){
  minus.matrix <- matrix(rep(seq(0, m*(n-1), by=m), each=m), nrow=m)
  subtract.matrix <- superposition.rank - minus.matrix
  
  return(subtract.matrix)
}



#' @@ sortMatrixByCol
#' @description Sort every column of a matrix. 
#' 
#' @param  raw.matrix A matrix consist of numeric value. 
#' @param  set.decreasing Logical, if true sort is decreasing, default by TRUE.

sortMatrixByCol <- function(raw.matrix, set.decreasing=TRUE){
  if(is.null(dim(raw.matrix))) {
    subtract.matrix <- matrix(raw.matrix, nrow = 1)
    
  }else{
    n.num <- ncol(raw.matrix)
    m.num <- nrow(raw.matrix)
    #superposition
    over.matrix <- matrix(rep(seq(0, max(raw.matrix)*(n.num-1), by = max(raw.matrix)), each=m.num), nrow=m.num)
    
    if(set.decreasing){
      superposition <- raw.matrix - over.matrix
      superposition.sort <- sort(superposition, decreasing=set.decreasing)
      subtract.matrix <- superposition.sort + over.matrix
    }else{ 
      superposition <- raw.matrix + over.matrix
      superposition.sort <- sort(superposition, decreasing=set.decreasing)
      subtract.matrix <- superposition.sort - over.matrix
    }
  }
  return(subtract.matrix)
}



#' @@ indicatorInsideGeneSetMatrix
#' @description With the rank of all values in the matrix(superposition.rank), 
#' get the matrix consist with 0 and 1, indicated the position of geneset, required for calculating ES score of the geneset.
#' 
#' @param  superposition.rank A matrix indicating rank of all values in the raw matrix.
#' @param  gSet.position A matrix consist of '0' and '1' show the position of signature genes in expression profile, 
#' with same size of expression profile, '1' stands for signature genes expression, '0' for other genes. 
#' @param  m Numeric value indicating number of row of matrix 'superposition.rank'.
#' @param  n Numeric value indicating number of column of matrix 'superposition.rank'.
#' @returnType matrix
#' @return 
 
indicatorInsideGeneSetMatrix <- function(superposition.rank, gSet.position, m, n){
  
  indicatorFunInsideGeneSet <- numeric(length = m*n)
  indicatorFunInsideGeneSet[superposition.rank] <- gSet.position
  
  indicatorFunInsideGeneSet <- matrix(indicatorFunInsideGeneSet, nrow = m)
  
  indicatorFunInsideGeneSet <- apply(indicatorFunInsideGeneSet, 2, rev)
  
  return(indicatorFunInsideGeneSet)
}



#' @@ oneSetssGSEA.ES
#' @description Calculating ES score of one geneset.
#' @param  exp.rank A matrix indicating rank of every element in column.
#' @param  superposition.rank.es A matrix indicating rank of all values in the matrix.
#' @param  gSet.position.matrix A matrix consist of '0' and '1' show the position of signature genes in expression profile, 
#' with same size of expression profile, '1' stands for signature genes expression, '0' for other genes. 
#' @param  alpha=0.25 Quantity alpha usually set to 0.25.
#' @param  m.num Numeric value indicating number of row of matrix 'exp.rank'.
#' @param  n.num Numeric value indicating number of column of matrix 'exp.rank'.
#' @returnType vector
#' @return vector of numeric values indicating ES score of one geneset for every sample
#' 
#' @author Liwen Xu
oneSetssGSEA.ES <- function(exp.rank, superposition.rank.es, gSet.position.matrix, alpha=0.25, m.num, n.num){

  indicatorInsideGeneSet.matrix <- indicatorInsideGeneSetMatrix(superposition.rank = superposition.rank.es, 
                                                                gSet.position = gSet.position.matrix, 
                                                                m=m.num, n=n.num)
  
  signature.temp.rank <- exp.rank * gSet.position.matrix
  signature.temp.rank <- signature.temp.rank[-which(rowSums(signature.temp.rank) == 0), ]
  signature.rank <- sortMatrixByCol(signature.temp.rank, set.decreasing=TRUE)
  
  score.for.oneset <- sapply(1:n.num, function(i){
    t4 <- indicatorInsideGeneSet.matrix[,i]
    t4[which(t4 != 0)] <- signature.rank[,i]
    stepCDFinGeneSet <- cumsum((abs(t4))^alpha)/ sum((abs(signature.rank[,i]))^alpha)
    t3 <- indicatorInsideGeneSet.matrix[,i]
    stepCDFoutGeneSet <- cumsum(!t3)/sum(!t3) 
    walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
    return(sum(walkStat))
  })
  
}
