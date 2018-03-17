#' @@ superpositionRank
#' @description get rank after superposition
superpositionRank <- function(raw.matrix, m, n){
  max.element <- max(raw.matrix)
  over.matrix <- matrix(rep(seq(0, max.element*(n-1), by = max.element), each=m), nrow=m)
  superposition <- raw.matrix + over.matrix
  
  superposition.rank <- matrix(rank(superposition, ties.method = "first"), nrow = m)
  
  return(superposition.rank)
}



#' @@ rankMatrixByCol
#' @description after superposition, get the rank of evert column of the raw matrix
rankMatrixByCol <- function(superposition.rank, m, n){
  minus.matrix <- matrix(rep(seq(0, m*(n-1), by=m), each=m), nrow=m)
  subtract.matrix <- superposition.rank - minus.matrix
  
  return(subtract.matrix)
}



#' @@ indicatorInsideGeneSetMatrix
#' @description after using rankMatrixByCol, get the matrix consist with 0 and 1, indicated the position of geneset
indicatorInsideGeneSetMatrix <- function(superposition.rank, gSet.position, m, n){
  
  indicatorFunInsideGeneSet <- numeric(length = m*n)
  indicatorFunInsideGeneSet[superposition.rank] <- gSet.position
  
  indicatorFunInsideGeneSet <- matrix(indicatorFunInsideGeneSet, nrow = m)
  
  indicatorFunInsideGeneSet <- apply(indicatorFunInsideGeneSet, 2, rev)
  
  return(indicatorFunInsideGeneSet)
}



#' @@ sortMatrixByCol
#' @description quick way to sort every column of a matrix 
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



#' @@ ssGSEA
#' @description ES score of one set
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
