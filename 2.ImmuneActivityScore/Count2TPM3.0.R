#' @title Transform the data type of RNA-seq from count to TPM
#' 
#' @description 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million(TPM).
#' 
#' @details 
#' We need 'genename_length2.0.Rdata' to caculate TPM data.
#' 
#' @param countName 
#' The name of counts file.
#' 
#' @param filePath 
#' The storage path of 'genename_length2.0.Rdata'.
#' 
#' @param sample 
#' Choose "multiple" or 'single' for sample number.


count2Tpm<-function(countName,filePath,sample="multiple"){
  load(sprintf("%s/genename_length3.0.Rdata",filePath));
  
  # Deal with sample matrix
  if(sample=="multiple"){
    
    count <- as.matrix(countName)
    inter<-intersect(rownames(count),genename_length[,2])
    count<-count[match(inter,rownames(count)),]
    Length<-genename_length[match(inter,genename_length[,2]),3]
    
  }else if(sample=="single"){
    colName<-colnames(countName)
    count <- as.matrix(countName)
    inter<-intersect(rownames(count),genename_length[,2])
    count<-matrix(count[match(inter,rownames(count)),],ncol = 1)
    rownames(count)<-inter
    colnames(count)<-colName
  }
  Length<-genename_length[match(inter,genename_length[,2]),3]
  
# function to transform data type.
countToTpm <- function(counts, effLen)
{
  
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
data_TPM<-apply(count,2,function(x){
  countDf <- data.frame(count = as.integer(x), length = as.integer(Length))
  countDf$effLength <- countDf$length
  countDf$tpm <- with(countDf, countToTpm(count, effLength))
  return(countDf$tpm)
})
rownames(data_TPM)<-rownames(count)
return(data_TPM)
}


