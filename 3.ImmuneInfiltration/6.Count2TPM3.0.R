#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'  Import 'the annotation file of gene length :'gene_length.Rdata'  
#'  Gene name in the first column ,and gene length in the second column.
#'  
#'
#' @param countName counts file name
#' @param sample "single" or "multiple" samples you upload
#' @param filePath The storage path of genename_length2.0.Rdata 


count2Tpm<-function(countName,filePath,sample="multiple"){
  
   # import gene length file
  load(sprintf("%s/genename_length3.0.Rdata",filePath));
  
  # Two processing path for the choice of parameter 'sample'
  
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
  
# The function to change
countToTpm <- function(counts, effLen)
{
  
  rate <- log(counts) - log(effLen)
  
  denom <- log(sum(exp(rate)))
  
  exp(rate - denom + log(1e6))
  
}
# produce results
data_TPM<-apply(count,2,function(x){
  
  countDf <- data.frame(count = as.integer(x), length = as.integer(Length))
  
  countDf$effLength <- countDf$length
  
  countDf$tpm <- with(countDf, countToTpm(count, effLength))
  
  return(countDf$tpm)
  
})

rownames(data_TPM)<-rownames(count)

return(data_TPM)

}
