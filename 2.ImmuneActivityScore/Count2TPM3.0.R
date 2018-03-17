#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'  导入gene_length.Rdata基因长度注释文件，第一列为基因名称，第二列为基因长度
#' @param countName counts file name
#' @param countPath count file address
#' @param genetype "gene symble" "entrez id" "ensemble id"
#' @param filePath genename_length2.0.Rdata 所在地址


count2Tpm<-function(countName,filePath,sample="multiple"){
  load(sprintf("%s/genename_length3.0.Rdata",filePath));
  
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
#运行


