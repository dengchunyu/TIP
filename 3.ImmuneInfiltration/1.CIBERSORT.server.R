#' Main functions
#' Author: DengChunyu, Haerbin Medical University, dengcyelena@gmail.com
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.packages('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script import 5 R files.
#' 
#'   Example In R:
#'       source('CIBERSORT.server.R')
#'       
#'       CIBERSORT.server(codePath="",
#'       filePath="",
#'       signaturePath="",
#'       saveDir="",
#'       perm = 100,
#'       CHIPorRNASEQ="Microarray",
#'       sample="multiple");
#' 
#' @param codePath: the storage path of import code files.
#' @param filePath ：the storage path of the users data.
#' @param fileName ：the name of the users file.this parameter is useless now after optimization.
#' @param signaturePath: the storage path of signature matrix files.
#' @param saveDir： the storage path of results.
#' @param perm： the NO. permutation.default=100
#' @param CHIPorRNASEQ: the type of data,RNAseq or Microarray.
#' @param sample： "multiple" or "single" sample for upload files.
#' @param dataType: "count" or "TPM" for upoad data.
#' library('e1071')
#' library('parallel')
#' library('preprocessCore')
#' @export 

CIBERSORT.server <- function(codePath, filePath, fileName= "",signaturePath, saveDir, perm = 100,CHIPorRNASEQ="RNA-seq", sample="multiple",dataType="TPM"){

  # Import the R files and data we next to use.
	
  source(sprintf("%s/2.CIBERSORT_main.R",codePath));
  source(sprintf("%s/3.CIBERSORT_func.R",codePath));
  source(sprintf("%s/4.CoreAlg.R",codePath));
  source(sprintf("%s/5.doPerm.R",codePath));
  source(sprintf("%s/6.Count2TPM3.0.R",codePath));
  
  print("The Infiltration source codes is import!");

	
	
  #Import the user's data,for RData format.because this file is from the previous section of Webserver.
	# matrix format
	
  mixture_file = get(load(paste0(filePath,"/expression.from.users.tpm.RData")))
  print("The user data is import!");
	
  #choose one type of signature matrix we next to use,RNA-seq data for LM14 matrix;Microarray data for LM22 matrix.
	
  if(CHIPorRNASEQ=="RNA-seq"){
    
    signature_name<-"LM14_name3.0.txt"
    
  }else if(CHIPorRNASEQ=="Microarray"){
    
    signature_name<-"LM22_name.txt"
    
    }
  
  #import the signature matrix.
  sig_matrix <- read.table(sprintf("%s/%s",signaturePath, signature_name), sep = "\t", stringsAsFactors = FALSE,header=T,row.names=1,check.names=F)
  
  print("The signature matrix is import!");
  
  # run the function to produce result files.
  Result <- CIBERSORT_server(mixture_file= mixture_file,
                             sig_matrix = sig_matrix ,
                             saveDir=saveDir,
                             perm=perm, 
                              CHIPorRNASEQ=CHIPorRNASEQ,
                             sample=sample)
  
  print("CIBERSORT.server is end");
	
  write("CIBERSORT.server", file=sprintf("%s/processResult.txt",saveDir),append=T);
}
