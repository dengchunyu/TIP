#' @@ immunityScore.server
#' 
#' @description The main function to calculating immune activity levels based on ssGSEA algorithm,   
#' both for multiple samples and single sample profile.
#' 
#' @param codePath Character represting the storage path of import code files and underlying .RData.
#' @param filePath Character represting the storage path of the users expression profile.
#' @param fileName Character representing the file name of expression profile.
#' @param saveDir Character represting the storage path of results.
#' @param sampleNumber Numeric value indicating number of samples in profile.
#' @param permTimes Numeric value indicating times of permutation, default by 100.
#' @param type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'. 
#' @param format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @note R package 'pheatmap' is required.
#' @returnType 
#' @return NULL
#' 
#' @author Liwen Xu

immunityScore.server <- function(codePath, filePath, fileName, saveDir, sampleNumber, permTimes, type.of.data, format.of.file){
	
	
	
  # Import the R files and data next to use.
	load(sprintf("%s/signature.GeneSymbol.list.RData",codePath));
	load(sprintf("%s/annotation.RData",codePath));
	load(sprintf("%s/genename_length3.0.Rdata",codePath)); #length of genes for normalization from gene count to TPM
	source(sprintf("%s/Count2TPM3.0.R",codePath));
  
	source(sprintf("%s/ProcessMultipleSample.R",codePath));#main function for analysising expression profile with multiple samples 
	source(sprintf("%s/ProcessSingleSample.R",codePath)); #main function for analysing expression profile with only one sample
	
	source(sprintf("%s/ssgsea.core.R",codePath)); 
	source(sprintf("%s/ssGSEAPermutation.R",codePath));
	
  source(sprintf("%s/makeHeatmapData.R",codePath)); #change the format of results for visualization
	source(sprintf("%s/PCAScatterPlot.R",codePath)); 
	source(sprintf("%s/score_boxplot.R",codePath));
	source(sprintf("%s/changeStepNames.R",codePath));
	
  source(sprintf("%s/filterOutliers.R",codePath)); 
  
	print("source is end");
	print("set.seed");
	set.seed(1:100);
	
	library(pheatmap);
	
  if(sampleNumber > 1){
	  example <- get(load(sprintf("%s/%s",saveDir, "expression.afterIDConvert.RData")))
    print("read down!")
	  process.check <- tryCatch({
	    example.result <- ProcessMultipleSample(expression.from.users = example, 
	                                            save.dir = saveDir, 
	                                            signatureList = signature.GeneSymbol.list, 
	                                            perm.times = permTimes, 
	                                            signature.annotation = annotation,
	                                            type.of.data=type.of.data,
	                                            format.of.file=format.of.file,
	                                            gene.length.path = codePath)
	  }, error=function(e) {
	    errorString <- "Other"
	    cat("Processing wrong...\n")
	    write(errorString, file = sprintf("%s/ErrorString.txt",saveDir))
	    #save(e, file = sprintf("%s/ErrorInssGSEA.RData", saveDir))
	    write(e$message, file = sprintf("%s/ErrorInssGSEA.txt", saveDir))   
	  })
	  						
	}else{
	  example <- get(load(sprintf("%s/%s",saveDir, "expression.afterIDConvert.RData")))
	  print("read down!")
	  process.check <- tryCatch({
	    example.result <- ProcessSingleSample(expression.from.users = example, 
	                                          save.dir = saveDir, 
	                                          signatureList = signature.GeneSymbol.list, 
	                                          perm.times = permTimes, 
	                                          signature.annotation = annotation,
	                                          type.of.data=type.of.data,
	                                          format.of.file=format.of.file,
	                                          gene.length.path=codePath)
	   
	  }, error=function(e) {
	    errorString <- "Other"
	    cat("ssGSEA processing wrong...\n")
	    write(errorString, file = sprintf("%s/ErrorString.txt",saveDir))
	    write(e$message, file = sprintf("%s/ErrorInssGSEA.txt", saveDir))
	       
	  })
	}
  print("immunityScore.server is end");
  write("immunityScore.server", file=sprintf("%s/processResult.txt",saveDir),append=T);
  write(format(Sys.time(), "%Y.%m.%d.%H.%M.%S"), file=sprintf("%s/processResult.txt",saveDir),append=T);
    
}