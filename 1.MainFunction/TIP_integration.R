#' @@ TIP_integration
#' @description Main function for communication between R and web interfaces.
#' @param  codePath Character represting the storage path of import code files and underlying .RData.
#' @param  filePath Character represting the storage path of the users expression profile.
#' @param  fileName Character representing the file name of expression profile.
#' @param  saveDir Character represting the storage path of results.
#' @param  sampleNumber Numeric value indicating number of samples in profile.
#' @param  permTimes Numeric value indicating times of permutation, default by 100.
#' @param  type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param  format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'. 
#' @param  sample Character indicating sample number, choose 'multiple' or 'single'. 
#' @param  CancerType The type of cancer selected by the user.
#' @param  Samples A string consists all names of samples seperated with tab.
#' @param  email E-mail address of user (optional).

TIP_integration <- function(
codePath,
filePath,
fileName,
saveDir,
sampleNumber,
permTimes,
type.of.data,
format.of.file,
sample,
CancerType, 
Samples, 
email){
  print("The filePath...")
  print(filePath)
  sink(file = sprintf("%s/%s",saveDir, "allPrint.txt"), append = TRUE)
  
  print("The filePath...")
  print(filePath)
  print("The saveDir...")
  print(saveDir)
  
  #1.checkerror
  print("start to run checkError!")
  process.check1 <- tryCatch({
    source(paste(codePath,"RCodeAndData/ErrorProcess.R",sep="/"))
    checkError(filePath = filePath, fileName = fileName, format.of.file = format.of.file, 
               type.of.data = type.of.data, codePath = codePath, saveDir = saveDir)
   
  }, error=function(e) {
    cat("checkError processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt",saveDir), append = TRUE)
    return("error")  
  })
  # process success
  if(!is.null(process.check1)){
    sink();
    stop();
  }else{
    cat("checkError processing over...\n")
    }
  
  
  #2.processResult.server
  process.check2 <- tryCatch({
    source(paste(codePath,"RCodeAndData/processResult.server.R",sep="/"))
    processResult.server(filePath, CancerType, type.of.data, format.of.file, sample,sampleNumber,Samples, email)
    
  }, error=function(e) {
    cat("processResult.server processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt", saveDir), append = TRUE)
    return("error")  
  })
  
  if(process.check2 == "error"){
    sink();
    stop();
  }else{
    cat("processResult.server over...\n")
  }
  
  #3.TIP_process
  process.check3 <- tryCatch({
    source(paste(codePath,"RCodeAndData/TIP_process.R",sep="/"))
    TIP_process(codePath, filePath, fileName, saveDir, sampleNumber, permTimes,
                type.of.data, format.of.file, sample)
    
  }, error=function(e) {
    cat("TIP_process processing wrong...\n")
    print(e$message)
    write(e$message, file = sprintf("%s/ErrorString.txt",saveDir),append = TRUE)
    return("error")  
  })
  if(process.check3 == "error"){
    sink();
    stop();
  }else{
    cat("TIP process over...\n")
  }
  
  sink()
}
  