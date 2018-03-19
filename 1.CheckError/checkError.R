#' @@ checkError
#' @description A main function of checking 
#' @param  filePath character representing the memory address of expression profile
#' @param  fileName character representing the file name of expression profile
#' @param  codePath character representing memory address of underlying .RData
#' @param  format.of.file character indicating the format of RNA-seq data, choose 'TPM' or 'Count' 
#' @param  saveDir character representing the memory address of the txt file contains error string return from 'DealWithError'.

checkError <- function(filePath, fileName, codePath, format.of.file, saveDir){
  result <- DealWithError(filePath = filePath, fileNam = fileNamee, codePath = codePath, format.of.file = format.of.file, saveDir=saveDir)
  print(result)
  if(!is.null(result)){
    if(!dir.exists(saveDir)){ dir.create(saveDir) }
    write(result, file = sprintf("%s/ErrorString.txt", saveDir))
  }
}