#' @@ checkError
#' invoking function 'DealWithError'   

checkError <- function(filePath = filePath, fileName = fileName, codePath=codePath, format.of.file=format.of.file, saveDir=saveDir){
  result <- DealWithError(filePath, fileName, codePath, format.of.file, saveDir=saveDir)
  print(result)
  if(!is.null(result)){
    if(!dir.exists(saveDir)){ dir.create(saveDir) }
    write(result, file = sprintf("%s/ErrorString.txt", saveDir))
  }
}