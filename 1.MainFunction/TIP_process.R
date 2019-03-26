#' @@ TIP_process
#' @description Entire calculation in TIP.
#' @param  codePath Character represting the storage path of import code files and underlying .RData.
#' @param  filePath Character represting the storage path of the users expression profile.
#' @param  fileName Character representing the file name of expression profile.
#' @param  saveDir Character represting the storage path of results.
#' @param  sampleNumber Numeric value indicating number of samples in profile.
#' @param  permTimes Numeric value indicating times of permutation, default by 100.
#' @param  type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param  format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'. 
#' @param  sample Character indicating sample number, choose 'multiple' or 'single'.

TIP_process <- function(
  codePath,
  filePath,
  fileName,
  saveDir,
  sampleNumber,
  permTimes,
  type.of.data,
  format.of.file,
  sample)
  {
  print("start to source function!")
  print(paste(codePath,"3.ImmuneInfiltration/CIBERSORT.server.R",sep="/"))
  print(paste(codePath,"2.ImmuneActivityScore/immunityScore.server.R",sep="/"))
  source(paste(codePath,"3.ImmuneInfiltration/CIBERSORT.server.R",sep="/"))
  source(paste(codePath,"2.ImmuneActivityScore/immunityScore.server.R",sep="/"))
  
  print("start to run immunityScore.server!")
  immunityScore.server(codePath =paste(codePath,"2.ImmuneActivityScore",sep="/"), 
                       filePath = filePath,
                       fileName = fileName,
                       saveDir = filePath,
                       sampleNumber = sampleNumber,
                       permTimes = permTimes,
                       type.of.data = type.of.data,
                       format.of.file = format.of.file);
  
  print("start to run CIBERSORT.server!")
  CIBERSORT.server(codePath =paste(codePath,"3.ImmuneInfiltration",sep="/"),
                   filePath = filePath,
                   fileName = fileName,
                   signaturePath = paste(codePath,"3.ImmuneInfiltration",sep="/"),
                   saveDir = filePath,
                   perm = permTimes,
                   CHIPorRNASEQ = type.of.data,
                   sample= sample,
                   dataType = format.of.file)
  
  print("the TIP_integration is end !")
}
