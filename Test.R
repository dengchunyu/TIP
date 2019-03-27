#' @title The test of running code of TIP 
#' 
data <- read.table(sprintf("%s/%s","./example_data", "RNA-seq_tpm_example_5.txt"), 
                   sep = "\t", stringsAsFactors = FALSE, header =TRUE,
                   check.names=F, na.strings = NULL, row.names = 1) 

source("./1.MainFunction/TIP_integration.R")
test <- TIP_integration(
  codePath = ".",
  filePath = "./Test",
  fileName = "RNA-seq_tpm_example_5.txt",
  saveDir = "./Test",
  sampleNumber = 5,
  permTimes = 100,
  type.of.data = "RNA-seq",
  format.of.file = "TPM",
  sample = "multiple",
  CancerType = "GBM",
  Samples=paste(colnames(data), collapse = "\t"), email="")
