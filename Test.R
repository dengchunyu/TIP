#This script is for quickly Test
source(file="./1.CheckError/checkError.R");
source(file="./1.CheckError/DealWithError.R");


#test
checkError(filePath = "./example_data", fileName = "RNA-seq_tpm_example_5.txt", format.of.file = "TPM", 
           codePath = "./2.ImmuneActivityScore", saveDir="./example_data_result")

example <- get(load("./example_data_result/expression.afterIDConvert.RData"))
sample.number <- ncol(example)

source("./2.ImmuneActivityScore/immunityScore.server.R")
immunityScore.server(codePath = "./2.ImmuneActivityScore", 
                     filePath = "./example_data_result", 
                     saveDir = "./example_data_result",
                     sampleNumber=sample.number,
                     permTimes = 10,
                     type.of.data = "RNA-seq",
                     format.of.file = "TPM");

source("./3.ImmuneInfiltration/1.CIBERSORT.server.R")
CIBERSORT.server(codePath="./3.ImmuneInfiltration",
                 filePath="./example_data_result",
                 signaturePath="./3.ImmuneInfiltration",
                 saveDir="./example_data_result",
                 perm = 10,
                 CHIPorRNASEQ="RNA-seq",
                 sample="multiple");
