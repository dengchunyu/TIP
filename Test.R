#' @title The test of running code of TIP 
#' 
#' @description 
#' You only should run this scripts if you want to test the code of TIP.
#' This scripts will show you how to run the function about different modules using example data.
#' There are three sections for TIP to calculate. We will run them in order:
#' 1.Check error
#' 2.Immunity Score
#' 3.Immune Infiltration(CIBERSORT)


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

#' @description 
#' This section is to calculate the proportion of infiltration immune cells use CIBERSORT.
#' This section will produce the file of infiltration proportion matrix including all samples.
#' For the operational efficiency of TIP,we set 100 times permutation, and have no options for other times, 
#' but we will improve it in the future.
#' For RNA-seq data, we will produce the results of 14 cell types, and 22 cell types for Micorarry, and
#' we hope we will provide more types of signature matrix for users' need in the future.
#' 
#' @details 
#' It should runs after 'immunityScore.R' because one of its results are the input data of 'CIBERSORT.server'.

source("./3.ImmuneInfiltration/1.CIBERSORT.server.R")
CIBERSORT.server(codePath="./3.ImmuneInfiltration",
                 filePath="./example_data_result",
                 signaturePath="./3.ImmuneInfiltration",
                 saveDir="./example_data_result",
                 perm = 10,
                 CHIPorRNASEQ="RNA-seq",
                 sample="multiple");
