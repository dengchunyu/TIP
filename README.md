# TIP: Tracking Tumor Immunophenotype  
## Discription: 
TIP is a user-friendly one-stop shop web tool to comprehensively resolve tumor immunophenotype. This repository deposits code of underlying softwares of TIP platform(http://biocc.hrbmu.edu.cn/TIP/).

## The code is organized into a few different directories, each with a theme:  
**[1.CheckError]**  
Contains two code files of first step of TIP web-server, checking whether the format of upload file is standard or not.  

**[2.ImmuneActivityScore]**  
Contains twelve core source code files and four underlying data files used for the calculation of immune activity levels.  

**[3.ImmuneInfiltration]**  
Contains five core source code files and two signature matrix files used for the calculation of immune infiltration.  

**[Test]**  
R scripts for quick run of total process, if you want to test the source code of TIP, run it first.  

**[example_data]**  
contains expression data from RNA-seq TPM for example.

**[example_data_result]**  
The running results for example data: RNA-seq_tpm_example_5.txt.

**[CIBERSORT_license.md]**   
Because we will utilize the source code of CIBERSORT to caculate immune cell infiltration proportion. We obtained the license from https://cibersort.stanford.edu/
