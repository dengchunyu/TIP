# TIP: Tracking Tumor Immunophenotype  
## Discription: 
TIP is a user-friendly one-stop shop web tool to comprehensively resolve tumor immunophenotype. This repository deposits code of underlying softwares of TIP platform(http://biocc.hrbmu.edu.cn/TIP/).

## The code is organized into a few different directories, each with a theme:  
**[1.MainFunction]**  
This folder contains the main function for entire calculation in TIP and integrated function for communication between R and web interfaces. There are also process for checking basic error and message logging.     

**[2.ImmuneActivityScore]**  
This folder contains twelve core source code files and four underlying data files used for the calculation of immune activity levels.  

**[3.ImmuneInfiltration]**  
This folder contains the core source code used for the calculation of immune infiltration.  

**[Test.R]**  
This script is for quick test of example data.  

**[Test]**  
This folder contains expression data from RNA-seq TPM for example and main calculation results. When you test by yourself, you will also get many other result files just for visualization on web, please ignore them.

## Citation  
- L. Xu, C. Deng, B. Pang, X. Zhang, W. Liu, G. Liao, H. Yuan, P. Cheng, F. Li, Z. Long, M. Yan, T. Zhao, Y. Xiao, and X. Li, *'Tip: A Web Server for Resolving Tumor Immunophenotype Profiling'*, Cancer Res, 78 (2018), 6575-80.  


