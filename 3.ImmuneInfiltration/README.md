# Immune Infiltration Component
## [Code]

**1.CIBERSORT_server.R**  
&ensp;&ensp;The final execute function of Infiltration. Reference all of other functions to produce results supporting web server.  

**2.CIBERSORT_main.R**  
&ensp;&ensp;Tranform the CIBERSORT infiltration results to interactive graphics data required by the TIP.  
&ensp;&ensp;It produces three types of files:   
&ensp;&ensp;*CIBER_barplot.txt*: Barplot interactive graphics format files.  
&ensp;&ensp;*CIBER_SampleName_lm14_pie.txt*: Pie interactive graphics format files.  
&ensp;&ensp;*CIBER_SampleName_lm14_radar.txt*: Radar plot interactive graphics format files.  

**3.CIBERSORT_func.R**    
&ensp;&ensp;The most important code about calculation of immune infiltration obtained from CIBERSORT v1.03.  
&ensp;&ensp;We downloaded the source code of CIBERSORT from https://cibersort.stanford.edu/   
&ensp;&ensp;Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)  
&ensp;&ensp;We make a slight modification on the source code of CIBERSORT, and we already obtain permission.  
&ensp;&ensp;License: http://cibersort.stanford.edu/CIBERSORT_License.txt  

&ensp;&ensp;When CIBERSORT has completed its run, a results matrix txt file will be produced in the file *CIBER_User_Allsample_Result.txt*.  
&ensp;&ensp;Columns represent cell types from the signature matrix file and rows represent deconvolution results for each sample. All results are reported as relative fractions normalized to 1 across all cell subsets.  
**The following metrics are provided for each sample:**  
&ensp;&ensp;*P-value*: Statistical significance of the deconvolution result across all cell subsets; p-value threshold default by 0.01. Increase the number of permutations (Basic Configuration Form) to increase the number of significant digits.  
&ensp;&ensp;*Correlation*: Pearson's correlation coefficient (R), generated from comparing the original mixture with the estimated mixture, the latter of which is calculated using imputed cell fractions and corresponding expression profiles from the signature genes file. Of note, the correlation is restricted to signature genes.  
&ensp;&ensp;*RMSE*: Root mean squared error between the original mixture and the imputed mixture, restricted to genes in the signature gene file.  

**4.CoreAlg.R**  
&ensp;&ensp;SVR core algorithm of CIBERSORT,obtained from CIBERSORT v1.03. NO change from source code. Referenced by CIBERSORT_func.R.   

**5.doPerm.R**  
&ensp;&ensp;Calculating permutations,obtained from CIBERSORT v1.03.NO change from source code. Referenced by CIBERSORT_func.R.    


## [File]
**1.LM14_name3.0.txt**  
&ensp;&ensp;14 immune cells signature matrix, background document used to caculate immune infiltration component of 14 immune cells. Calculated by single cell data of PBMC from GENOMIX 10X.  

**2.LM22_name.txt**  
&ensp;&ensp;22 immune cells signature matrix, background document used to caculate immune infiltration component of 22 immune cells. Obtained from CIBERSORT website.  
