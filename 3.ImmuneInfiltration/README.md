# TIP
## Ｃode
The TIP web server source code including two sections,immune activity score and CIBERSORT immune infiltration,the code in this project is about CIBERSORT immune infiltration. 
### Ｔhe most important code about immune infiltration come from CIBERSORT v1.03,little revise
1.CoreAlg.R  
Core algorithm  

2.doPerm.R   
do permutations  

3.CIBERSORT_func.R   
Main functions for CIBERSORT computation.   

### Ｔhe code for transforming  
4.Count2TPM3.0.R   
Transform Count data to TPM data.    

### Ｔhe code for webserver results arrangement 
5.CIBERSORT_main.R  
change the CIBERSORT infiltration results to interactive graphics data that the TIP need.  

6.CIBERSORT_server.R  
the final execute function, all code files are on tap. export final results.  

## Ｆile
1.LM14_name3.0.txt  
14 immune cells signature matrix  

2.LM22_name.txt  
22 immune cells signature matrix  

3.genename_length3.0.Rdata  
gene length file for Count2TPM3.0.R  



