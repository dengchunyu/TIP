`[scripts]`  
**DealWithError.R**  is the core function to judge the possible error within users uploaded file.  
**checkError.R**  is main function and export a txt contains error string return from 'DealWithError'.   

**Note:**  
	If encounter the warning message like 'NAs introduced by coercion', it doesn't means that there is something wrong in code，we designed to transfrom the colnames into numeric to check if there is no sample names.
