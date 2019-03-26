#' @@ processResult.server
#' @description Message logging of users uploaded expression profile.
#' @param  filePath Character represting the storage path of the users expression profile.
#' @param  CancerType The type of cancer selected by the user.
#' @param  type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param  format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @param  sample Character indicating sample number, choose 'multiple' or 'single'.
#' @param  sampleNumber Numeric value indicating number of samples in profile.
#' @param  Samples A string consists all names of samples seperated with tab.
#' @param  email E-mail address of user (optional).

processResult.server <- function(filePath, CancerType, type.of.data, format.of.file, sample, sampleNumber, Samples, email){

	if(email==""){email = "no Email";}
	inputInfor <- c(CancerType, type.of.data, format.of.file, sample,sampleNumber,Samples, email, format(Sys.time(), "%Y.%m.%d.%H.%M.%S"));
	
	write.table(inputInfor, file=sprintf("%s/processResult.txt",filePath),quote = FALSE, row.names = FALSE, col.names= FALSE);
	print("All R compute finishes successfully.");
}
