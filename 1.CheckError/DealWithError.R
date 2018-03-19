#' @@ DealWithError
#' @description A function for checking users uploaded expression profile, when the format isn't standard, 
#' return a string indicating the corresbonding error type. 
#' Required to function checkError.
#' @param  filePath character representing the memory address of expression profile
#' @param  fileName character representing the file name of expression profile
#' @param  codePath character representing memory address of underlying .RData
#' @param  format.of.file character indicating the format of RNA-seq data, choose "TPM" or "Count" 
#' @param  saveDir If the gene identifier were thansformed, the new expression profile will be save in the directory.
#' @return NULL or character

DealWithError <- function(filePath, fileName, codePath, format.of.file, 
                          saveDir){
  #error1: InputError --------------------------------------------------------------------
  input.check <- tryCatch({
    expression.profile <- read.table(sprintf("%s/%s",filePath, fileName), sep = "\t", stringsAsFactors = FALSE, header =TRUE, 
                                     row.names = 1, check.names=F)
   # return("success")
  }, error=function(e){
    print(e$message)
    
    return(NULL)    
  })
  
  
  ##Input fail.
  if(is.null(input.check)){
    #Try again.
    input.check2 <- tryCatch({
      expression.profile <- read.table(sprintf("%s/%s",filePath, fileName), sep = "\t", stringsAsFactors = FALSE, header =TRUE, 
                                       check.names=F)
      if(sum(duplicated(expression.profile[,1])) > 0){errorString <- "DuplicateNameError"; return(errorString);}
    }, error=function(e){
      print(e$message)
      
      return(NULL)    
    })
    
    if(is.null(input.check2)){errorString <- "InputError"; return(errorString);}
    
   
   }else{
       SampleNumber.check <- ncol(expression.profile)
       
  #error2: SeparatorError Separators are not "\t"--------------------------------------------------------------
       if(SampleNumber.check == 0){errorString <- "SeparatorError"; return(errorString);}
  #error3: geneNameError Gene identifiers don't match----------------------------------------------------------
       #Transform ID
       gene.info <- get(load(sprintf("%s/human_gene2ensembl2symbol_list.Rdata",codePath)));
       id <- as.character(rownames(expression.profile))
       names(id) <- id
       #
       inter_id <- intersect(id, gene.info[[1]][,1])
       if(length(inter_id) > 0){
         id[inter_id] <- gene.info[[1]][inter_id, 2]
       }
       
       inter_id<-intersect(id, gene.info[[2]][,1])
       if(length(inter_id) > 0){
         id[inter_id] <- gene.info[[2]][inter_id, 2]
       }
       
       rownames(expression.profile) <- id
       print("Rownames to geneSymbol done!")
       
       load(sprintf("%s/annotation.RData",codePath));
       overlap.num <- length(intersect(rownames(expression.profile), annotation[,3]))
       if(overlap.num < 20){errorString <- "GeneNameError"; return(errorString);}
  #error4: HeaderNameError No colnames ----------------------------------------------------------------------------
      if(all(!is.na(as.numeric(colnames(expression.profile)[-1])))){
        errorString <- "HeaderNameError"; return(errorString);
  #error5: DefaultValueError If there is any default value in profile----------------------------------------------
      }else if(sum(complete.cases(expression.profile))!=nrow(expression.profile) | sum(is.na(expression.profile))!=0 | 
         length(which(expression.profile==""))!=0 | length(which(expression.profile==" "))!=0){
        errorString <- "DefaultValueError"; return(errorString);
      }
  #error6: RNASeqCountError If the 'format.of.file' is set to 'Count', but there is decimals------------------------
      Count.check <- all((round(expression.profile) == expression.profile) == TRUE)
      if(format.of.file == "Count" & !Count.check){errorString <- "RNASeqCountError";return(errorString);}
  
  }
  if(!dir.exists(saveDir)){ dir.create(saveDir) }
  save(expression.profile, file = paste0(saveDir, "/expression.afterIDConvert.RData"))
  return(NULL)
}
