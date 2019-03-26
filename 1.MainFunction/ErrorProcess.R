#' @@ DealWithError
#' @description A function for checking some basic error of users uploaded expression profile.
#' When the format isn't standard, return a string indicating the corresbonding error type. 
#' Required to function checkError.
#' 
#' @param  filePath Character represting the storage path of the users expression profile.
#' @param  fileName Character representing the file name of expression profile.
#' @param  codePath Character represting the storage path of import code files and underlying .RData.
#' @param  format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @param  type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param  saveDir If the gene identifier were thansformed, the new expression profile will be save in the directory.
#' @return NULL or character

DealWithError <- function(filePath, fileName, codePath, format.of.file, type.of.data, saveDir){
  #1、error1：InputError --------------------------------------------------------------------
  input.check <- tryCatch({
    expression.profile <- read.table(sprintf("%s/%s",filePath, fileName), sep = "\t", stringsAsFactors = FALSE, header =TRUE, 
                                     row.names = 1, check.names=F, na.strings = NULL)
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
                                       check.names=F, na.strings = NULL)
      if(sum(duplicated(expression.profile[,1])) > 0){errorString <- "DuplicateNameError"; return(errorString);}
    }, error=function(e){
      print(e$message)
      
      return(NULL)    
    })
    
    if(is.null(input.check2)){errorString <- "InputError"; return(errorString);}
    
    
  }else{
       SampleNumber.check <- ncol(expression.profile)
       
  #error2: SeparatorError Separators are not "\t" ---------------------------------------------------------
       if(SampleNumber.check == 0){errorString <- "SeparatorError"; return(errorString);}
  #error3: geneNameError Gene identifiers don't match -----------------------------------------------------
       #Transform gene identifer
       gene.info <- get(load(sprintf("%s/2.ImmuneActivityScore/human_gene2ensembl2symbol_list.RData", codePath)));
       id <- as.character(rownames(expression.profile))
       names(id) <- id
       
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
       
       load(sprintf("%s/2.ImmuneActivityScore/annotation.RData",codePath));
       overlap.num <- length(intersect(rownames(expression.profile), annotation[,3]))
       if(overlap.num < 20){errorString <- "GeneNameError"; return(errorString);}
  #error4: HeaderNameError No colnames --------------------------------------------------------------------
      if(all(!is.na(as.numeric(colnames(expression.profile)[-1])))){
        errorString <- "HeaderNameError"; return(errorString);
  #error5: DefaultValueError If there is any default value in profile -------------------------------------
      }else if(sum(complete.cases(expression.profile))!=nrow(expression.profile) | sum(is.na(expression.profile))!=0 | 
         length(which(expression.profile==""))!=0 | length(which(expression.profile==" "))!=0){
        errorString <- "DefaultValueError"; return(errorString);
      }
  }
  if(!dir.exists(saveDir)){ dir.create(saveDir) }
  save(expression.profile, file = paste0(saveDir, "/expression.afterIDConvert.RData"))
  return(NULL)
}



#' @@ checkError
#' @description A main function of checking error.
#' @param  filePath Character represting the storage path of the users expression profile.
#' @param  fileName Character representing the file name of expression profile.
#' @param  codePath Character represting the storage path of import code files and underlying .RData.
#' @param  format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @param  type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param  saveDir Character representing the storage path of the txt file contains error string return from 'DealWithError'.

checkError <- function(filePath, fileName, codePath, format.of.file, type.of.data, saveDir){
  result <- DealWithError(filePath = filePath, fileName = fileName, codePath = codePath, 
                          format.of.file = format.of.file, type.of.data = type.of.data, saveDir=saveDir)
  print(result)
  if(!is.null(result)){
    if(!dir.exists(saveDir)){ dir.create(saveDir) }
    write(result, file = sprintf("%s/ErrorString.txt",saveDir))
  }
}
