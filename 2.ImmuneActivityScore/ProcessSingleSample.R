#' @@ ProcessSingleSample
#' 
#' @description Several steps to deal with users uploaded single-sample expression profile:    
#' calculating and normalizing ES score of the sample, combine the normalized ES score from   
#' positive set and negative set of genes belong to the same step by calculating the difference between them,    
#' then change the format of results for visualization on web.
#' 
#' @param expression.from.users A data.frame contains expression data from users.  
#' @param save.dir Character represting the storage path of results.  
#' @param signatureList List of signature genes, default by 'signature.GeneSymbol.list'.
#' @param perm.times Numeric value indicating times of permutation, default by 100.  
#' @param signature.annotation A data.frame contains step and function direction information about signature genes,   
#' default by 'annotation'.
#' @param type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.  
#' @param format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.  
#' @param gene.length.path Character represting the storage path of the gene length file. 

ProcessSingleSample <- function(expression.from.users, save.dir, 
                                signatureList = signature.GeneSymbol.list, 
                                perm.times = 100, 
                                signature.annotation = annotation,
                                type.of.data = c("Microarray", "RNA-seq"), 
                                format.of.file = c("TPM","Count"),
                                gene.length.path){
  
  if(!dir.exists(save.dir)){ dir.create(save.dir) }
  ##########
  ##---------1. signature match(at least 1----------------------------------------------------------------------
  ##########
  #1. intersect number of every signature sets
  separate.remain <- sapply(signatureList, function(i){
    num <- length(intersect(i, rownames(expression.from.users)))
    return(num)
  })
  
  show.num <- paste0("(", separate.remain, "/", sapply(signatureList, length), ")")
  names(show.num) <- names(signatureList)
  
  #2. combine the positive genes and negative genes
  positive.l <- grep("positive", names(signatureList))
  #mix.show.num is for visualization
  mix.show.num <- data.frame(separate.remain = separate.remain, signature.num = sapply(signatureList, length),
                             step = c(1,2,2,3,3,4:20,21,21,22,22,23,23))
  
  #3.if one direction of a step matched 0 genes, the other direction will be also discarded
  l0 <- which(separate.remain == 0)
  a <- which(l0 %in% positive.l) 
  b <- which(l0 %in% (positive.l-1)) 
  if(length(a > 0) | length(b > 0)){
    delete.signature <- c(l0[-c(a,b)],l0[a],l0[a]-1,l0[b],l0[b]+1)
  }else{
    delete.signature <- l0
  }
  
  if(length(delete.signature) > 0){
    signatureList <- signatureList[-delete.signature]
    
    mix.show.num <- mix.show.num[-delete.signature,]
   
    get.separate <- unlist(strsplit(names(delete.signature),"\\."))
    delete.step <- intersect(get.separate[grep("Step", get.separate)], c("Step2","Step3","Step5","Step6","Step7"))
    delete.step <- gsub("Step","",delete.step)
    #4.delete needless siganture sets
    if(length(delete.step) > 0){
      signature.annotation <- signature.annotation[-which(signature.annotation[,4] %in% delete.step), ]
    }
    
  } 
  print(delete.signature)
  print(show.num)
  
  mix.show.num$step <- as.factor(mix.show.num$step)
  mix.show.num <- paste0("(", tapply(mix.show.num[,1], mix.show.num$step, sum), "/", 
                         tapply(mix.show.num[,2], mix.show.num$step, sum), ")")
  positive.l <- grep("positive", names(signatureList))
  if(length(positive.l) > 0){
    names(mix.show.num) <- names(signatureList)[-positive.l]
    names(mix.show.num) <- gsub(".negative", "", names(mix.show.num))
  }else{
    names(mix.show.num) <- names(signatureList)
  }
  
  mix.show.num <- paste(names(mix.show.num), mix.show.num, sep = " ")
  
  
  
  ##########
  ##---------2. ssGSEA process---------------------------------------------------------------------------------------------------
  ##########
  #1. Data preprocessing 
  #log2 transformation of expression data, and plus 1 for all values to make accurate rank later
  if(type.of.data == "RNA-seq" & format.of.file == "TPM"){
    save(expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)
    
  }else if(type.of.data == "RNA-seq" &  format.of.file == "Count"){
    #transform count to tpm
    normalized.expression.from.users <- count2Tpm(countName = expression.from.users, filePath = gene.length.path, sample="single")
    print(paste("count to TPM down---dim(normalized.expression.from.users)", dim(normalized.expression.from.users)))
    save(normalized.expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    normalized.expression.from.users <- as.matrix(log2(normalized.expression.from.users+1)+1)
    
  }else{  #if the microarray data was log-transformed
    save(expression.from.users, file = paste0(save.dir, "/expression.from.users.tpm.RData"))
    if(max(expression.from.users) < 50){
      normalized.expression.from.users <- as.matrix(expression.from.users+1)
    }else{
      normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)
    }
  }
  
  #2. build the random matrix 'permutation.exp', permutation of matrix 'expression.from.users'
  permutation.exp <- matrix(rep(normalized.expression.from.users, times = (perm.times+1)), 
                            nrow = nrow(normalized.expression.from.users))
  for(i in 1:perm.times){
    permutation.exp[, i+1] <- sample(normalized.expression.from.users, size = nrow(normalized.expression.from.users), 
                                     replace = FALSE)
  }
  rownames(permutation.exp) <- rownames(normalized.expression.from.users)
  
  #3.ssGSEA for the whole matrix 'permutation.exp'
  exp.n <- ncol(permutation.exp)
  exp.m <- nrow(permutation.exp)
  #(1)get the superimposed matrix
  superposition.rank.matrix <- superpositionRank(raw.matrix = permutation.exp, m = exp.m, n = exp.n)
  #(2)get rank by column of matrix 'permutation.exp'
  permutation.exp.rank <- rankMatrixByCol(superposition.rank = superposition.rank.matrix, m = exp.m, n = exp.n)
  
  #activity score of signature sets one by one
  permutation.score <- t(sapply(signatureList, function(sig){
    gSetIdx <- which(rownames(permutation.exp) %in% sig)
    #(3)vector of '0'and '1' show the position of signature genes in 'permutation.exp'
    t1 <- numeric(exp.m)
    t1[gSetIdx] <- 1
    geneSet.position <- matrix(rep(t1, times = exp.n), ncol = exp.n)
    #(4)activity score of one signature set
    score.for.oneset <- oneSetssGSEA.ES(exp.rank = permutation.exp.rank, superposition.rank.es = superposition.rank.matrix, 
                                        gSet.position.matrix = geneSet.position, alpha=0.25, m.num = exp.m, n.num = exp.n)
    
    return(score.for.oneset)
  }))
  print("permutation.score over!")
  
  
  ##########
  ##---------3.normalization of activity score--------------------------------------------------------
  ##########
  #1. whether activity score of real sample and permutation samples(N+1，2N+1 ... perm.times*N+1) are of the same sign
  #if different, set the score of permutation samples to 'NA'
  permutation.score <- t(apply(permutation.score, 1, function(row){
    row[which(row[1] * row < 0)] <- NA
    return(row)
  }))
  
  #2. z-score normalization
  zScore <- apply(permutation.score, 1, function(x) {
    if(length(which(!is.na(x))) < 2){
      z.score <- x[1]
    }else{
      z.score <- scale(x)[1]
    }
    return(z.score)
  })
  print("zScore over!")
  
  #3. step 2,3,5,6,7 z-score positive - negative
  positive.l <- grep("positive", names(signatureList))
  if(length(positive.l) > 0){
    #change the sign of score of negative genes 
    zScore[positive.l-1] <- zScore[positive.l-1]*(-1)
    direction <- 1:length(signatureList)
    direction[positive.l] <- direction[positive.l-1]
  }else{
    direction <- 1:length(signatureList)
  }
  print("direction over!")
  
  #conbine positive genes score and negative genes score
  ssGSEA.normalized.score <- as.vector(tapply(zScore, factor(direction), sum))
  if(length(positive.l) > 0){
    names(ssGSEA.normalized.score) <- names(signatureList)[-positive.l]
    names(ssGSEA.normalized.score) <- gsub(".negative", "", names(ssGSEA.normalized.score))
  }else{
    names(ssGSEA.normalized.score) <- names(signatureList)
  }
  
  ssGSEA.normalized.score <- round(ssGSEA.normalized.score, 3)
  print("ssGSEA.normalized.score over!")
  save(ssGSEA.normalized.score, file = paste0(save.dir, "/ssGSEA.normalized.score.RData"))
  
  
  
  ##########
  ##---------4. Change the format of results for visualization------------------------------------------------------------
  ##########
  #1.ssGSEA.normalized.score----------------------------------------------------------------------------------------------
  #(1)show.num
  final.mix.show.num <- paste('["',paste(rev(mix.show.num), collapse='","'), '"]', sep="")
  write(final.mix.show.num, file = paste0(save.dir, "/final.mix.show.num.txt"))
  print("final.mix.show.num over!")
  #(2)ssGSEA.normalized.score txt
  ssGSEA <- data.frame(Steps = names(ssGSEA.normalized.score), ssGSEA.normalized.score)
  colnames(ssGSEA)[-1] <- colnames(expression.from.users)
  write.table(ssGSEA, file = paste0(save.dir, "/ssGSEA.normalized.score.txt"), quote = FALSE, 
              row.names = FALSE, col.names = TRUE, sep = "\t")
  print("ssGSEA over!")
  ##2. visualization from an individual perspective-----------------------------------------------------------------------
  #(1)line chart
  #combine step4
  step4.lo <- grep("Step4", names(ssGSEA.normalized.score))
  if(length(step4.lo) > 0){
    indicator <- 1:length(ssGSEA.normalized.score)
    indicator[step4.lo] <- indicator[step4.lo[1]]
    step4.com.score <- tapply(ssGSEA.normalized.score, indicator, mean)
    step4.com.score <- as.matrix(round(step4.com.score, 3))
    rownames(step4.com.score) <- sort(c(names(ssGSEA.normalized.score)[-step4.lo], "Step4"))
  }else{
    step4.com.score <- ssGSEA.normalized.score
  }
  
  rownames(step4.com.score) <- changeStepNames(rownames(step4.com.score))
  colnames(step4.com.score) <- colnames(expression.from.users)
  nu <- BarPlotFormat(step4.com.score, save.dir = save.dir)
  print("BarPlotFormat over!")
  
  #(2)scores of every steps except step4
  a <- sapply(setdiff(1:length(ssGSEA.normalized.score), step4.lo), function(step){
    write.table(ssGSEA.normalized.score[step], file = paste0(save.dir, "/", colnames(expression.from.users), "_", 
                                                             names(ssGSEA.normalized.score)[step], ".txt"), 
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  })
  print("step123567 over!")
  #(3)radar chart
  if(length(step4.lo) > 0){
    radar <- data.frame(name = gsub("Step4.", "", names(ssGSEA.normalized.score)[step4.lo]),
                        sample = ssGSEA.normalized.score[step4.lo])
    write.table(t(radar), file = paste0(save.dir, "/", colnames(expression.from.users), "_Step4.txt"), 
                quote = FALSE, row.names = TRUE, col.names = FALSE, sep = ",")
  }
  
  print("radar over!")
  
  #3. expression of signature expression---------------------------------------------------------------------------
  anno <- unique(signature.annotation[,3:4])
  inter <- intersect(signature.annotation[,3], rownames(expression.from.users))
  print(paste0("intersect number: ", length(inter)))
  anno <- anno[which(anno[,1] %in% inter),]
  print(paste("dim(anno)---", dim(anno)))
  
  annotation_row_length <- tapply(anno[,2], anno[,2], length)
  print(paste("sum of annotation_row_length: ", sum(annotation_row_length)))
  
  com.expression <- cbind(expression.from.users, expression.from.users)
  expression.profile <- com.expression[anno[,1], ]
  expression.profile <- as.matrix(expression.profile[,1])
  colnames(expression.profile) <- colnames(expression.from.users)
  print(paste("dim(expression.profile): ",dim(expression.profile)))
  storage.mode(expression.profile) <- "numeric"
  save(expression.profile, file = paste0(save.dir, "/expression.profile.RData"))
  
  #(1)signature.expression txt
  heatmap.expression <- expression.profile
  rownames(heatmap.expression) = 1:nrow(heatmap.expression)
  
  signature.expression <- data.frame(GeneSymbol = anno[,1], Steps = anno[,2], round(heatmap.expression, 3))
  signature.expression[,1] <- as.character(signature.expression[,1])
  colnames(signature.expression)[-(1:2)] <- colnames(expression.from.users)
  write.table(signature.expression, file = paste0(save.dir, "/SignatureGenes.Expression.txt"), quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep = "\t")
  print("signature.expression over!")

  #(2)signature expression for individual sample(box plot)
  box.expression <- signature.expression
  box.expression[, -(1:2)] <- log2(box.expression[, -(1:2)]+1)
  Score_boxplot(Exp.class.matrix = box.expression, saveDir = save.dir)
  print("boxplot over!")
  

  ##########
  ##---------5.return---------------------------------------------------------------------------------------------------
  ##########
  result <- list(ssGSEA.normalized.score = ssGSEA.normalized.score, signature.expression = expression.profile, 
               separate.remain = separate.remain, 
               show.num = show.num, delete.signature = delete.signature, mix.show.num = mix.show.num)
  save(result, file = paste0(save.dir, "/multiSample.result.RData"))
  return(list(separate.remain = separate.remain, show.num = show.num, delete.signature = delete.signature,
            ssGSEA.normalized.score = ssGSEA.normalized.score, signature.expression = expression.profile))
}


