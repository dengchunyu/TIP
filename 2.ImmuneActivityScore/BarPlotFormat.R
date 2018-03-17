BarPlotFormat <- function(score.matrix, save.dir){
  #7 step scores of one sample save into a file
  for(col.num in 1:ncol(score.matrix)){
    separate.step <- paste0('{"name":"',rownames(score.matrix),'","value":', score.matrix[, col.num], "}")
    onesample <- paste0("[", paste(separate.step, collapse = ","), "]")
    write(onesample, file = paste0(save.dir,"/",colnames(score.matrix)[col.num], "_BarPlotFormat.txt"))
  }
}