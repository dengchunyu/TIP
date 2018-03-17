makeHeatmapData <- function(raw.matrix){
  m <- nrow(raw.matrix)
  n<- ncol(raw.matrix)
  heatmapData <- c()
  for(col in 1:n){
    for(row in 1:m){
      coordinate <- paste('[', col-1, ',' , m- row, ',', 
                          raw.matrix[row,col], ']', sep = "")
      heatmapData <- c(heatmapData, coordinate) 
    }
  }
  final <- paste('[',paste(heatmapData, collapse=','), ']', sep="")
  return(final)
}