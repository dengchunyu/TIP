filterOutliers <- function(raw.matrix){
  
  print("raw matrix summary before filter outliers: ")
  print(summary(as.vector(raw.matrix)))
  
  out <- boxplot.stats(as.vector(raw.matrix))
  print("stats:")
  print(out$stats)
  
  threshold_up <- out$stats[5]
  threshold_low <- out$stats[1]
  
  filter.number <- length(which(raw.matrix > threshold_up))
  filter.number <- filter.number + length(which(raw.matrix < threshold_low))

  print(paste0("filter.number:  ", filter.number/length(raw.matrix)))
  
  raw.matrix[which(raw.matrix > threshold_up)] <- threshold_up
  raw.matrix[which(raw.matrix < threshold_low)] <- threshold_low
  
  return(raw.matrix)
}