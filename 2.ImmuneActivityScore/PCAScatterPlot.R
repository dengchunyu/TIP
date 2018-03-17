PCASatterPlot <- function(pca.matrix, save.file){
  each.spot <- paste0('{"x":',pca.matrix[,1], ',"y":', pca.matrix[,2],
                      ',"name":"', rownames(pca.matrix),'"}')
  total.spots <- paste0("[", paste(each.spot, collapse = ","), "]")
  write(total.spots, file = save.file)
}
