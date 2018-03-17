ssGSEA.MultipleSamples.Permutation <- function(exp.profile){
  permutation.geneTags <- sample(rownames(exp.profile), size = nrow(exp.profile), replace = FALSE)
  permutation.profile <- exp.profile[permutation.geneTags, ]
  return(permutation.profile)
}