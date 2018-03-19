
#' @@ ssGSEA.MultipleSamples.Permutation
#' @description  Perturb the rownames of an expression profile in order to obtain a random expression profile.
#'
#' @param  exp.profile A data.frame indicating the real expression profile.
#' @returnType data.frame
#' @return data.frame of random expression profile

ssGSEA.MultipleSamples.Permutation <- function(exp.profile){
  permutation.geneTags <- sample(rownames(exp.profile), size = nrow(exp.profile), replace = FALSE)
  permutation.profile <- exp.profile[permutation.geneTags, ]
  return(permutation.profile)
}
