changeStepNames <- function(name.vector){
  name.vector <- gsub(".*tep1", "Step1:Release of cancer antigens", name.vector)
  name.vector <- gsub(".*tep2", "Step2:Cancer antigen presentation", name.vector)
  name.vector <- gsub(".*tep3", "Step3:Priming and activation", name.vector)
  name.vector <- gsub(".*tep4", "Step4:Trafficking of immune cells to tumors", name.vector)
  name.vector <- gsub(".*tep5", "Step5:Infiltration of immune cells into tumors", name.vector)
  name.vector <- gsub(".*tep6", "Step6:Recognition of cancer cells by T cells", name.vector)
  name.vector <- gsub(".*tep7", "Step7:Killing of cancer cells", name.vector)
  return(name.vector)
 
}