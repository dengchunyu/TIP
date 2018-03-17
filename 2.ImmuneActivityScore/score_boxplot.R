Score_boxplot<-function(Exp.class.matrix,saveDir){
  
  step_num <- unique(Exp.class.matrix[,2])
  stage_name <- c("Step1:Release of cancer antigens","Step2:Cancer antigen presentation",
                  "Step3:Priming and activation","Step4:Trafficking of immune cells to tumors",
                  "Step5:Infiltration of immune cells into tumors","Step6:Recognition of cancer cells by T cells",
                  "Step7:Killing of cancer cells")
  stage_name <- stage_name[step_num]
  
  for(j in 3:ncol(Exp.class.matrix)){
    result_1<-Exp.class.matrix[,c(1,2,3)]
    name_score<-unlist(apply(result_1,1,function(x){
    R<-round(runif(1,min=0.1,max=0.9),1)
    paste0('{"name":"',x[1],'","x":',R,',"y":',x[3],'}')
  }))
  
    result_1_merge<-tapply(name_score,result_1[,2],function(y){paste0(y,collapse=",")})
    m1_all<-list()
    for(i in 1:length(step_num)){
      m1<-paste0('{"stageName":"',stage_name[i],'","value":[',result_1_merge[i],"]}")
      m1_all[i]<-m1
    }
    
    m1_all<-paste0(unlist(m1_all),collapse=",")
    m1_all<-paste0("[",m1_all,"]")
    filename<-paste0(saveDir,"/",colnames(Exp.class.matrix)[j],"_score_boxplot.txt")
    write.table(m1_all,filename,col.names = F,row.names = F,quote = F,sep = "\t")       
  }
}
