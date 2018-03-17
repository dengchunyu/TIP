
#' The function to produce results file of interactive picture for web server
#' 
#' @param mixture_file input data of user data
#' @param sig_matrix input data of signature matrix
#' @param perm =100 by defaut ,the NO. permutations
#' @param QN = TRUE by defaut ,Perform quantile normalization or not
#' @param CHIPorRNASEQ the mixture_file are "Microarray" data or "RNA-seq" data
#' @param saveDir the storage path of results
#' @param sample 'single' and 'multiple' samples to choose
#' @export



CIBERSORT_server<-function(mixture_file,sig_matrix,perm=100, QN=TRUE,CHIPorRNASEQ="RNA-seq",saveDir="",sample="multiple"){
  
  #Two processing path for the choice of parameter 'CHIPorRNASEQ'
  
  if(CHIPorRNASEQ=="Microarray"){
    
    # The result of infiltration proportion.
    
    Results<-CIBERSORT(sig_matrix=sig_matrix,mixture_file=mixture_file,perm=perm,QN=T,saveDir=saveDir,sigMat="LM22",sample=sample)
    
    # Change the sort of cell types, for the Aesthetic quality of radar plot.
    
    cell_type<-colnames(Results)[1:(ncol(Results)-3)]
    
    sort_name_lm22<-c(16,4,21,5,22,10,2,7,11,20,18,17,13,9,1,8,12,15,3,19,14,6)
    
    cell_type<- cell_type[sort_name_lm22]
    
    n1<-paste0('["',rownames(Results)[1],'"')
    
    # Color distribution of 22 immune cell type
    color_lm22<-c("#D03D00","#C3627F","#C84A51","#ACCC69","#969F39","#B3C44C","#7EA33F","#477034","#77B35C","#57B473","#6A84BA","#6067A8","#4BB7A1","#8F88B8","#D2A6CB","#BAB8CC","#5B9EB7","#A7D1E4","#E9D9C3","#C1DDD9","#4FAC6F","#405C70")

    #n m :radar plot; n1 m1 :stackplot
    
    #prepare for radar plot
    Results_radar<-Results[,sort_name_lm22]

        n<-"name"
   
         for(j in 1:22){
      n<- paste(n,cell_type[j],sep = ",")
    }
    
        #Two processing path for the choice of parameter 'sample'
    if(sample=="multiple"){
      
      #radar
      for(i in 1:dim(Results_radar)[1]){
        
        m<-"Sample"
        for(j in 1:22)
        {  
          m<- paste(m,Results_radar[i,j],sep = ",")  
            
        }
        #radar output
        z<-matrix(c(n,m),ncol=1)
        filename<-paste0(saveDir,"/CIBER_",rownames(Results_radar)[i],"_lm22_radar.txt")
        write.table(z,filename,col.names = F,row.names = F,quote = F) 
      }
      
      #stackplot
      for(i in 2:dim(Results)[1]){
        n1<-paste0(n1,',"',rownames(Results)[i],'"')
      }
      n1<-paste0(n1,"]")
      m1_all<-list()
      for(i in 1:22){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm22[i],'","data":[',Results[1,i])
        for(j in 2:dim(Results)[1]){
          m1<-paste(m1,Results[j,i],sep = ",")
        }
        m1<-paste0(m1,"]}")
        m1_all[i]<-m1
      }
      #stackplot output
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)
      filename<-paste0(saveDir,"/CIBER_barplot.txt")
      write.table(z1,filename,col.names = F,row.names = F,quote = F,sep = "\t")
      
      #pie 
      pie_1<-"B cells naive,B cells memory,Plasma cells,T cells CD8,T cells CD4 naive,T cells CD4 memory resting,T cells CD4 memory activated,T cells follicular helper,T cells regulatory (Tregs),T cells gamma delta,NK cells resting,NK cells activated,Monocytes,Macrophages M0,Macrophages M1,Macrophages M2,Dendritic cells resting,Dendritic cells activated,Mast cells resting,Mast cells activated,Eosinophils,Neutrophils"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm22,collapse = ",")
      for(i in 1:dim(Results)[1]){
        pie_2<-paste0(Results[i,1:22],collapse = ",")
        pie_2<-paste0("[",pie_2,"]")
        
        z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)
        filename<-paste0(saveDir,"/CIBER_",rownames(Results)[i],"_lm22_pie.txt")
        write.table(z1,file = filename,col.names = F,row.names = F,quote = F,sep = "\t")
        
      }
     
    }else if(sample=="single"){
      
      #radar one file
      m<-"Sample"
      for(j in 1:22)
      {  
        m<- paste(m,Results_radar[j],sep = ",")  
      }
      z<-matrix(c(n,m),ncol=1)
      filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm22_radar.txt")
      write.table(z,filename,col.names = F,row.names = F,quote = F,sep = "\t")
      
      #stackplot
      n1<-paste0('["',rownames(Results),'"]')
      m1_all<-list()
      for(i in 1:22){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm22[i],'","data":[',Results[1,i],"]}")
        m1_all[i]<-m1
      }
      #stackplot outpur
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)
      filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm22_barplot.txt")
      write.table(z1,filename,col.names = F,row.names = F,quote = F,sep = "\t")
      
      #pie
      pie_1<-"B cells naive,B cells memory,Plasma cells,T cells CD8,T cells CD4 naive,T cells CD4 memory resting,T cells CD4 memory activated,T cells follicular helper,T cells regulatory (Tregs),T cells gamma delta,NK cells resting,NK cells activated,Monocytes,Macrophages M0,Macrophages M1,Macrophages M2,Dendritic cells resting,Dendritic cells activated,Mast cells resting,Mast cells activated,Eosinophils,Neutrophils"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm22,collapse = ",")
      pie_2<-paste0(Results[1:22],collapse = ",")
      pie_2<-paste0("[",pie_2,"]")
      
      z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)
      filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm22_pie.txt")
      write.table(z1,file = filename,col.names = F,row.names = F,quote = F,sep = "\t")
    }
    
  }else if(CHIPorRNASEQ=="RNA-seq"){
      
    sort_name_lm14<-c(13,3,7,4,14,9,10,6,5,1,11,12,2,8)
    Results<-CIBERSORT(sig_matrix=sig_matrix,mixture_file=mixture_file,perm=perm,QN=F,saveDir=saveDir,sigMat="LM14",sample=sample)
    
    cell_type<-colnames(Results)[1:(ncol(Results)-3)]
    
        n1<-paste0('["',rownames(Results)[1],'"')
        color_lm14<-c("#928DBE","#77B35C","#E7D429","#F49E14","#A7D2E6","#5B9EB7","#504A8C","#E9D9C3","#EE8887","#CE482B","#E3C047","#3AB9DD","#4FAC6F","#405C70")
        

   
        Results_radar<-Results[,sort_name_lm14]
        cell_type<- cell_type[sort_name_lm14]
        
        n<-"name"
        for(j in 1:14){
          n<- paste(n,cell_type[j],sep = ",")
        }

        
        #Two processing path for parameter 'sample'
        if(sample=="multiple"){
          #radar
          for(i in 1:dim(Results_radar)[1]){
            m<-"Sample"
          
            for(j in 1:14)
            {  
              m<- paste(m,Results_radar[i,j],sep = ",") 
              
            }
            #radar output
            z<-matrix(c(n,m),ncol=1)
            filename<-paste0(saveDir,"/CIBER_",rownames(Results_radar)[i],"_lm14_radar.txt")
            write.table(z,filename,col.names = F,row.names = F,quote = F) 
            }

          #stackplot
          for(i in 2:dim(Results)[1]){
            n1<-paste0(n1,',"',rownames(Results)[i],'"')
          }
          n1<-paste0(n1,"]")
          m1_all<-list()
          for(i in 1:14){
            m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm14[i],'","data":[',Results[1,i])
            for(j in 2:dim(Results)[1]){
              m1<-paste(m1,Results[j,i],sep = ",")
            }
            m1<-paste0(m1,"]}")
            m1_all[i]<-m1
          }
          #stackplot output
          m1_all<-paste0(unlist(m1_all),collapse=",")
          m1_all<-paste0("[",m1_all,"]")
          z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)
          filename<-paste0(saveDir,"/CIBER_barplot.txt")
          write.table(z1,filename,col.names = F,row.names = F,quote = F,sep = "\t")
          
          #pie
          pie_1<-"B cells,CD4 Naive,CD4 Memory,CD8 Naive,CD8 Memory,CD8 Effector,Treg cell,Th cell,Monocytes CD16,Monocytes CD14,DC,pDC,NK,Plasma"
          pie_1<-paste0("[",pie_1,"]")
          pie_3<-paste0(color_lm14,collapse = ",")
          for(i in 1:dim(Results)[1]){
            pie_2<-paste0(Results[i,1:14],collapse = ",")
            pie_2<-paste0("[",pie_2,"]")
            
            z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)
            filename<-paste0(saveDir,"/CIBER_",rownames(Results)[i],"_lm14_pie.txt")
            write.table(z1,file = filename,col.names = F,row.names = F,quote = F,sep = "\t")
            
          }
        
      }else if(sample=="single"){
        
       #radar
        m<-"Sample"
    
        for(j in 1:14)
        {  
          m<- paste(m,Results_radar[j],sep = ",")  
         
        }
        z<-matrix(c(n,m),ncol=1)
        filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm14_radar.txt")
        write.table(z,filename,col.names = F,row.names = F,quote = F,sep = "\t")
        
        #stackplot
          n1<-paste0('["',rownames(Results),'"]')
        m1_all<-list()
        for(i in 1:14){
          m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm14[i],'","data":[',Results[1,i],"]}")
          m1_all[i]<-m1
        }
        #stackplot output
        m1_all<-paste0(unlist(m1_all),collapse=",")
        m1_all<-paste0("[",m1_all,"]")
        z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)
        filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm14_barplot.txt")
        write.table(z1,filename,col.names = F,row.names = F,quote = F,sep = "\t") 
        
        #pie
        pie_1<-"B cells,CD4 Naive,CD4 Memory,CD8 Naive,CD8 Memory,CD8 Effector,Treg cell,Th cell,Monocytes CD16,Monocytes CD14,DC,pDC,NK,Plasma"
        pie_1<-paste0("[",pie_1,"]")
        pie_3<-paste0(color_lm14,collapse = ",")
          pie_2<-paste0(Results[1:14],collapse = ",")
          pie_2<-paste0("[",pie_2,"]")
          
          z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)
          filename<-paste0(saveDir,"/CIBER_",rownames(Results),"_lm14_pie.txt")
          write.table(z1,file = filename,col.names = F,row.names = F,quote = F,sep = "\t")
        
      }
   }
}

  
  
