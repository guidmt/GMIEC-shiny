##
## Plot heatmap function
##

###
### Heatmap summary results
###
plot_heatmap_report_gmiec<-function(input_for_gmiec,type){

      #check if it the ouput of GMIEC MD or GMIEC-FD
      if(type=="genes" & length(grep(colnames(input_for_gmiec),pattern="s_score"))!=0){
        colToselect<-"s_score"
      } 
  
     if(type=="genes" & length(grep(colnames(input_for_gmiec),pattern="s_score"))==0) {
        colToselect<-paste(c("score","genes"),collapse="_")
      }
        
      if(type=="drugs"|type=="sad"){
      colToselect<-paste(c("score",type),collapse="_")
      }
  
      idx<-grep(colnames(input_for_gmiec),pattern=colToselect)
      input_for_gmiec_select<-input_for_gmiec[,idx]
      rownames(input_for_gmiec_select)<-input_for_gmiec[,1]
      colnames(input_for_gmiec_select)<-paste("Module",1:ncol(input_for_gmiec_select))
      
      if(type=="drugs"){
      ht1<-heatmaply(input_for_gmiec_select,colors=PRGn,Colv=FALSE)
      }
      
      if(type=="genes"){
      ht1<-heatmaply(input_for_gmiec_select,colors=BrBG,Colv=FALSE)
      }

      if(type=="sad"){
      ht1<-heatmaply(input_for_gmiec_select,colors=RdYlBu,Colv=FALSE)
      }
      
      return(ht1)
      
}
###
### Table summary modules for patient
###

plotTable<-function(input_for_gmiec,current_patient){
  
  input_for_gmiec_select<-input_for_gmiec[input_for_gmiec[,1]==current_patient,]
  
  idx_number_drugs<-grep(colnames(input_for_gmiec_select),pattern="number_drugs_in_module")
  idx_number_genes<-grep(colnames(input_for_gmiec_select),pattern="number_genes_in_module")
  
  ###
  ### Score
  ###
  
  idx_number_sad<-grep(colnames(input_for_gmiec_select),pattern="sad")
  
  if(length(grep(colnames(input_for_gmiec_select),pattern="s_score"))!=0){
    
    idx_number_sam<-grep(colnames(input_for_gmiec_select),pattern="s_score")
    idx_number_sadrugs<-grep(colnames(input_for_gmiec_select),pattern="score_drugs_rdg")
    
  } else {
    idx_number_sam<-grep(colnames(input_for_gmiec_select),pattern="score_genes")
    idx_number_sadrugs<-grep(colnames(input_for_gmiec_select),pattern="score_drugs_rdg")
  }
  
  
  ##create the data.frame for the statistics about the number of patient for module
  number_drugs_table<-as.numeric(input_for_gmiec_select[,idx_number_drugs])
  number_drugs_table[is.na(number_drugs_table)]<-0
  number_genes_table<-as.numeric(input_for_gmiec_select[,idx_number_genes])
  number_genes_table[is.na(number_genes_table)]<-0
  
  tab_ndrugs_ngenes<-data.frame(Module = paste("Module",seq=1:length(number_drugs_table)),
                                number_of_drugs=as.numeric(number_drugs_table),
                                number_of_genes=as.numeric(number_genes_table),
                                scores_drugs=round(as.numeric(input_for_gmiec_select[,idx_number_sadrugs]),3),
                                scores_genes=round(as.numeric(input_for_gmiec_select[,idx_number_sam]),3),
                                sad=round(as.numeric(input_for_gmiec_select[,idx_number_sad]),3))
  

  
  tab_ndrugs_ngenes %>%
    mutate(
      sad=cell_spec(
        sad, 
        color="white",bold=T,
        background = ifelse(sad >=0, "blue", "red"))
    ,
    scores_genes=cell_spec(
      scores_genes, 
      color="white",bold=T,
      background = ifelse(scores_genes >=0, "cyan", "brown"))
    ,
    scores_drugs=cell_spec(
      scores_drugs, 
      color="white",bold=T,
      background = ifelse(scores_drugs >=0, "green", "violet")),
      number_of_drugs = color_bar("lightgreen")(number_of_drugs),
      number_of_genes = color_bar("lightblue")(number_of_genes)
    ) %>%
    kable(escape = F) %>%
    kable_styling("hover", full_width = F) %>%
    column_spec(c(1), width = "3cm") %>%
    column_spec(c(2,3,4,5,6), width = "3cm") %>%
    row_spec(1:nrow(tab_ndrugs_ngenes),hline_after=F,align='center',color="black")%>%
    row_spec(0,color="black",align='center')
  
}

plot_summary_genes_drugs<-function(input_for_gmiec,current_patient,type,module){
  
  module_to_select<-paste("_",module,sep="")
  
  tab_patient<-input_for_gmiec[input_for_gmiec[,1]==current_patient,]
  
  type_to_select<-paste(type,"in_module",sep="_")
  
  idx_module_select<-grep(grep(grep(colnames(tab_patient),pattern=type_to_select,value=T),pattern="number",invert=T,value=T),pattern=module_to_select,value=T)
  
  df_genes_drugs<-data.frame(tab_patient[,idx_module_select])
  
    select_in_current_module<-df_genes_drugs[,1]
    
    if(type=="genes"){
    
      split_string<-strsplit(as.character(select_in_current_module),split=",")[[1]]
    
      }else{
      
        split_string<-strsplit(as.character(select_in_current_module),split="#")[[1]]
      
      }
    
    if(!is.null(split_string)){
      #i add na to reduce probles in creation matrix
      length(split_string) <- prod(dim(matrix(split_string, ncol = 4)))
      
      if(type=="genes"){
        links_output<-NULL
      #create a vector with the urls gene from Ncbi
      for(current_select in split_string){
        lg<-paste("www.ncbi.nlm.nih.gov/gene/?term=",current_select,sep="")
        lg2<-paste0("https://", lg)
        links_output<-c(links_output,lg2)
      }
      } else{
        
        links_output<-NULL
        split_string[is.na(split_string)]<-0
        #create a vector with the urls gene from Ncbi
        for(current_select in split_string){
          print(current_select)
            if(current_select!=0){
              
          ld<-paste("www.dgidb.org/drugs/",current_select,"#_interactions",sep="")
          ld2<-paste0("http://", ld)
          links_output<-c(links_output,ld2) 
            }
        }
      }
      
    }
    

    dfTable<-data.frame(matrix(split_string, ncol = 1, byrow = TRUE))
    colnames(dfTable)<-"variable"
    
    #finish to create the table with the links
    dfTable %>% mutate(variable=cell_spec(dfTable[,1],"html",link=links_output))%>%
    kable(format="html",escape=F) %>%
      kable_styling(bootstrap_options = "striped", full_width = F)
}

