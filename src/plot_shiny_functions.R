##
## Plot heatmap function
##

###
### Heatmap summary results
###
plot_heatmap_report_gmiec<-function(input_for_gmiec,type){

      #check if it the ouput of GMIEC MD or GMIEC-FD
      if(type=="genes" & length(grep(colnames(input_for_gmiec),pattern="s_score"))!=0){
        type="s_score"
        colToselect<-type
      }else{
        colToselect<-paste(c("score_sad",type),collapse="_")
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
  
  idx_number_drugs<-grep(colnames(input_for_gmiec_select),pattern="drugs_in_module")
  idx_number_genes<-grep(colnames(input_for_gmiec_select),pattern="genes_in_module")
  
  ###
  ### Score
  ###
  
  idx_number_sad<-grep(colnames(input_for_gmiec_select),pattern="sad")
  
  if(length(grep(colnames(input_for_gmiec_select),pattern="s_score"))!=0){
    
    idx_number_sam<-grep(colnames(input_for_gmiec_select),pattern="s_score")
    idx_number_sadrugs<-grep(colnames(input_for_gmiec_select),pattern="rdg")
    
  } else {
    idx_number_sam<-grep(colnames(input_for_gmiec_select),pattern="score_alteration_module")
    idx_number_sadrugs<-grep(colnameS(input_for_gmiec_select),pattern="score_alteration_drugs")
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
      background = ifelse(scores_genes >=0, "cyan2", "brown"))
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

