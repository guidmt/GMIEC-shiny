create_output_gmiec_parse<-function(input_for_gmiec,type_input,type_of_output){
          
          #if the analysis was performed with RFK
          if(type_input=="Random forest + k-means"){
          
          if(type_of_output=="Module active"){
            
            column_to_select<-"s_onc" #get s_onc
            
          } else {
            
            column_to_select<-"s_sup" #get s_sup
            
          }
          
          #otherwise if the user selected to consider the drugs select sad score
          if(type_of_output=="Module active with drugs" | type_of_output=="Module inactive with drugs"){
    
            column_to_select<-"sad"
    
          }
            
          } else {
            
            if(type_of_output=="Module active" | type_of_output=="Module inactive"){
              
              column_to_select<-"score_alteration_module"
              
            }
          
                    # if the type of 
          if(type_of_output=="Module active with drugs" | type_of_output =="Module inactive with drugs"){
    
           column_to_select<-"sad"
    
        }

            
        }
  
  
          #########
          modules_to_select<-grep(colnames(input_for_gmiec),pattern=column_to_select)

          if((type_of_output=="Module active" | type_of_output=="Module inactive") & type_input=="Random forest + k-means"){
          #also for s_sup the algorithm must select the maximum value
          idx_cell_with_score_selected<-apply(input_for_gmiec[,modules_to_select],1,which.max)
          
          } 
        
          if(type_of_output=="Module active" & type_input=="Logic rules + k-means"){
            
            idx_cell_with_score_selected<-apply(input_for_gmiec[,modules_to_select],1,which.max)
            
          } 
          
          if(type_of_output=="Module inactive" & type_input=="Logic rules + k-means"){
            
            idx_cell_with_score_selected<-apply(input_for_gmiec[,modules_to_select],1,which.min)
            
          } 
          
          if(type_of_output=="Module active with drugs"){
            
            idx_cell_with_score_selected<-apply(input_for_gmiec[,modules_to_select],1,which.max)
            
          } 
        
          if(type_of_output=="Module inactive with drugs"){
            
            idx_cell_with_score_selected<-apply(input_for_gmiec[,modules_to_select],1,which.min)
            
          } 
          
          matrix_scores<-input_for_gmiec[,modules_to_select]
          
          idx_genes<-grep(grep(colnames(input_for_gmiec),pattern="genes_in_module",value=T),pattern="number",inver=T,value=T)
          matrix_genes<-input_for_gmiec[,idx_genes]
          
          idx_drugs<-grep(grep(colnames(input_for_gmiec),pattern="drugs_in_module",value=T),pattern="number",inver=T,value=T)
          matrix_drugs<-input_for_gmiec[,idx_drugs]
          
          FINAL_OUTPUT<-data.frame()
          
          for(r in 1:nrow(input_for_gmiec)){
            #the number of rows it is equal to the number of index identified
            
            selected_modules<-idx_cell_with_score_selected[r]
            
            current_patient<-as.character(input_for_gmiec[r,1])
            score_n_module_current_patient<-as.numeric(matrix_scores[r,selected_modules])
            genes_n_module_current_patient<-as.character(matrix_genes[r,selected_modules])
            drugs_n_module_current_patient<-as.character(matrix_drugs[r,selected_modules])
           
            df<-data.frame(patientID=current_patient,
                       score=score_n_module_current_patient,
                       genes_in_module=genes_n_module_current_patient,
                       drugs_in_module=drugs_n_module_current_patient,
                       stringsAsFactors=F) 
            
            FINAL_OUTPUT<-rbind(FINAL_OUTPUT,df)
            
          }

          return(FINAL_OUTPUT)
}