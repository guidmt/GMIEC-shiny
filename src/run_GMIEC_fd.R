run_GMIEC_fd<-function(input_dataset1,input_dataset2,check_exp,check_exp2,check_cnv,check_cnv2,check_meth,check_meth2,check_mut,check_mut2,tabDrugs,input_clinical,parameter_discr,genes_annotated_TF_fv,input_GE_tf=NA,input_CNV_tf=NA,input_MUTATION_tf=NA,input_METH_tf=NA,clusters){
  
  withProgress(message="Start analysis!",min=0,max=1,{
    
    Sys.sleep(0.25)
    
    #this is for TCGA data
    parseID<-function(x){
      
      x2<- gsub(x,pattern="-",replacement=".") #replace all - with . 
      x3<- substr(x2,0,15)
      
      return(x3)
      
    }

    #check columns of the genes
    colnames(input_dataset1)[1]<-"genesID"
    colnames(input_dataset2)[1]<-"genesID"

    colnames(tabDrugs)[1]<-"genes"
    #change id clinical data
    colnames(input_clinical)[1]<-"SAMPLE_ID"
    
    
    pts_exp_dataset1<-parseID(colnames(input_dataset1)[-1]) #the first column is the gene-name
    pts_exp_dataset2<-parseID(colnames(input_dataset2)[-1])
    
    print("Do")
    list_pts_experiments<-list(dataset1=pts_exp_dataset1,dataset2=pts_exp_dataset2)
    
    ###
    ### Start the analysis
    ###
    
    print("Step1: Run the analysis")
    
    #ALL_samples_UNIQUE<-common_patient_GE_CNV_METH_MUT
    ALL_samples_UNIQUE<-unique(as.character(unlist(list_pts_experiments)))
    # Increment the progress bar, and update the detail text.
    
    incProgress(0.15, detail = "Step1: Find common patients between all datasets")
    
    print("Step2: Get the data for each patient")
    
    incProgress(0.15, detail = "Step2: Get the data for each patient")
    
    RES_ENGINE<-list()
    
    for(asu in 1:length(ALL_samples_UNIQUE)){
      
      se_patient_selection<-ALL_samples_UNIQUE[asu]
      
      print("Step3: Extract the data from the patient")
      print(se_patient_selection)
      
      input_dataset1_current_patient<-input_dataset1[,c("genesID",se_patient_selection)]
      input_dataset2_current_patient<-input_dataset2[,c("genesID",se_patient_selection)]

      #it is possible that 1 patient does not have a gene
     # if(isTRUE(input_MUTATION_selected$Tumor_Sample_Barcode==se_patient_selection)){
      #  MUT_current_patient<-input_MUTATION_selected[which(input_MUTATION_selected$Tumor_Sample_Barcode==se_patient_selection),]
        
      ##} else {
      #  print("i do not have a patient with mutations")
       # 
       ## dfparse<-matrix(nrow=1,ncol=ncol(input_MUTATION_selected))
      # # dfparse[1]<-"NA" #false gene
       # #dfparse[2:ncol(input_MUTATION_selected)]<-as.numeric(0)
        #colnames(dfparse)<-colnames(input_MUTATION_selected)
        #
      # # MUT_current_patient<-dfparse
      #}
      
      list_DF<-c("input_dataset1_current_patient","input_dataset2_current_patient")
      #check the presence of empty data.frame
      res_nrow<-NULL
      
      for(dfe in list_DF){
        
        current_df<-get(dfe)
        
        res_nrow<-c(res_nrow,dim(current_df)[1])
        
      }
      
      #check2: control which samples does not have genes
      nogenesindataset<-which(res_nrow==0)
      
      if(length(nogenesindataset)==0){
        
        list_DF_clean<-list_DF
        
        DF_notpresent<-"AlldataAvailable"
        
      } else{
        
        list_DF_clean<-list_DF[-nogenesindataset]
        
        DF_notpresent<-list_DF[nogenesindataset]
        
      }
      
      list_df_patients<-list()
      
      for(ldfc in 1:length(list_DF_clean)){
        
        list_df_patients[[ldfc]]<-get(list_DF_clean[ldfc])
        
      }
      
      ###
      ### merge the different experiments for the same patient
      ###
      merge_experiment_patient_df<-Reduce(function(...) merge(...,by="genesID",all=T),list_df_patients)
      
      ##
      ## Start to parse the object merge_experiment_patient_df
      ##
      
      #find the index of columns with the string of patients: Select Data of Copy Number. -> Mutation does not have a column with the name of sample
      index_mepd<-grep(colnames(merge_experiment_patient_df),pattern=se_patient_selection)
      
      #create a new data.frame with hugo symbol and entrez and the values of experiments
      dfPatientForAnalysis<-cbind(merge_experiment_patient_df[,1],merge_experiment_patient_df[index_mepd],DF_notpresent=rep(0,nrow(merge_experiment_patient_df)))
      
      colnames(dfPatientForAnalysis)<-c("genesID",list_DF_clean)
      
      #Change the last column, the experiment without data is always in the last column, see previously line of code
      colnames(dfPatientForAnalysis)[ncol(dfPatientForAnalysis)]<-DF_notpresent
      
      ###
      ### just for test! dfPatientForAnalysis can be used for other analysis not Aprior.
      ###
      print("Step4: ")
      
      incProgress(0.15, detail = "Step3: Find the properties of genes")
      
      #remove NA values 
      dfPatientForAnalysis[is.na(dfPatientForAnalysis)]<-0
      
      rownames(dfPatientForAnalysis)<-paste(dfPatientForAnalysis[,1],seq(1:nrow(dfPatientForAnalysis)))
      
      ##
      ##apply the rules
      ##
      
      #create a data.frame with the update data
      ###
      ### Step1: Identify if the gene expression, cnv, methylation of TFs predicted expression of genes of interest
      ###
      parameter_discr_unlist<-as.numeric(unlist(strsplit(parameter_discr,split=";")))
      ge_d<-parameter_discr_unlist[1]
      cnv_d<-parameter_discr_unlist[2]
      meth_d<-parameter_discr_unlist[3]
      
      if(genes_annotated_TF_fv==TRUE){
        
        print("Step5: Find rules for patient")
        
        
        results_for_tf2<-rules_for_tf_fd(dfPatientForAnalysis=dfPatientForAnalysis,
                                         se_patient_selection=se_patient_selection,
                                         ge_d=ge_d,
                                         cnv_d=cnv_d,
                                         meth_d=meth_d,
                                         check_exp=check_exp,
                                         check_cnv=check_cnv,
                                         check_meth=check_meth,
                                         check_mut=check_mut,
                                         check_exp2=check_exp2,
                                         check_cnv2=check_cnv2,
                                         check_meth2=check_meth2,
                                         check_mut2=check_mut2,
                                         input_GE_tf=input_GE_tf,
                                         input_CNV_tf=input_CNV_tf,
                                         input_METH_tf=input_METH_tf,
                                         input_MUTATION_tf=input_MUTATION_tf,
                                         MUT_current_patient=MUT_current_patient)
        
        dfPatientForAnalysis_GAC<-results_for_tf2[[1]]
        col_relTF<-results_for_tf2[[2]]
        
      } else {
        
        results_for_tf2<-rules_notfor_tf_fd(dfPatientForAnalysis=dfPatientForAnalysis,
                                            se_patient_selection=se_patient_selection,
                                            ge_d=ge_d,
                                            cnv_d=cnv_d,
                                            meth_d=meth_d,
                                            MUT_current_patient=MUT_current_patient,
                                            check_exp=check_exp,
                                            check_cnv=check_cnv,
                                            check_meth=check_meth,
                                            check_mut=check_mut,
                                            check_exp2=check_exp2,
                                            check_cnv2=check_cnv2,
                                            check_meth2=check_meth2,
                                            check_mut2=check_mut2)
        
        dfPatientForAnalysis_GAC<-results_for_tf2[[1]] #this a data.frame with the values of omics experiments and binary rules
        col_relTF<-results_for_tf2[[2]] 
        #col_relTF
        #[1] "genesID"                   "FC_GE_TF"                  "Genes_overexpressed"      
        #[4] "Genes_underexpressed"      "CNV_EC_gain"               "CNV_EC_depletion"         
        #[7] "CNV_gain"      
      }
      
      dfPatientForAnalysis_GAC_rel_TF<-dfPatientForAnalysis_GAC[,col_relTF]
      rownames(dfPatientForAnalysis_GAC_rel_TF)<-dfPatientForAnalysis_GAC_rel_TF[,1]
      #do a control, if the properties of the genes are always equal to 0 it you can remove these columns
      #this step it is important to reduce the computational cost.
      resSumControl<-apply(dfPatientForAnalysis_GAC_rel_TF[,2:ncol(dfPatientForAnalysis_GAC_rel_TF)],2,sum)
      names.good.properties<-names(resSumControl[which(resSumControl!=0)])
      
      input_for_klar<-cbind(genesID=dfPatientForAnalysis_GAC_rel_TF[,1],dfPatientForAnalysis_GAC_rel_TF[,names.good.properties])
      input_for_klar2 <- data.frame(sapply(input_for_klar,as.factor))
      
      ###
      ### run the engine for the analysis 
      ###
      
      print("Step6: Analysis")
      incProgress(0.15, detail = "Step5: Run the engine")
      
      mergeGAC_COM_res_K_2<-engine_all_dataset(input_for_klar2,dfPatientForAnalysis_GAC=dfPatientForAnalysis_GAC,clusters)
      
      print(colnames(mergeGAC_COM_res_K_2))
      #search in which rows are present the genes in the table of drugs-genes interactions
      #two columns in the tables of database must be present: genes and drug_primary_name
      indexSubDrug<-which(tabDrugs$genes%in%mergeGAC_COM_res_K_2[,1])
      subtabDrugs<-tabDrugs[indexSubDrug,]
      
      #i use the symbol "#" to concatenate the strings, because when i will count the number of drugs for gene if are present "," inside the name of
      #drugs i will obtain a mistake number of genes (es. 1-2ethyl,diol,benze)
      collapseDrugTable<-aggregate(drug_primary_name ~ genes, data = subtabDrugs, paste,collapse = "#")
      
      mergeGAC_COM_res_K_2_drugs<-merge(mergeGAC_COM_res_K_2,collapseDrugTable,by.x="genesID",by.y="genes",all.x=T)
      
      countDrugsFunc<-function(x){
        #drugs name are repeated for the same genes, i used unique to manage this issue: the reason is that different database have the same drugs.
        ld<-as.numeric(length(unique(unlist(strsplit(x,split="#")))))
        
        return(ld)
      }
      resCountDrugs<-as.numeric(sapply(mergeGAC_COM_res_K_2_drugs$drug_primary_name,countDrugsFunc))
      resCountDrugs[which(is.na(mergeGAC_COM_res_K_2_drugs$drug_primary_name))]<-0
      
      #test the druggability of a modules
      mergeGAC_COM_res_K_2_drugs<-cbind(mergeGAC_COM_res_K_2_drugs,Count_Drugs_For_Gene=resCountDrugs)
      print(colnames(mergeGAC_COM_res_K_2_drugs))
      
      ### estimate the druggability of the modules
      
      incProgress(0.15, detail = "Step4: Compute the scores")
      
      TOTAL_score_module<-NULL
      TOTAL_score_module_drugs<-NULL
      mergeGAC_COM_res_K_2_drugs$clusters<-as.numeric(mergeGAC_COM_res_K_2_drugs$clusters)
      print(colnames(mergeGAC_COM_res_K_2_drugs))
      
      mergeGAC_COM_res_K_2_drugs<-mergeGAC_COM_res_K_2_drugs[order(mergeGAC_COM_res_K_2_drugs$clusters),]
      print(colnames(mergeGAC_COM_res_K_2_drugs))
      
      uniqGA<-as.numeric(unique(mergeGAC_COM_res_K_2_drugs$clusters))
      
      for(mga in uniqGA){
        
        print(mga)
        
        smgcrkl<- mergeGAC_COM_res_K_2_drugs[mergeGAC_COM_res_K_2_drugs$clusters==mga,]
        
        # $ genesID                  : Factor w/ 22569 levels "1-Dec","1-Mar",..: 22 72 108 152 242 284 293 294 315 368 ...
        # $ FC_GE_TF                 : num  1 1 1 1 1 1 1 1 1 1 ...
        # $ Genes_overexpressed      : num  0 0 0 0 0 0 0 0 0 0 ...
        # $ Genes_underexpressed     : num  0 0 0 0 0 0 0 0 0 0 ...
        
        #check the number of modules alterated: the first column is hugo symbol, the second is the fold-change between the expression (is not an alteration)
        #of tf and target genes. I remove these columns because are not useful in the categorization of alterated and not alterated genes
        
        # col_relTF
        # [1] "genesID"                   "FC_GE_TF"                  "Genes_overexpressed"      
        # [4] "Genes_underexpressed"      "CNV_EC_gain"               "CNV_EC_depletion"         
        # [7] "CNV_gain"                  "CNV_depletion"             "CNV_TF_categorization_TF" 
        # [10] "METH_EC_hyper"             "METH_EC_hypo"              "METH_hyper"               
        # [13] "METH_hypo"                 "METH_TF_categorization_TF" "MUT_genes"                
        # [16] "MUT_TF"
        
        alterated_genes_in_module<-apply(smgcrkl[,col_relTF][,-c(1:2)],1,FUN=function(x){ifelse(sum(x)==0,"Not_altered","Alterated")})
        
        #count the number of not alterated genes in modules
        check_alteration<-cbind(smgcrkl,alterated_genes_in_module)
        
        nrow_module<-nrow(smgcrkl[,col_relTF][,-c(1:2)])
        
        ncol_module<-ncol(smgcrkl[,col_relTF][,-c(1:2)])
        
        #total size module
        total_cell<-nrow_module*ncol_module
        
        number_not_alteration_modules<-length(which(smgcrkl[,col_relTF][,-c(1:2)]==0))
        
        number_alteration_modules<-length(which(smgcrkl[,col_relTF][,-c(1:2)]==1))
        
        #estimate LNAM
        ratio_na_inmodule_with_totsize<-number_not_alteration_modules/total_cell
        #estimate LAM
        ratio_alt_inmodule_with_totsize<-number_alteration_modules/total_cell
        
        #estimate DELTA-A
        
        #estimate the difference between the number of cell without alteration and with alteration in modules
        #positive values indicate that the modules is again integrate otherwise not. [range-1,1]
        
        deltamodulesinalt<-ratio_na_inmodule_with_totsize-ratio_alt_inmodule_with_totsize
        scorestatusmodule<-as.numeric(rep(deltamodulesinalt,nrow_module))
        
        TOTAL_score_module<-c(TOTAL_score_module,scorestatusmodule)
        
        #estimate rdg
        genes_with_drugs<-length(which(smgcrkl[,"Count_Drugs_For_Gene"]>=1))/nrow_module
        
        #estimate NRDG
        genes_without_drugs<-length(which(smgcrkl[,"Count_Drugs_For_Gene"]==0))/nrow_module
        
        #estimate DELTA-D
        deltamodulegenesdrugs<-genes_with_drugs-genes_without_drugs
        deltamodulegenesdrugs <-rep(deltamodulegenesdrugs,nrow_module)
        
        TOTAL_score_module_drugs<-c(TOTAL_score_module_drugs,deltamodulegenesdrugs)
        
      }
      
      #estimate SAD
      mergeGAC_COM_res_K_2_drugs <-cbind(cbind(cbind(mergeGAC_COM_res_K_2_drugs,scorescorestatusmodule=TOTAL_score_module),TOTAL_score_module_drugs),combinedscore=TOTAL_score_module_drugs+TOTAL_score_module)
      
      # > colnames(mergeGAC_COM_res_K_2_drugs)
      # [1] "genesID"                   "FC_GE_TF"                  "Genes_overexpressed"      
      # [4] "Genes_underexpressed"      "CNV_depletion"
      # scorescorestatusmodule
      # TOTAL_score_module_drugs
      # combinedscore
      
      
      RES_ENGINE[[asu]]<-mergeGAC_COM_res_K_2_drugs
      
    }
    names(RES_ENGINE)<-as.character(ALL_samples_UNIQUE)
    
    #res_analysis_each_patient<-grep(ls(),pattern=".analysisGMIEC",value=T,fixed=T)
    
    print("Step5: Save output")
    incProgress(0.15, detail = "Step6: Create the output")
    
    print(length(RES_ENGINE))
    MATRIX_RESULTS_ALL_CLINICAL<-create_output(res_analysis_each_patient=RES_ENGINE,input_clinical=input_clinical)
    #write.table(t(MATRIX_RESULTS_ALL_CLINICAL[-1,]),file="Analysis_GMIEC_main_results.txt",sep="\t",row.names=T,col.names=F,quote=F) # the first row is always empty
    
    print("Step8: send results")
    
    incProgress(0.15, detail = "Step7: Analysis finished, click on download button!")
    
    return(MATRIX_RESULTS_ALL_CLINICAL)
  }) #close the progressbar
}