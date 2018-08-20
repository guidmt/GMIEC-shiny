rules_for_tf_fd<-function(dfPatientForAnalysis,se_patient_selection,ge_d,cnv_d,meth_d,input_dataset1_tf,input_dataset2_tf,check_exp=FALSE, check_cnv=FALSE, check_meth=FALSE,check_mut=FALSE,check_exp2=FALSE, check_cnv2=FALSE, check_meth2=FALSE,check_mut2=FALSE,input_GE_tf,input_CNV_tf,input_METH_tf,input_MUTATION_tf,MUT_current_patient) {
  
    results_rules_TF<-list()
  
    print(check_exp)
    print(check_exp2)
    print(check_cnv)
    print(check_cnv2)
    print(check_meth)
    print(check_meth2)
    print(check_mut)
    print(check_mut2)
    
    #check the presence of gene-expression data for TF
    if(check_exp == TRUE | check_exp2 == TRUE ){
      
    ge_TF_current_patient<-as.numeric(input_GE_tf[se_patient_selection]) 
    
    }
  
    #check the presence of cnv data for TF
    if(check_cnv==TRUE| check_cnv2 == TRUE){
    cnv_TF_current_patient<-as.numeric(input_CNV_tf[se_patient_selection])
    cnv_TF_current_patient[is.na(cnv_TF_current_patient)]<-0 #if NA the data is not available
    }
  
    #check the presence of methylation data for TF
    if(check_meth==TRUE|check_meth2==TRUE){
      
    meth_TF_current_patient<-as.numeric(input_METH_tf[se_patient_selection])
    meth_TF_current_patient[is.na(meth_TF_current_patient)]<-0 #if NA the data is not available
    
    } else {
      
      meth_TF_current_patient<-0
    }
  
    #check the presence of mutation data for TF
    if(check_mut==TRUE|check_mut2==TRUE){
      
    if(nrow(input_MUTATION_tf[input_MUTATION_tf$Tumor_Sample_Barcode==se_patient_selection,])){
      
    mutation_TF_current_patient_variant<-0 # no mutation for TF 0
    
    } else {
      
    mutation_TF_current_patient_variant<-1 #yes mutation of TF 1
    
    }
    } else {
      
      mutation_TF_current_patient_variant<-1 #yes mutation of TF 1
      
    }
  

    ##
    ## Step 1.1: categorize the genes associated with the expression or not of the tf using a fold-change values cut-off
    ##
    
    if(check_exp == TRUE|check_exp2 == TRUE){ #if the datasets are from gene-expression
      
      
      if(check_exp == TRUE){ #verify which is the dataset of gene-expression (first or second)
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
     
      if(check_exp2 == TRUE){
        
        dataset_colnames<-"input_dataset2_current_patient"
      }
      
      FC_GE_TF<-abs(dfPatientForAnalysis[,dataset_colnames]-ge_TF_current_patient)
      
      FC_GE_TF_categorization<-rep(0,length(FC_GE_TF))
    
      #where the fold-change is less than to a thresholds then those genes are related with the transcriptional factors
      FC_GE_TF_categorization[which(FC_GE_TF<=ge_d)]<-1 #absolute values
      FC_GE_TF_categorization[which(FC_GE_TF>ge_d)]<-0 
      
    } else { #if are not gene-expression data
        
      FC_GE_TF_categorization<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
      
    }
    

    ### Discreterization of ghe gene expression values according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC151169/
    ### works for that median-centered at the levels of genes and z-score
    
  if(check_exp == TRUE|check_exp2 == TRUE){
    
      if(check_exp == TRUE){
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
      
      if(check_exp2 == TRUE){ 
        
        dataset_colnames<-"input_dataset2_current_patient"
        
      }
      
    genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    names(genes_overexpressed)<-as.character(dfPatientForAnalysis[,1])#give the names to the vector
    rows_to_selected<-round((length(dfPatientForAnalysis[,1])*5)/100)
    #select the over expressed genes in the current_patient
    genes.overexpressed.strings<-dfPatientForAnalysis[order(dfPatientForAnalysis[,dataset_colnames],decreasing=T),][1:rows_to_selected,1]
    genes_overexpressed[which(names(genes_overexpressed) %in% genes.overexpressed.strings)]<-1
    
    genes_underexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    names(genes_underexpressed)<-as.character(dfPatientForAnalysis[,1])#give the names to the vector
    rows_to_selected<-round((length(dfPatientForAnalysis[,1])*5)/100)
    #select the under expressed genes in the current_patient
    genes.underexpressed.strings<-dfPatientForAnalysis[order(dfPatientForAnalysis[,dataset_colnames]),][1:rows_to_selected,1]
    genes_underexpressed[which(names(genes_underexpressed) %in% genes.underexpressed.strings)]<-1
    
    } else {
      genes_underexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
      genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
      
    }
    
    ##
    ## Copy number alteration
    ##
    
    if(check_cnv == TRUE|check_cnv2 == TRUE){
      

    if(check_cnv == TRUE){
      
      dataset_colnames<-"input_dataset1_current_patient"
      
    } 
      
    if(check_cnv2 == TRUE){ 
      
      dataset_colnames<-"input_dataset2_current_patient"
      
    }
    
    CNV_TF_categorization<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    CNV_TF_gain<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tcnv<-which(dfPatientForAnalysis[,dataset_colnames] >= cnv_d |  dfPatientForAnalysis[,dataset_colnames]<= -cnv_d)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_gain #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_gain[which(dfPatientForAnalysis[,dataset_colnames] >= cnv_d)]<-1
      CNV_TF_gain[which(dfPatientForAnalysis[,dataset_colnames] <= -cnv_d)]<-0
      
    }

    ###
    ### Step3: Find Depletion
    ###
    
    CNV_TF_depletion<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tcnv<-which(dfPatientForAnalysis[,dataset_colnames] >= cnv_d |  dfPatientForAnalysis[,dataset_colnames]<= -cnv_d)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_depletion #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_depletion[which(dfPatientForAnalysis[,dataset_colnames] <= -cnv_d)]<-1
      CNV_TF_depletion[which(dfPatientForAnalysis[,dataset_colnames] >= cnv_d)]<-0
      
    }
    
    } else {
  
      CNV_TF_gain<-rep(0,length(dfPatientForAnalysis[,3]))
      CNV_TF_depletion<-rep(0,length(dfPatientForAnalysis[,3]))
      
    }
    
    #now i can test in which case the copy-number is greater or less than the copy number of TFs
    #0 The CNV of the genes is greater than the copy-number variation of TF
    #1 The CNV of the genes is less than the copy-number variation of the TF, interest this case because allow to identify TFs that are alterated
    #and that regulates the genes
    if(check_cnv == TRUE|check_cnv2 == TRUE){ 
      
      if(check_exp == TRUE){
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
      
      if(check_exp2 == TRUE){ 
        
        dataset_colnames<-"input_dataset2_current_patient"
        
      }
      
    if(cnv_TF_current_patient>=cnv_d | cnv_TF_current_patient<=-cnv_d){
      #create a new object to fill
      CNV_TF_categorization_TF<-CNV_TF_categorization
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,dataset_colnames] >= cnv_TF_current_patient)]<-0
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,dataset_colnames] < cnv_TF_current_patient)]<-1
      
    } else {
      
      CNV_TF_categorization_TF<-CNV_TF_categorization
      
    }
      
    } else {
      
      CNV_TF_categorization_TF<-CNV_TF_categorization
    }
    
    ##
    ## Step 1.3: categorize the methylation
    ##
    
if(check_meth == TRUE | check_meth2 == TRUE){
  
      output_meth<-list(1:2)
      
      if(check_meth == TRUE){
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
      
      if(check_meth2 == TRUE){ 
        
        dataset_colnames<-"input_dataset2_current_patient"
        
      }
      
      
    METH_TF_categorization<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    ###
    ### Step 1.4: categorize the hyper-hypomethylation
    ###
    METH_TF_hyper<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tmeth<-which(dfPatientForAnalysis[,dataset_colnames] >= meth_d |  dfPatientForAnalysis[,dataset_colnames]< meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_TF_hyper[which(dfPatientForAnalysis[,dataset_colnames] >= meth_d)]<-1
      METH_TF_hyper[which(dfPatientForAnalysis[,dataset_colnames] < meth_d)]<-0
      
    }
    
    METH_TF_hypo<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tmeth<-which(dfPatientForAnalysis[,dataset_colnames] >=meth_d |dfPatientForAnalysis[,dataset_colnames] < meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_TF_hypo[which(dfPatientForAnalysis[,dataset_colnames] <  meth_d)]<-1
      METH_TF_hypo[which(dfPatientForAnalysis[,dataset_colnames] >= meth_d)]<-0
      
    }
  

    } else {
      
      METH_TF_hypo<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
      METH_TF_hyper<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
      
    }
    
    #now i can test in which case the methylation is greater or less than the methylation of TFs
    #0 The METH of the genes is greater than the methylation variation of TF
    #1 The METH of the genes is less than the methylation variation of the TF, interest this case because allow to identify TFs that are alterated
    #and that regulates the genes
    
    if(check_meth == TRUE | check_meth2 == TRUE){
      
      if(check_meth == TRUE){
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
      
      if(check_meth2 == TRUE){ 
        
        dataset_colnames<-"input_dataset2_current_patient"
        
      }
      
    if(meth_TF_current_patient>=meth_TF_current_patient | meth_TF_current_patient< meth_TF_current_patient){
      
      METH_TF_categorization_TF<-METH_TF_categorization
      
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,dataset_colnames] >= meth_TF_current_patient)]<-0
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,dataset_colnames] < meth_TF_current_patient)]<-1
      
    } else {
      
      
      METH_TF_categorization_TF<-rep(0,length(dfPatientForAnalysis))
      
    }
      
    } else {
      
      
      METH_TF_categorization_TF<-rep(0,length(dfPatientForAnalysis))
      
    }
    
    ##
    ## Step 1.4: categorize the MUTATION
    ##
    
    if(check_mut == TRUE|check_mut2 == TRUE){
      
    MUT_TF_categorization<-rep(0,nrow(dfPatientForAnalysis))

    umutationallpatient<-unique(MUT_current_patient[,"genesID"])
    #find which genes are mutated 
    index_MUT_genes_allPatients<-which(dfPatientForAnalysis[,"genesID"] %in% umutationallpatient)
    MUT_TF_categorization[index_MUT_genes_allPatients]<-1 
    } else {
      MUT_TF_categorization<-rep(0,nrow(dfPatientForAnalysis))
    }
    #create a data.frame with the update data
    dfPatientForAnalysis_GAC<-cbind(dfPatientForAnalysis,
                                    
                                    GexpTF=rep(ge_TF_current_patient,nrow(dfPatientForAnalysis)),
                                    FC_GE_TF=FC_GE_TF_categorization,
                                    Genes_overexpressed=genes_overexpressed,
                                    Genes_underexpressed=genes_underexpressed,
                                    
                                    CNV_TF=rep(cnv_TF_current_patient,nrow(dfPatientForAnalysis)),
                                    CNV_EC_gain=if(cnv_TF_current_patient>=cnv_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    CNV_EC_depletion=if(cnv_TF_current_patient<=-cnv_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    
                                    CNV_gain=CNV_TF_gain,
                                    CNV_depletion=CNV_TF_depletion,
                                    CNV_TF_categorization_TF=CNV_TF_categorization_TF,
                                    
                                    METH_TF=rep(meth_TF_current_patient,nrow(dfPatientForAnalysis)),
                                    METH_EC_hyper=if(meth_TF_current_patient>=meth_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    METH_EC_hypo=if(meth_TF_current_patient<meth_d){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    
                                    METH_hyper=METH_TF_hyper,
                                    METH_hypo=METH_TF_hypo,
                                    METH_TF_categorization_TF=METH_TF_categorization_TF,
                                    
                                    MUT_genes=MUT_TF_categorization,
                                    MUT_TF=rep(mutation_TF_current_patient_variant,nrow(dfPatientForAnalysis)))
    print(dim(dfPatientForAnalysis_GAC))
    #columns that describe relation between genes and TF considering other data
    col_relTF<-c("genesID","FC_GE_TF",
                 "Genes_overexpressed","Genes_underexpressed",
                 "CNV_EC_gain","CNV_EC_depletion",
                 "CNV_gain","CNV_depletion",
                 "CNV_TF_categorization_TF",
                 "METH_EC_hyper","METH_EC_hypo",
                 "METH_hyper","METH_hypo",
                 "METH_TF_categorization_TF","MUT_genes","MUT_TF")
    
    results_rules_TF[[1]]<-dfPatientForAnalysis_GAC
    results_rules_TF[[2]]<-col_relTF
    return(results_rules_TF)
}