rules_for_tf_fd<-function(dfPatientForAnalysis,se_patient_selection,ge_d,cnv_d,meth_d,input_dataset1_tf,input_dataset2_tf,check_exp=FALSE, check_cnv=FALSE, check_meth=FALSE,check_mut=FALSE,check_exp2=FALSE, check_cnv2=FALSE, check_meth2=FALSE,check_mut2=FALSE,input_GE_tf,input_CNV_tf,input_METH_tf,input_MUTATION_tf,MUT_current_patient)
  {
  
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
    

  if(check_exp == TRUE|check_exp2 == TRUE){
    
      if(check_exp == TRUE){
        
        dataset_colnames<-"input_dataset1_current_patient"
        
      } 
      
      if(check_exp2 == TRUE){ 
        
        dataset_colnames<-"input_dataset2_current_patient"
        
      }
      
    
    zscore<-dfPatientForAnalysis[,"GE_current_patient"]
    
    genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_overexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]>= ge_d)]<-1
    
    genes_downexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_downexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]<= -ge_d)]<-1
    
    genes_lowexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_lowexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]> -ge_d & dfPatientForAnalysis[,"GE_current_patient"]<= -ge_d/2)]<-1
    
    genes_expressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_expressed[which(dfPatientForAnalysis[,"GE_current_patient"]>= ge_d/2 & dfPatientForAnalysis[,"GE_current_patient"]< ge_d)]<-1
    
    genes_otherexp<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_otherexp[which(dfPatientForAnalysis[,"GE_current_patient"]>-ge_d/2 & dfPatientForAnalysis[,"GE_current_patient"]<ge_d/2)]<-1
    
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
    
      
      CNV_TF_gain<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
      
      tcnv<-which(dfPatientForAnalysis[,dataset_colnames] >= 1 & dfPatientForAnalysis[,dataset_colnames] < 2)
      
      if(length(tcnv)==0){
        
        tcnv<-CNV_TF_gain #all genes are not alterated by the copy number variation 
        
      } else {
        
        #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
        CNV_TF_gain[which(dfPatientForAnalysis[,dataset_colnames] >= 1 & dfPatientForAnalysis[,dataset_colnames]< 2)]<-1
        CNV_TF_gain[which(dfPatientForAnalysis[,dataset_colnames] < 1)]<-0
        
      }
      
      ## find amplification
      
      CNV_TF_amplification<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
      
      tcnv<-which(dfPatientForAnalysis[,dataset_colnames] >= 2)
      
      if(length(tcnv)==0){
        
        tcnv<-CNV_TF_amplification #all genes are not alterated by the copy number variation 
        
      } else {
        
        #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
        CNV_TF_amplification[which(dfPatientForAnalysis[,dataset_colnames] >= 2)]<-1
        CNV_TF_amplification[which(dfPatientForAnalysis[,dataset_colnames] < 1)]<-0
        
      }
    }

    ###
    ### Step3: Find Depletion
    ###
    
    CNV_TF_loss<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tcnv<-which(dfPatientForAnalysis[,dataset_colnames] > -2 |  dfPatientForAnalysis[,dataset_colnames]<= -1)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_loss #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_loss[which(dfPatientForAnalysis[,dataset_colnames] > -2 |  dfPatientForAnalysis[,dataset_colnames]<= -1)]<-1
      CNV_TF_loss[which(dfPatientForAnalysis[,dataset_colnames] > -1)]<-0
      
    }
    
    
    #find depletion
    CNV_TF_depletion<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tcnv<-which(dfPatientForAnalysis[,dataset_colnames] <= -2)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_depletion #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_depletion[which(dfPatientForAnalysis[,dataset_colnames] > -2 |  dfPatientForAnalysis[,dataset_colnames]<= -1)]<-1
      CNV_TF_depletion[which(dfPatientForAnalysis[,dataset_colnames] > -1)]<-0
      
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
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_TF_current_patient)]<-1
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] < cnv_TF_current_patient)]<-1
      
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
      
      METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-1
      METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_d)]<-0
      
    }
    
    METH_TF_hypo<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tmeth<-which(dfPatientForAnalysis[,dataset_colnames] >=meth_d |dfPatientForAnalysis[,dataset_colnames] < meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] <  meth_d)]<-1
      METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-0
      
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
      
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_TF_current_patient)]<-1
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_TF_current_patient)]<-1
      
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
                                    genes_overexpressed=genes_overexpressed,
                                    genes_downexpressed=genes_downexpressed,
                                    genes_lowexpressed=genes_lowexpressed,
                                    genes_expressed=genes_expressed,
                                    genes_otherexp=genes_otherexp,
                                    
                                    CNV_TF=rep(cnv_TF_current_patient,nrow(dfPatientForAnalysis)),
                                    CNV_EC_gain=if(cnv_TF_current_patient>=1 & cnv_TF_current_patient<2){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    CNV_EC_amplified=if(cnv_TF_current_patient>=2){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    
                                    CNV_EC_loss=if(cnv_TF_current_patient>=-2 & cnv_TF_current_patient<=-1){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    CNV_EC_depletion=if(cnv_TF_current_patient<=-2){rep(1,nrow(dfPatientForAnalysis))}else{rep(0,nrow(dfPatientForAnalysis))},
                                    
                                    CNV_gain=CNV_TF_gain,
                                    CNV_amplification=CNV_TF_amplification,
                                    CNV_loss=CNV_TF_loss,
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
                     "genes_overexpressed","genes_downexpressed","genes_lowexpressed","genes_expressed","genes_otherexp",
                     "CNV_EC_gain","CNV_EC_amplified","CNV_EC_loss","CNV_EC_depletion",
                     "CNV_gain","CNV_amplification","CNV_loss","CNV_depletion",
                     "CNV_TF_categorization_TF",
                     "METH_EC_hyper","METH_EC_hypo",
                     "METH_hyper","METH_hypo",
                     "METH_TF_categorization_TF","MUT_genes","MUT_TF")
        
    results_rules_TF[[1]]<-dfPatientForAnalysis_GAC
    results_rules_TF[[2]]<-col_relTF
    return(results_rules_TF)
}