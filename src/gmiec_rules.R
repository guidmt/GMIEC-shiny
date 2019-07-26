GMIEC_RULES<-function(input_GE_selected,input_CNV_selected,input_METH_selected,input_MUTATION_selected,drugs_for_analysis2,input_clinical,k_user,all_genes_for_analysis){
  
    results_rules_TF<-list()
     
    #IF GMIEC wants analyze TF use the data of TF
    if(analysis_tf==TRUE){
      
    ge_TF_current_patient<-as.numeric(input_GE_tf[se_patient_selection]) 
    cnv_TF_current_patient<-as.numeric(input_CNV_tf[se_patient_selection])
    cnv_TF_current_patient[is.na(cnv_TF_current_patient)]<-0 #if NA the data is not available
    meth_TF_current_patient<-as.numeric(input_METH_tf[se_patient_selection])
    meth_TF_current_patient[is.na(meth_TF_current_patient)]<-0 #if NA the data is not available
    
    if(nrow(input_MUTATION_tf[input_MUTATION_tf$Tumor_Sample_Barcode==se_patient_selection,])){
      
    mutation_TF_current_patient_variant<-0 # no mutation for TF 0
    
    } else {
      
    mutation_TF_current_patient_variant<-1 #yes mutation of TF 1
    
    }
    
    } else {    #IF GMIEC no analyze TF create empty vector

      
      ge_TF_current_patient<-0
      cnv_TF_current_patient<-0
      meth_TF_current_patient <-0
      mutation_TF_current_patient_variant <-0 
      
    }
    
    
    
    #extract the other data of experiment
    TF_ge_rep<-rep(as.numeric(ge_TF_current_patient),length(dfPatientForAnalysis[,"GE_current_patient"]))
    TF_CNV_rep<-rep(as.numeric(cnv_TF_current_patient),length(dfPatientForAnalysis[,"CNV_current_patient"]))
    TF_METH_rep<-rep(as.numeric(meth_TF_current_patient),length(dfPatientForAnalysis[,"METH_current_patient"]))
      

    ##
    ## Step 1.1: categorize the genes associated with the expression or not of the tf using a fold-change values cut-off
    ##
    FC_GE_TF<-abs(dfPatientForAnalysis[,"GE_current_patient"]-ge_TF_current_patient)
    FC_GE_TF_categorization<-rep(0,length(FC_GE_TF))
    
    #where the fold-change is less than to a thresholds then those genes are related with the transcriptional factors
    FC_GE_TF_categorization[which(FC_GE_TF<=ge_d)]<-1 #absolute values
    FC_GE_TF_categorization[which(FC_GE_TF>ge_d)]<-0 
    
    
    ### Discreterization of ghe gene expression values 
    
    zscore<-dfPatientForAnalysis[,"GE_current_patient"]
    
    genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_overexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]>= ge_d)]<-1
    
    genes_downexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_downexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]<= -ge_d)]<-1

    ##
    ## Step 1.2: categorize the copy-number alteration
    ##
    
    #here i have a different problem, if one gene has a copy number greater than 
    # 
    # #first: test in which genes the copy-number is alterated respect with a copy-number thr
    CNV_TF_categorization<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
    
    ###
    ### Step3: Find Gain-amplication
    ###
    CNV_TF_gain<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
    
    tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] >= 1 & dfPatientForAnalysis[,"CNV_current_patient"] < 2)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_gain #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_gain[which(dfPatientForAnalysis[,"CNV_current_patient"] >= 1 & dfPatientForAnalysis[,"CNV_current_patient"]< 2)]<-1
      CNV_TF_gain[which(dfPatientForAnalysis[,"CNV_current_patient"] < 1)]<-0
      
    }
    
    ## find amplification
    
    CNV_TF_amplification<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
    
    tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] >= 2)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_amplification #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_amplification[which(dfPatientForAnalysis[,"CNV_current_patient"] >= 2)]<-1
      CNV_TF_amplification[which(dfPatientForAnalysis[,"CNV_current_patient"] < 1)]<-0
      
    }
    
    ###
    ### Step3: Find loss
    ###
    
    CNV_TF_loss<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
    
    tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] > -2 |  dfPatientForAnalysis[,"CNV_current_patient"]<= -1)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_loss #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_loss[which(dfPatientForAnalysis[,"CNV_current_patient"] > -2 |  dfPatientForAnalysis[,"CNV_current_patient"]<= -1)]<-1
      CNV_TF_loss[which(dfPatientForAnalysis[,"CNV_current_patient"] > -1)]<-0
      
    }
    
    
    #find depletion
    CNV_TF_depletion<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
    
    tcnv<-which(dfPatientForAnalysis[,"CNV_current_patient"] <= -2)
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_TF_depletion #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_TF_depletion[which(dfPatientForAnalysis[,"CNV_current_patient"] > -2 |  dfPatientForAnalysis[,"CNV_current_patient"]<= -1)]<-1
      CNV_TF_depletion[which(dfPatientForAnalysis[,"CNV_current_patient"] > -1)]<-0
      
    }
    
    #now i can test in which case the copy-number is greater or less than the copy number of TFs
    #0 The CNV of the genes is greater than the copy-number variation of TF
    #1 The CNV of the genes is less than the copy-number variation of the TF, interest this case because allow to identify TFs that are alterated
    #and that regulates the genes
    
    
    if(cnv_TF_current_patient>=cnv_d | cnv_TF_current_patient<=-cnv_d){
      #create a new object to fill
      CNV_TF_categorization_TF<-CNV_TF_categorization
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] >= cnv_TF_current_patient)]<-1
      CNV_TF_categorization_TF[which(dfPatientForAnalysis[,"CNV_current_patient"] < cnv_TF_current_patient)]<-1
      
    } else {
      
      CNV_TF_categorization_TF<-CNV_TF_categorization
      
    }    
    # ##
    # ## Step 1.3: categorize the methylation
    # ##
    # 
    METH_TF_categorization<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))
    
    ###
    ### Step 1.4: categorize the hyper-hypomethylation
    ###
    
    METH_TF_hyper<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))
    
    tmeth<-which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d |  dfPatientForAnalysis[,"METH_current_patient"]< meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-1
      METH_TF_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_d)]<-0
      
    }
    
    METH_TF_hypo<-rep(0,length(dfPatientForAnalysis[,"METH_current_patient"]))
    
    tmeth<-which(dfPatientForAnalysis[,"METH_current_patient"] >=meth_d |dfPatientForAnalysis[,"METH_current_patient"] < meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_TF_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] <  meth_d)]<-1
      METH_TF_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-0
      
    }
    
    #now i can test in which case the methylation is greater or less than the methylation of TFs
    #0 The METH of the genes is greater than the methylation variation of TF
    #1 The METH of the genes is less than the methylation variation of the TF, interest this case because allow to identify TFs that are alterated
    #and that regulates the genes
    
    if(meth_TF_current_patient>=meth_TF_current_patient | meth_TF_current_patient< meth_TF_current_patient){
      
      METH_TF_categorization_TF<-METH_TF_categorization
      
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_TF_current_patient)]<-1
      METH_TF_categorization_TF[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_TF_current_patient)]<-1
      
    } else {
      
      
      METH_TF_categorization_TF<-METH_TF_categorization
      
    }
    ##
    ## Step 1.4: categorize the MUTATION
    ##
    
    MUT_TF_categorization<-rep(0,nrow(dfPatientForAnalysis))
    print("test")
    print(MUT_current_patient)
    umutationallpatient<-unique(MUT_current_patient[,"genesID"])
    #find which genes are mutated 
    index_MUT_genes_allPatients<-which(dfPatientForAnalysis[,"genesID"] %in% umutationallpatient)
    MUT_TF_categorization[index_MUT_genes_allPatients]<-1 
  
    #create a data.frame with the update data
    dfPatientForAnalysis_GAC<-cbind(dfPatientForAnalysis,
                                    
                                    GexpTF=rep(ge_TF_current_patient,nrow(dfPatientForAnalysis)),
                                    FC_GE_TF=FC_GE_TF_categorization,
                                    genes_overexpressed=genes_overexpressed,
                                    genes_downexpressed=genes_downexpressed,
                       
                                    
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
    
    #columns that describe relation between genes and TF considering other data
    col_relTF<-c("genesID","FC_GE_TF",
                 "genes_overexpressed","genes_downexpressed",
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