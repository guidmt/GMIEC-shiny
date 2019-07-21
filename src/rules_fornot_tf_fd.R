rules_notfor_tf_fd<-function(dfPatientForAnalysis=dfPatientForAnalysis, se_patient_selection=se_patient_selection, ge_d=ge_d, cnv_d=cnv_d, meth_d=meth_d, MUT_current_patient=MUT_current_patient, check_exp=check_exp, check_cnv=check_cnv,check_meth=check_meth, check_mut=check_mut, check_exp2=check_exp2, check_cnv2=check_cnv2, check_meth2=check_meth2, check_mut2=check_mut2){


  results_rules_TF<-list()
  
  ge_TF_current_patient<-0
  cnv_TF_current_patient<-0
  meth_TF_current_patient <-0
  mutation_TF_current_patient_variant <-0 
    
  #extract the other data of experiment
  TF_ge_rep<-rep(as.numeric(ge_TF_current_patient),nrow(dfPatientForAnalysis))
  TF_CNV_rep<-rep(as.numeric(cnv_TF_current_patient),nrow(dfPatientForAnalysis))
  TF_METH_rep<-rep(as.numeric(meth_TF_current_patient),nrow(dfPatientForAnalysis))
  CNV_TF_categorization<-rep(as.numeric(0),length(dfPatientForAnalysis))
  
  ### Discreterization of the gene expression values
  if(check_exp == TRUE|check_exp2 == TRUE){
    
    if(check_exp == TRUE){
      
      dataset_colnames<-"GE_current_patient"
      
    } 
    
    if(check_exp2 == TRUE){ 
      
      dataset_colnames<-"GE_current_patient"
      
    }
    
    zscore<-dfPatientForAnalysis[,"GE_current_patient"]
    
    genes_overexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_overexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]>= ge_d)]<-1
    
    genes_underexpressed<-rep(0,length(dfPatientForAnalysis[,3]))
    genes_underexpressed[which(dfPatientForAnalysis[,"GE_current_patient"]<= -ge_d)]<-1
  }
  
  ##
  ## Copy number alteration
  ##
  
  if(check_cnv == TRUE|check_cnv2 == TRUE){
    
    
    if(check_cnv == TRUE){
      
      dataset_colnames<-"CNV_current_patient"
      
    } 
    
    if(check_cnv2 == TRUE){ 
      
      dataset_colnames<-"CNV_current_patient"
      
    }
    
    
    CNV_gain<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tcnv<-which(dfPatientForAnalysis[,dataset_colnames] >= 1 & dfPatientForAnalysis[,dataset_colnames])
    
    if(length(tcnv)==0){
      
      tcnv<-CNV_gain #all genes are not alterated by the copy number variation 
      
    } else {
      
      #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
      CNV_gain[which(dfPatientForAnalysis[,dataset_colnames] >= 1)]<-1
      CNV_gain[which(dfPatientForAnalysis[,dataset_colnames] < 1)]<-0
      
    }
    
  ###
  ### Step3: Find Depletion
  ###
  
  #find depletion
  CNV_depletion<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
  
  tcnv<-which(dfPatientForAnalysis[,dataset_colnames] <= -2)
  
  if(length(tcnv)==0){
    
    tcnv<-CNV_depletion #all genes are not alterated by the copy number variation 
    
  } else {
    
    #if is false the condition tcnv ==0, then this mean that all genes have a copy number alteration less than -1 or greater than 1
    CNV_depletion[which(dfPatientForAnalysis[,dataset_colnames] > -2 |  dfPatientForAnalysis[,dataset_colnames]<= -1)]<-1
    CNV_depletion[which(dfPatientForAnalysis[,dataset_colnames] > -1)]<-0
    
  }
  
  

  ##
  ## Step 1.3: categorize the methylation
  ##
  
  if(check_meth == TRUE | check_meth2 == TRUE){
    
    output_meth<-list(1:2)
    
    if(check_meth == TRUE){
      
      dataset_colnames<-"METH_current_patient"
      
    } 
    
    if(check_meth2 == TRUE){ 
      
      dataset_colnames<-"METH_current_patient"
      
    }
    
    
    METH_TF_categorization<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    ###
    ### Step 1.4: categorize the hyper-hypomethylation
    ###
    METH_hyper<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tmeth<-which(dfPatientForAnalysis[,dataset_colnames] >= meth_d |  dfPatientForAnalysis[,dataset_colnames]< meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-1
      METH_hyper[which(dfPatientForAnalysis[,"METH_current_patient"] < meth_d)]<-0
      
    }
    
    METH_hypo<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
    tmeth<-which(dfPatientForAnalysis[,dataset_colnames] >=meth_d |dfPatientForAnalysis[,dataset_colnames] < meth_d)
    
    if(length(tmeth)==0){
      
      tmeth<-METH_hyper #all genes are not alterated by the copy number variation 
      
    } else {
      
      METH_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] <  meth_d)]<-1
      METH_hypo[which(dfPatientForAnalysis[,"METH_current_patient"] >= meth_d)]<-0
      
    }
    
    
  } else {
    
    METH_hypo<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    METH_hyper<-rep(0,length(dfPatientForAnalysis[,dataset_colnames]))
    
  }
  
  ##
  ## Step 1.4: categorize the MUTATION
  ##
  
  if(check_mut == TRUE|check_mut2 == TRUE){
    
    MUT_categorization<-rep(0,nrow(dfPatientForAnalysis))
    
    umutationallpatient<-unique(MUT_current_patient[,"genesID"])
    #find which genes are mutated 
    index_MUT_genes_allPatients<-which(dfPatientForAnalysis[,"genesID"] %in% umutationallpatient)
    MUT_categorization[index_MUT_genes_allPatients]<-1 
  } else {
    MUT_categorization<-rep(0,nrow(dfPatientForAnalysis))
  }
  
  #create a data.frame with the update data
  dfPatientForAnalysis_GAC<-cbind(dfPatientForAnalysis,
                                  
                                  Genes_overexpressed=genes_overexpressed,
                                  Genes_underexpressed=genes_underexpressed,
                                  CNV_gain=CNV_gain,
                                  CNV_depletion=CNV_depletion,
                                  METH_hyper=METH_hyper,
                                  METH_hypo=METH_hypo,
                                  MUT_genes=MUT_categorization)
  
  #columns that describe relation between genes and TF considering other data
  col_relTF<-c("genesID","Genes_overexpressed","Genes_underexpressed",
               "CNV_gain","CNV_depletion","METH_hyper","METH_hypo","MUT_genes")
  
  results_rules_TF[[1]]<-dfPatientForAnalysis_GAC
  results_rules_TF[[2]]<-col_relTF

  }
  return(results_rules_TF)
}