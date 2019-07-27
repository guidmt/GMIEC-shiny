rules_notfor_tf<-function(dfPatientForAnalysis=dfPatientForAnalysis,se_patient_selection=se_patient_selection,ge_d=ge_d,cnv_d=cnv_d,meth_d=meth_d,MUT_current_patient) {
  
  results_rules_TF<-list()
  ge_TF_current_patient<-0
  cnv_TF_current_patient<-0
  meth_TF_current_patient <-0
  mutation_TF_current_patient_variant <-0 
    
  #extract the other data of experiment
  TF_ge_rep<-rep(as.numeric(ge_TF_current_patient),length(dfPatientForAnalysis[,"GE_current_patient"]))
  TF_CNV_rep<-rep(as.numeric(cnv_TF_current_patient),length(dfPatientForAnalysis[,"CNV_current_patient"]))
  TF_METH_rep<-rep(as.numeric(meth_TF_current_patient),length(dfPatientForAnalysis[,"METH_current_patient"]))

  
  ### Discreterization of the gene expression values 
  
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
  
  ##
  ## Step 1.2: categorize the copy-number alteration
  ##
  
  #here i have a different problem, if one gene has a copy number greater than 
  # 
  # #first: test in which genes the copy-number is alterated respect with a copy-number thr
  CNV_TF_categorization<-rep(0,length(dfPatientForAnalysis[,"CNV_current_patient"]))
  
  ###
  ### Step3: Find Gain
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
  
  ##
  ## Step 1.4: categorize the MUTATION
  ##
  
  MUT_TF_categorization<-rep(0,nrow(dfPatientForAnalysis))
  umutationallpatient<-(unique(MUT_current_patient[,1]))
  #find which genes are mutated 
  index_MUT_genes_allPatients<-which(dfPatientForAnalysis[,1] %in% umutationallpatient)
  MUT_TF_categorization[index_MUT_genes_allPatients]<-1
  
  #create a data.frame with the update data
  dfPatientForAnalysis_GAC<-cbind(dfPatientForAnalysis,
                
                                  genes_overexpressed=genes_overexpressed,
                                  genes_downexpressed=genes_downexpressed,

                                  CNV_gain=CNV_TF_gain,
                                  CNV_amplification=CNV_TF_amplification,
                                  
                                  CNV_loss=CNV_TF_loss,
                                  CNV_depletion=CNV_TF_depletion,
                                  
                                  METH_hyper=METH_TF_hyper,
                                  METH_hypo=METH_TF_hypo,
                                  
                                  MUT_genes=MUT_TF_categorization)
  
  #columns that describe relation between genes and TF considering other data
  col_relTF<-c("genesID",
               "genes_overexpressed","genes_downexpressed",
               "CNV_gain","CNV_amplification","CNV_loss","CNV_depletion",
               "METH_hyper","METH_hypo",
               "MUT_genes")
  
  results_rules_TF[[1]]<-dfPatientForAnalysis_GAC
  results_rules_TF[[2]]<-col_relTF
  return(results_rules_TF)
}