filter_gene_expression<-function(input_GE,genes_for_analysis,mad_thr){

  input_GE_selectedgenes<-input_GE[input_GE[,1]%in%genes_for_analysis,]
  
  #remove possible NA values
  input_GE_selectedgenes[is.na(input_GE_selectedgenes)]<-0
  
  
  ### The genes derived in the first step are first selected by expression value. Filter out the genes with low expression
  
  #estimate variance of selected genes
  sdevestimation<-apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,sd)
  #estimate mean of selected genes
  mean_genes<-apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,mean)
  
  #create a new.data.frame with the coefficient of variation
  input_GE_selectedgenes<-data.frame(input_GE_selectedgenes,mean=mean_genes,mad =apply(input_GE_selectedgenes[,3:ncol(input_GE_selectedgenes)],1,mad))
  #remove na
  input_GE_selectedgenes_filter<-input_GE_selectedgenes[which(!is.na(input_GE_selectedgenes$mad)),]
  #sort by cv
  input_GE_selectedgenes2<-input_GE_selectedgenes_filter[order(input_GE_selectedgenes_filter$mad,decreasing=TRUE),]
  
  if(isTRUE(mad_behavior=="AUTO")){
    
    #find the flex point of the distribution of CV
    bcp_results<-bcp(input_GE_selectedgenes2$mad,mcmc = 1000)
    
    #find where is the minimum break points (this correspond with the flex point)
    min_bp<-min(bcp_results$blocks)
    #find the corresponding values of c.v. in min_bp
    min_cv_automatic_thresholds<-input_GE_selectedgenes2$mad[min_bp]
    
    #i selected all genes with a thresholds less than or greater to the min cv defined after break points changes identification
    input_GE_selectedgenes2_mad_select_afterBCP<-input_GE_selectedgenes2[which(input_GE_selectedgenes2$mad>=min_cv_automatic_thresholds),]
    
  } else {
    
    input_GE_selectedgenes2_mad_select_afterBCP<-input_GE_selectedgenes2[which(input_GE_selectedgenes2$mad>=as.numeric(mad_thr)),]
    
  }
  
  
  input_GE_afterCVBCP_filtering<-input_GE_selectedgenes2_mad_select_afterBCP
  
  check_ge_for_patients<-input_GE_afterCVBCP_filtering

return(check_ge_for_patients)

}