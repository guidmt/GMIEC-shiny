internal_annotation<-function(input_TF_fr,input_annotation,distance_max){
  
  
  #remove the chr string from the columns of chromosome, only chromosome wihtout chr are allowed
  input_TF_fr[,1]<-gsub(input_TF_fr[,1],pattern="chr",replacement="")
  colnames(input_TF_fr)<-c("seqnames","start","end")
  peaks <- toGRanges(input_TF_fr[,c(1:3)], format="broadPeak") #of default use only three columns
  
  colnames(input_annotation)<-c("seqnames","start","end","gene_name")
  input_annotation[,1]<-gsub(input_annotation[,1],pattern="chr",replacement="")
  
  rownames(input_annotation)<-paste(input_annotation[,4],seq(1,nrow(input_annotation)),sep="_")
  
  annotationGR<-makeGRangesFromDataFrame(input_annotation)
  
  anno <- annotatePeakInBatch(peaks, AnnotationData=annotationGR,output="nearestLocation",
                              maxgap=distance_max, ignore.strand=TRUE) #time consuming with huge number of peaks
  
  annotation_results_current_dist<-as.data.frame(anno)
  
  subsetAnnoDF<-annotation_results_current_dist[,c("feature","distancetoFeature")]
  subsetAnnoDF[,1]<-sapply(strsplit((subsetAnnoDF$feature),split="_"),"[[",1)
  
  #first column genes #second column distance
  subsetAnnoDF_unique<-aggregate(distancetoFeature ~ feature, data=subsetAnnoDF, FUN=mean)[,1]#return only the gene-symbol
  return(subsetAnnoDF_unique)
}