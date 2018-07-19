
engine_all_dataset<-function(input_for_klar2,dfPatientForAnalysis_GAC){
  
  print("Step 6: Run klaR")
  

   resKLAR = kmodes(input_for_klar2[,-1], 5)
   
   ifkl<-cbind(clusters=resKLAR$cluster,input_for_klar2)  
   
   list_rules<-list()
   #Find the rules for each clusters
   unique_clusters<-unique(ifkl[,"clusters"])
   print(unique_clusters)
   for(cksearch in 1:length(unique_clusters)){
     ck<-unique_clusters[cksearch]
     current_rule_value<-apply(ifkl[ifkl[,"clusters"]==ck,-c(1,2)],2,FUN=function(X){sum(as.numeric(X))})#computer number of genes with a properties
     current_rule_names<-names(apply(ifkl[ifkl[,"clusters"]==ck,-c(1,2)],2,FUN=function(X){sum(as.numeric(X))}))
     current_rules_defined<-paste(current_rule_names,current_rule_value,sep=":",collapse=";")
     list_rules[[cksearch]]<-current_rules_defined
   }
   
   ###create a data.frame in with two columns: clusters, rule
   df_clusters_rules<-data.frame(clusters=unique_clusters,rule=unlist(list_rules))
   
   mergeGAC_COM<-merge(ifkl,df_clusters_rules,by="clusters")

   mergeGAC_COM_sort<-mergeGAC_COM[order(mergeGAC_COM$clusters,decreasing=F),]

   rb1<-merge(dfPatientForAnalysis_GAC,mergeGAC_COM_sort) #of default use the same columns
   
   mergeGAC_COM_res_K_2<-rb1

  return(mergeGAC_COM_res_K_2)
}