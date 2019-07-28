create_output<-function(res_analysis_each_patient,input_clinical,...){
  
###extract number of modules in total

GPM_TOT<-NULL


#find the modules in total
for(bpm in 1:length(res_analysis_each_patient)){
  gpm<-unique(res_analysis_each_patient[[bpm]][,"clusters"])
  GPM_TOT<-c(GPM_TOT,gpm)
}

#create the colnames to save the results
patient_id<-"patient_id"
total_genes_patients<-"total_genes_patients"
number_modules<-"number_modules"
total_modules<-unique(GPM_TOT)

number_genes_for_module_size<-paste("#genes_module",total_modules,sep="")
genes_in_each_module<-paste("genes_in_",paste("module",total_modules,sep=""),sep="")
number_drugs_for_module_size<-paste("#drugs_module",total_modules,sep="")
drugs_in_each_module<-paste("drugs_in_",paste("module",total_modules,sep=""),sep="")
scores_module_unique<-paste("score_alteration_genes",total_modules,sep="")
scores_module_drugs_unique<-paste("score_alteration_drugs_rdg",total_modules,sep="")
scores_module_sad_unique<-paste("score_sad",total_modules,sep="")

ALL_colnames<-c(patient_id,total_genes_patients,number_modules
                ,number_genes_for_module_size,genes_in_each_module,
                number_drugs_for_module_size,
                drugs_in_each_module,
                scores_module_unique,
                scores_module_drugs_unique,
                scores_module_sad_unique
)

MATRIX_RESULTS_ALL<-data.frame(matrix(,ncol=length(ALL_colnames)))
colnames(MATRIX_RESULTS_ALL)<-ALL_colnames

print("check bpm")

for(bpm in 1:length(res_analysis_each_patient)){
  print(bpm)
  data_current_patient<-res_analysis_each_patient[[bpm]]
  patient_id<- names(res_analysis_each_patient)[[bpm]]
  names(patient_id)<-"sample_id"
  
  genes_patients<-data_current_patient[,1] 
  total_genes_patients<-length(data_current_patient[,2])
  names(total_genes_patients)<-"total_genes_patients"
  
  number_modules<-length(unique(data_current_patient[,"clusters"]))
  names(number_modules)<-"number_modules"
  
  #extract the number of genes for modules and the genes inside each modules
  number_genes_for_module<-aggregate(genesID ~ clusters, data =  data_current_patient[,c("clusters","genesID")], length)
  
  number_genes_for_module_size<-number_genes_for_module[,2]
  names(number_genes_for_module_size)<-paste("#genes_in_module",number_genes_for_module[,1],sep="")
  
  genes_in_each_module<-aggregate(genesID ~ clusters, data =  data_current_patient[,c("clusters","genesID")], paste,collapse=",")
  genes_in_each_module<-genes_in_each_module[,2]
  names(genes_in_each_module)<-paste("genes_in_",paste("module",number_genes_for_module[,1],sep=""),sep="")
  
  #extract the number of drugs for modules and the genes inside each modules
  number_drugs_for_module<-aggregate(drug_primary_name ~ clusters, data =  data_current_patient[,c("clusters","drug_primary_name")], paste,collapse=",")
  number_drugs_for_module[,2]<-apply(number_drugs_for_module,1,FUN=function(x){length(unlist(strsplit(x[2],split="#")))})
  
  number_drugs_for_module_size<-number_drugs_for_module[,2]
  names(number_drugs_for_module_size)<-paste("#drugs_in_module",number_drugs_for_module[,1],sep="")
  
  drugs_in_each_module<-aggregate(drug_primary_name ~ clusters, data =  data_current_patient[,c("clusters","drug_primary_name")], paste,collapse="@") #use this character to merge drugs because it is important in the next step.
  drugs_in_each_module<-drugs_in_each_module[,2]
  names(drugs_in_each_module)<-paste("drugs_in_",paste("module",number_drugs_for_module[,1],sep=""),sep="")
  
  #extract the alteration scores for each module
  scores_module<-aggregate(scorescorestatusmodule ~ clusters, data =  data_current_patient[,c("clusters","scorescorestatusmodule")], unique)
  scores_module_unique<-scores_module[,2]
  names(scores_module_unique)<-paste("score_alteration_module",scores_module[,1],sep="")
  
  #extract the drugs scores for each module
  scores_module_drugs<-aggregate(TOTAL_score_module_drugs ~ clusters, data =  data_current_patient[,c("clusters","TOTAL_score_module_drugs")], unique)
  scores_module_drugs_unique<-scores_module_drugs[,2]
  names(scores_module_drugs_unique)<-paste("score_alteration_drugs",scores_module_drugs[,1],sep="")
  
  #extract the sad scores for each module
  scores_module_sad<-aggregate(combinedscore ~ clusters, data =  data_current_patient[,c("clusters","combinedscore")], unique)
  scores_module_sad_unique<-scores_module_sad[,2]
  names(scores_module_sad_unique)<-paste("score_sad",scores_module_sad[,1],sep="")
  
  #extract rule for each module 
  rule_in_each_module<-aggregate(clusters ~ rule, data =  data_current_patient[,c("clusters","rule")], unique)
  rule_save<-as.character(rule_in_each_module[,1])
  names(rule_save)<-paste(paste("rule_module",rule_in_each_module[,2],sep=""))
  
  row_for_each_patient<-c(patient_id,
                          total_genes_patients,
                          number_modules,
                          number_genes_for_module_size,
                          genes_in_each_module,
                          number_drugs_for_module_size,
                          drugs_in_each_module,
                          scores_module_unique,
                          scores_module_drugs_unique,
                          scores_module_sad_unique,
                          rule_save
  )
  row_for_each_patient_t<-t(data.frame(row_for_each_patient))
  colnames(row_for_each_patient_t)<-names(row_for_each_patient)
  
  dfrow<-data.frame(row_for_each_patient_t,stringsAsFactors=F)
  colnames(dfrow)<-as.character(names(row_for_each_patient))
  # MATRIX_RESULTS_ALL[[bpm]]
  
  MATRIX_RESULTS_ALL<-rbind.fill(MATRIX_RESULTS_ALL,dfrow)
}

print(dim(MATRIX_RESULTS_ALL))

MATRIX_RESULTS_ALL$sample_id<-gsub(MATRIX_RESULTS_ALL$sample_id,pattern=".analysisGMIEC",replacement="")

MATRIX_RESULTS_ALL_CLINICAL<-merge(MATRIX_RESULTS_ALL,input_clinical,by.x="sample_id",by.y="sample_id")

return(MATRIX_RESULTS_ALL_CLINICAL)

}
