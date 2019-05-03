GMIEC_MLK<-function(input_GE_selected,input_CNV_selected,input_METH_selected,input_MUTATION_selected,drugs_for_analysis2,input_clinical,k_user,all_genes_for_analysis){
  
  colGE<-colnames(input_GE_selected)[2:ncol(input_GE_selected)]
  colnames(input_GE_selected)[1]<-"genes"
  colCNV<-colnames(input_CNV_selected)[2:ncol(input_CNV_selected)]
  colnames(input_CNV_selected)[1]<-"genes"
  colMETH<-colnames(input_METH_selected)[2:ncol(input_METH_selected)]
  colnames(input_METH_selected)[1]<-"genes"
  colMUT<-input_MUTATION_selected$Tumor_Sample_Barcode
  
  colnames(drugs_for_analysis2)[1:2]<-c("genes","drugs")   
  
  unique_samples<-colGE
  
  FINAL_DF<-data.frame()
  
  for(se_patient_selection in unique_samples){
    
    if(isTRUE(input_MUTATION_selected$Tumor_Sample_Barcode==se_patient_selection)){
      
      #select the mutation of the current patients                   
      MUT_current_patient_tp<-input_MUTATION_selected[which(input_MUTATION_selected$Tumor_Sample_Barcode==se_patient_selection),]
      #genes Tumor_Sample_Barcode
      #HIST1H2BA      TCGA.A2.A0D2.01
      #CENPF      TCGA.A2.A0ST.01
      #BRCA2      TCGA.A2.A0T0.01
      
      #create a data.frame with 1)name genes of current patient 2) mark mutations as 1 3)  barcode
      MUT_current_patient_temp<-data.frame(genes=MUT_current_patient_tp,mutation=1)   
      
    } else {
      
      print("i do not have a patient with mutations")
      
      print(length(all_genes_for_analysis))
      save.image("C:/Users/guida/Desktop/test.RData")
      dfparse<-matrix(nrow=length(all_genes_for_analysis),ncol=3)
      
      dfparse[,1]<-0 #false gene
      
      dfparse[2:3]<-as.numeric(0)
      
      colnames(dfparse)<-c("genes","mutation","barcode")
      
      dfparse[,"genes"]<-as.character(all_genes_for_analysis)
      dfparse[,"mutation"]<-as.numeric(0)
      
      MUT_current_patient_temp<-data.frame(dfparse[,c(1:2)])
      
    }
    
    
    GE_pzt<-input_GE_selected[,c("genes",se_patient_selection)]
    CNV_pzt<-input_CNV_selected[,c("genes",se_patient_selection)]
    METH_pzt<-input_METH_selected[,c("genes",se_patient_selection)]
    
    genes_common<-Reduce(intersect, list(all_genes_for_analysis,
                                         GE_pzt[,1],
                                         CNV_pzt[,1],
                                         METH_pzt[,1]))
    
    ge_pzt2<-GE_pzt[which(GE_pzt[,1]%in%genes_common),]
    cnv_pzt2<-CNV_pzt[which(CNV_pzt[,1]%in%genes_common),]
    meth_pzt2<-METH_pzt[which(METH_pzt[,1]%in%genes_common),]
    
    
    df_GCM<-data.frame(genes=ge_pzt2[,1],
                       ge=ge_pzt2[,2],
                       cnv=cnv_pzt2[,2],
                       meth=meth_pzt2[,2],
                       stringsAsFactors = F)
    
    INPUT_FOR_RF<-merge(df_GCM,MUT_current_patient_temp,by="genes",stringsAsFactors=F)
    INPUT_FOR_RF$mutation<-as.numeric(as.character(INPUT_FOR_RF$mutation))
    
    rf.pzt.clust <- randomForest(x = INPUT_FOR_RF[,-1], y = NULL, proximity = TRUE, oob.prox = TRUE)
    prox_rf<-rf.pzt.clust$proximity
    
    #users definened
    resKmeans<-kmeans(prox_rf,k_user)$cluster
    RFKR<-data.frame(INPUT_FOR_RF,resKmeans,stringsAsFactors=F)
    colnames(drugs_for_analysis2)[1]<-"genes"
    RFKR_drugs<-merge(RFKR,drugs_for_analysis2,by="genes",all.x=T)
    RFKR_drugs[is.na(RFKR_drugs)]<-0
    
    #group the genes by clusters
    GENES_GROUPED_FOR_K<-aggregate(genes ~ resKmeans, data = unique(RFKR_drugs[,c("genes","resKmeans")]), paste,collapse="#")
    GGFK<-data.frame(t(GENES_GROUPED_FOR_K))[2,]
    colnames(GGFK)<-paste("genes_modules",1:k_user,sep="_")
    #group the drugs by clusters
    DRUGS_GROUPED_FOR_K<-aggregate(drugs ~ resKmeans, data = unique(RFKR_drugs[,c("drugs","resKmeans")]), paste,collapse="#")
    DGFK<-data.frame(t(DRUGS_GROUPED_FOR_K))[2,]
    colnames(DGFK)<-paste("drugs_modules",1:k_user,sep="_")
    
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3978018/
    mu_pos<-function(x){
      vector_x<-x
      x_norm_mu<-mean(abs(vector_x)/max(abs(vector_x)))
      #vector_x_pos<-vector_x[idx_pos]
      #mu_pos<-mean(vector_x_pos)
      return(x_norm_mu)
    }
    
    #count number of amplified genes
    cnv_amp<-function(x){
      vector_x<-x
      idx_cnv_amp<-length(which(vector_x >= 1))
      
      return((idx_cnv_amp/length(vector_x)))
    }
    #count number of deleted genes
    cnv_del<-function(x){
      vector_x<-x
      idx_cnv_del<-length(which(vector_x <= -1))
      
      return((idx_cnv_del/length(vector_x)))
    }
    
    #count drugs score
    countDrugs<-function(sub_df_drugs){
      #genes target by at least one drug
      score_genes_with_drugs<-length(unique(sub_df_drugs[sub_df_drugs[,2]!=0,1]))/length(unique(sub_df_drugs[,1]))
      return(score_genes_with_drugs)
    }
    
    ge_score<-data.frame(t(aggregate(ge ~ resKmeans, data = RFKR,mu_pos)))[2,]
    colnames(ge_score)<-paste("ge_scores",1:k_user,sep="_")
    
    cnv_score_amp<-data.frame(t(aggregate(cnv ~ resKmeans, data = RFKR, cnv_amp)))[2,]
    colnames(cnv_score_amp)<-paste("cnv_score_amp",1:k_user,sep="_")
    
    cnv_score_del<-data.frame(t(aggregate(cnv ~ resKmeans, data = RFKR, cnv_del)))[2,]
    colnames(cnv_score_del)<-paste("cnv_score_dep",1:k_user,sep="_")
    
    meth_score<-data.frame(t(aggregate(meth ~ resKmeans, data = RFKR, mu_pos)))[2,]
    colnames(meth_score)<-paste("meth_score",1:k_user,sep="_")
    
    mut_score<-data.frame(t(aggregate(mutation ~ resKmeans, data = RFKR, FUN=function(x){sum(x)/length(x)})))[2,]
    colnames(mut_score)<-paste("mut_score",1:k_user,sep="_")
    
    tab_SCORES<-cbind(t(ge_score),t(cnv_score_amp),t(cnv_score_del),t(meth_score),t(mut_score))
    colnames(tab_SCORES)<-c("gene_expression","cnv_amp","cnv_del","meth","mutation")
    tab_SCORES[is.na(tab_SCORES)]<-0
    
    GLOBAL_GENETIC_SCORE<-t(data.frame(rowSums(tab_SCORES)))
    colnames(GLOBAL_GENETIC_SCORE)<-paste("gg_score",1:k_user,sep="_")
    
    #create the ouput
    number_genes=t(data.frame(apply(GGFK,2,FUN=function(x){unlist(lapply(strsplit(x,split="#"),length))})))
    number_drugs=t(data.frame(apply(DGFK,2,FUN=function(x){unlist(lapply(strsplit(x,split="#"),length))})))
    
    #count number of genes with drugs
    SCORES_DRUGS<-NULL
    
    for(i in 1:k_user){
      
      res_drug<-countDrugs(RFKR_drugs[RFKR_drugs$resKmeans==i,c("genes","drugs")])
      
      SCORES_DRUGS<-c(SCORES_DRUGS,res_drug)
      
    }
    
    
    SCORES_DRUGS_DF<-t(data.frame(SCORES_DRUGS))
    colnames(SCORES_DRUGS_DF)<-paste("number_of_genes_with_drugs_on_total",1:k_user,sep="_")
    #################
    
    
    #################
    #combine GLOBAL GENETIC SCORE WITH SCORES DRUGS
    
    COMBINED_SCORES=data.frame(GLOBAL_GENETIC_SCORE*SCORES_DRUGS_DF)
    colnames(COMBINED_SCORES)<-paste("combined_score",1:k_user,sep="_")
    #################
    TEMP_DF<-data.frame(sample_id=se_patient_selection,
                        GGFK,
                        number_genes=number_genes,
                        DGFK,
                        number_drugs=number_drugs,
                        ge_score,
                        cnv_score_amp,
                        cnv_score_del,
                        meth_score,
                        mut_score,
                        GLOBAL_GENETIC_SCORE,
                        SCORES_DRUGS_DF,
                        COMBINED_SCORES
    )
    
    FINAL_DF<-rbind(FINAL_DF,TEMP_DF)
    
  }
  colnames(input_clinical)[1]<-"sample_id"
  
FINAL_DF_WITH_CLINIC<-merge(FINAL_DF,input_clinical,by="sample_id",all.x=T)
 
return(FINAL_DF_WITH_CLINIC) 
}