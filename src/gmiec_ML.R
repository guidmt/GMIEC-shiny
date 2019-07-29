GMIEC_MLK<-function(input_GE_selected,input_CNV_selected,input_METH_selected,input_MUTATION_selected,drugs_for_analysis2,input_clinical,k_user,all_genes_for_analysis){
  
  colGE<-colnames(input_GE_selected)[2:ncol(input_GE_selected)]
  #compute the z-score between the patients
  input_GE_selected[,-1]<-log(input_GE_selected[,-1]+1,2)
  
  input_GE_selected[,-1]<-apply(X=input_GE_selected[,-1],1,FUN=function(X){(X-mean(X))/sd(X)})
  
  colnames(input_GE_selected)[1]<-"genes"
  colCNV<-colnames(input_CNV_selected)[2:ncol(input_CNV_selected)]
  colnames(input_CNV_selected)[1]<-"genes"
  colMETH<-colnames(input_METH_selected)[2:ncol(input_METH_selected)]
  colnames(input_METH_selected)[1]<-"genes"
  #the column of the samples with mutations must be only the second
  colMUT<-input_MUTATION_selected[,2]
  
  colnames(drugs_for_analysis2)[1:2]<-c("genes","drugs")   
  
  #find all common patients
  unique_samples<-Reduce(intersect, list(colGE[-1],colCNV[-1],colMETH[-1],colMUT))

  FINAL_DF<-data.frame()
  
  for(se_patient_selection in unique_samples){
    
    if(length(grep(input_MUTATION_selected[,2],pattern=se_patient_selection))!=0){
      
      #select the mutation of the current patients                   
      MUT_current_patient_tp<-input_MUTATION_selected[which(input_MUTATION_selected[,2]==se_patient_selection),]
      #genes Tumor_Sample_Barcode
      #HIST1H2BA      TCGA.A2.A0D2.01
      #CENPF      TCGA.A2.A0ST.01
      #BRCA2      TCGA.A2.A0T0.01
      
      #create a data.frame with 1)name genes of current patient 2) mark mutations as 1 3)  barcode
      MUT_current_patient_temp<-data.frame(genes=MUT_current_patient_tp,mutation=1)   
      
    } else {
      
      print("i do not have a patient with mutations")
      
      print(length(all_genes_for_analysis))

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
    
    zscore_ge<-GE_pzt 
    
    #genes_common<-all_genes_for_analysis deprecated, add step to check overlap genes
    #mutation data are add next with a match
    genes_common<-Reduce(intersect,list(GE_pzt[,1],CNV_pzt[,1],METH_pzt[,1]))
    
    ge_pzt2<-GE_pzt[which(GE_pzt[,1]%in%genes_common),]
    cnv_pzt2<-CNV_pzt[which(CNV_pzt[,1]%in%genes_common),]
    meth_pzt2<-METH_pzt[which(METH_pzt[,1]%in%genes_common),]
    meth_pzt2[is.na(meth_pzt2)]<-0
    
    df_temp<-data.frame(matrix(0,ncol=5,nrow=length(genes_common)))
    df_temp[,1]<-genes_common
    
    #fill the data frame with the omics data
    df_temp[match(ge_pzt2[,1],df_temp[,1]),2]<-ge_pzt2[,2]
    df_temp[match(cnv_pzt2[,1],df_temp[,1]),3]<-cnv_pzt2[,2]
    df_temp[match(meth_pzt2[,1],df_temp[,1]),4]<-meth_pzt2[,2]
    
    #check if for the common genes exists a mutation
    if(!is.na(match(MUT_current_patient_temp[,1],df_temp[,1]))){
      
    #select only the genes with mutation in common with the other datatasets
    select_gene_mut<-which(MUT_current_patient_temp[,1] %in% df_temp[,1])
      
    df_temp[match(MUT_current_patient_temp[select_gene_mut,1],df_temp[,1]),5]<-1 #if there is a mutation is always 1
    
    } #else the data.fram contains already 0 values
    
    colnames(df_temp)<-c("genes","ge","cnv","meth","mutation")
    
    INPUT_FOR_RF<-df_temp
    INPUT_FOR_RF[is.na(INPUT_FOR_RF)]<-0
  
    #run randomForest
    rf.pzt.clust <- randomForest(x = INPUT_FOR_RF[,-1], y = NULL, proximity = TRUE, oob.prox = TRUE)
    prox_rf<-rf.pzt.clust$proximity
    
    #k-means using the clusters defined by the user
    resKmeans<-kmeans(prox_rf,k_user)$cluster
    RFKR<-data.frame(INPUT_FOR_RF,resKmeans,stringsAsFactors=F)
    colnames(drugs_for_analysis2)[1]<-"genes"
    
    #merge with the data.frame genes-drugs
    RFKR_drugs<-merge(RFKR,drugs_for_analysis2,by="genes",all.x=T)
    RFKR_drugs[is.na(RFKR_drugs)]<-0
            
    #find hyper-methylated genes
    genes_ismix<-meth_pzt2[which(meth_pzt2[,2]>=0.7),1]
    
    #create a column to save which are the genes hyper-methylated
    RFKR<-data.frame(RFKR,meth_expressed=0)  
    RFKR[match(genes_ismix,RFKR$gene),"meth_expressed"]<-1

    #find up-regulated and down-regulated genes
    zscore<-ge_pzt2[,2]
    
    idx_up_genes<-which(zscore >= 2)
    idx_down_genes<-which(zscore <= -2)
    
    genes_up<-ge_pzt2[idx_up_genes,1]
    genes_down<-ge_pzt2[idx_down_genes,1]
    
    #create a column to save which are the genes over-expressed
    RFKR<-data.frame(RFKR,genes_up=0)  
    RFKR[match(genes_up,RFKR$genes),"genes_up"]<-1
    
    #create a column to save which are the genes under-expressed
    RFKR<-data.frame(RFKR,genes_down=0)  
    RFKR[match(genes_down,RFKR$genes),"genes_down"]<-1
    
    #group the genes by clusters
    GENES_GROUPED_FOR_K<-aggregate(genes ~ resKmeans, data = unique(RFKR_drugs[,c("genes","resKmeans")]), paste,collapse=",")
    # GGFK contains for each module the genes
    GGFK<-data.frame(t(GENES_GROUPED_FOR_K))[2,]
    colnames(GGFK)<-paste("genes_in_module",1:k_user,sep="_")
    
    #group the drugs by clusters
    DRUGS_GROUPED_FOR_K<-aggregate(drugs ~ resKmeans, data = unique(RFKR_drugs[,c("drugs","resKmeans")]), paste,collapse="#")
    DGFK<-data.frame(t(DRUGS_GROUPED_FOR_K))[2,]
    colnames(DGFK)<-paste("drugs_in_module",1:k_user,sep="_")
    
    #
    # Computation S-SCORE, definition of the functions
    #

    # for S-score: gene-expression, amplication, up-regulated genes, down-regulated genes 
    # count the number of 1 in one vector
    
    s_onc<-function(x){
      vector_x<-x
      idx_s_onc<-length(which(vector_x >= 1))
      return((idx_s_onc/length(vector_x)))
    }
    
    #count number of deleted genes
    cnv_del<-function(x){
      vector_x<-x
      idx_cnv_del<-length(which(vector_x <= -1))
      return((idx_cnv_del/length(vector_x)))
    }
    
    #count the number of mutations
    mut_count<-function(x){
      vector_x<-x
      idx_mut<-length(which(vector_x%in%1))
      return((idx_mut/length(vector_x)))
    }
    
    #count drugs score
    countDrugs<-function(sub_df_drugs){
      #genes target by at least one drug
      score_genes_with_drugs<-length(unique(sub_df_drugs[sub_df_drugs[,2]!=0,1]))/length(unique(sub_df_drugs[,1]))
      return(score_genes_with_drugs)
    }
    
    ########################################################################
    
    #1) compute overxpression scores S-ONC
    ge_score_up<-data.frame(t(aggregate(genes_up ~ resKmeans, data = RFKR,s_onc)))[2,]
    colnames(ge_score_up)<-paste("ge_score_up",1:k_user,sep="_")
    
    #2) compute amplification scores  S-ONC
    cnv_score_amp<-data.frame(t(aggregate(cnv ~ resKmeans, data = RFKR, s_onc)))[2,]
    colnames(cnv_score_amp)<-paste("cnv_score_amp",1:k_user,sep="_")
    
    #######################################################################
    
    #3) compute deletion scores S-SUP 
    cnv_score_del<-data.frame(t(aggregate(cnv ~ resKmeans, data = RFKR, cnv_del)))[2,]
    colnames(cnv_score_del)<-paste("cnv_score_dep",1:k_user,sep="_")

    #4) compute methylation scores S-SUP
    meth_score<-data.frame(t(aggregate(meth_expressed ~ resKmeans, data = RFKR, s_onc)))[2,]
    colnames(meth_score)<-paste("meth_score_exp",1:k_user,sep="_")
    
    #5) compute mutations scores S-SUP
    mut_score<-data.frame(t(aggregate(mutation ~ resKmeans, data = RFKR, FUN=mut_count)))[2,]
    colnames(mut_score)<-paste("mut_score",1:k_user,sep="_")
    
    #6) compute low expression scores S-SUP
    ge_score_down<-data.frame(t(aggregate(genes_down ~ resKmeans, data = RFKR, s_onc)))[2,]
    colnames(ge_score_down)<-paste("ge_score_down",1:k_user,sep="_")
    
    ############################################################################
    
    
    # mean of gene-expression
    ge_mean<-data.frame(t(aggregate(ge ~ resKmeans, data = RFKR,mean)))[2,]
    colnames(ge_mean)<-paste("ge_mean",1:k_user,sep="_")
    
    # mean of methylation values
    meth_mean<-data.frame(t(aggregate(meth ~ resKmeans, data = RFKR, mean)))[2,]
    colnames(meth_mean)<-paste("meth_mean",1:k_user,sep="_")
    
    #create a data.frame with hte scores computed
    tab_SCORES<-cbind(t(ge_mean),t(ge_score_up),t(ge_score_down),t(cnv_score_amp),t(cnv_score_del),t(meth_mean),t(meth_score),t(mut_score))
    colnames(tab_SCORES)<-c("mean_ge",
                            "genes_up",
                            "genes_down",
                            "cnv_amp",
                            "cnv_del",
                            "meth_mean",
                            "meth_exp",
                            "mut")
    
    tab_SCORES[is.na(tab_SCORES)]<-0
    
    # compute S-ONC
    S_ONC=(100*tab_SCORES[,"cnv_amp"])+(100*tab_SCORES[,"genes_up"])
    # compute S-SUP
    S_SUP=tab_SCORES[,"mut"]+(100*tab_SCORES[,"meth_exp"])+(100*tab_SCORES[,"cnv_del"])+(100*tab_SCORES[,"genes_down"])
    
    
    S_SUP_DF<-t(data.frame(S_ONC))
    colnames(S_SUP_DF)<-paste("s_sup",1:k_user,sep="_")
    
    S_ONC_DF<-t(data.frame(S_SUP))
    colnames(S_ONC_DF)<-paste("s_onc",1:k_user,sep="_")
    
    S = log(S_ONC/S_SUP,2)
    S[grep(S,pattern="Inf")]<-0
    
    # create a data.frame with the S-SCORES
    GLOBAL_GENETIC_SCORE<-data.frame(t(S))
    # remove inf values
    GLOBAL_GENETIC_SCORE[mapply(is.infinite, GLOBAL_GENETIC_SCORE)] <- 0
    
    # save for which modules GMIEC-MD compute correctly the S-SCORE, true, when the numerator or denominator
    # are 0 
    TRUE_MODULE_INDEX=paste(which(GLOBAL_GENETIC_SCORE!=0),collapse=",")
    
    #WARNING!! Here GMIEC-MD save the values of S_ONC or S_UP nevertheless are 0, these are not S-SCORE are S-ONC or S-SUP
    GLOBAL_GENETIC_SCORE[S_ONC==0]<--log(S_SUP[S_ONC==0],2)
    GLOBAL_GENETIC_SCORE[S_SUP==0]<-log(S_ONC[S_SUP==0],2)
    GLOBAL_GENETIC_SCORE[mapply(is.infinite, GLOBAL_GENETIC_SCORE)] <- 0
    
    colnames(GLOBAL_GENETIC_SCORE)<-paste("s_score",1:k_user,sep="_")
    
    #create the ouput
    number_genes=t(data.frame(apply(GGFK,2,FUN=function(x){unlist(lapply(strsplit(x,split=","),length))})))
    number_drugs=t(data.frame(apply(DGFK,2,FUN=function(x){unlist(lapply(strsplit(x,split="#"),length))})))
    
    #count number of genes with drugs
    SCORES_DRUGS<-NULL
    
    for(i in 1:k_user){
      
      res_drug<-countDrugs(RFKR_drugs[RFKR_drugs$resKmeans==i,c("genes","drugs")])
      
      SCORES_DRUGS<-c(SCORES_DRUGS,res_drug)
      
    }
    
    
    SCORES_DRUGS_DF<-t(data.frame(SCORES_DRUGS))
    colnames(SCORES_DRUGS_DF)<-paste("score_drugs_rdg",1:k_user,sep="_")
    #################
    
    
    #################
    # SAD SCORES
    COMBINED_SCORES=data.frame(GLOBAL_GENETIC_SCORE*SCORES_DRUGS_DF)
    
    colnames(COMBINED_SCORES)<-paste("score_sad",1:k_user,sep="_")
    #################
    TEMP_DF<-data.frame(sample_id=se_patient_selection,
                        GGFK,
                        number_genes=number_genes,
                        DGFK,
                        number_drugs=number_drugs,
                        ge_mean,
                        ge_score_up,
                        cnv_score_amp,
                        cnv_score_del,
                        ge_score_down,
                        meth_mean,
                        meth_score,
                        mut_score,
                        S_SUP_DF,
                        S_ONC_DF,
                        GLOBAL_GENETIC_SCORE,
                        SCORES_DRUGS_DF,
                        COMBINED_SCORES,
                        TRUE_MODULE_INDEX
    )
    
    FINAL_DF<-rbind(FINAL_DF,TEMP_DF)
    
  }
  colnames(input_clinical)[1]<-"sample_id"
  
FINAL_DF_WITH_CLINIC<-merge(FINAL_DF,input_clinical,by="sample_id",all.x=T)
 
return(FINAL_DF_WITH_CLINIC) 
}
