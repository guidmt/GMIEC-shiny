
engine_all_dataset<-function(input_for_apriori2,dfPatientForAnalysis_GAC){

    
    print("Step 6: Run apriori")
    #a confidence of 0.05 is a good solution
    # support: a numeric value for the minimal support of an item set (default: 0.1)
    rules = apriori(input_for_apriori2, parameter=list(support=0.0001, confidence=0.90,minlen=3,maxlen=ncol(input_for_apriori2),
                                                       target = "closed frequent itemset"))
    
    # support: a numeric value for the minimal support of an item set (default: 0.1)
    # rules = apriori(trans , parameter=list(support=0.80, confidence=0.80,
    #                                        target = "closed frequent itemset"))
    
    ruledf = DATAFRAME(rules)
    subrules.sort<-ruledf[order(ruledf$support,decreasing=T),]
   
    #print(dim(subrules.sort))
    #print(head(subrules.sort))
    
    subrules.sort2<-subrules.sort[grep(x=subrules.sort$items,pattern="genesID"),]
    if(nrow(subrules.sort2)>20000){
      subrules.sort2<-head(subrules.sort2,20000)
    }
    ###
    ### Filtering of the rules
    ###
    print("Step 7: Run filtering")
    
    list_rules<-unique(unlist(lapply(X=strsplit(as.character(subrules.sort2[,1]),split=","),FUN=function(X){paste(X[2:length(X)],collapse=",")})))
    
    genes_in_rules<-sapply(strsplit(sapply(strsplit(as.character(subrules.sort2[,1]),split=","),"[[",1),split="="),"[[",2)#
    
    if(length(genes_in_rules)==nrow(input_for_apriori2)){ #if the number of genes detected in rules is the same of the initial rules do the follow processing
      
      
      # Create a list in which the number of element of the list correspond with the number of rule, the content of each elements are the genes with a given rule
      
      association_rules_genes<-sapply(list_rules,FUN=function(x){
        
        gsub(sapply(strsplit(unlist(strsplit(grep(pattern=x,x=subrules.sort2[,1],value=T),split="genesID=")),split=","),"[[",1),pattern="\\{",replacement="")
        
      })
      
      print("Do find association_rules_genes")
      
      #parse the previously results in a data.frame
      ARG<-list(1:length(association_rules_genes))
      
      for(l in 1:length(association_rules_genes)){
        
        current_rule<-names(association_rules_genes[l])
        number_of_genes<-length(association_rules_genes[[l]])
        data_arg_module<-data.frame(rule=rep(current_rule,number_of_genes),genes=association_rules_genes[[l]])
        ARG[[l]]<-data_arg_module
      }
      
      ARG2<-do.call(rbind,ARG)
      ARG2clean<-ARG2[ARG2[,2]!="",]
      ARG2clean[,1]<-as.character(ARG2clean[,1])
      
      #assign for each rule a number (module) that identify same rules with one id
      association_rule_module<-data.frame(rule=names(association_rules_genes),Groups_Apriori=1:length(association_rules_genes))
      association_rule_module[,1]<-as.character(association_rule_module[,1])
      
      ARM<-merge(x=association_rule_module,y=ARG2clean,by="rule")
      
    } else { #statement if, code to decide how to treat the genes with multiple rules: slow to run
      
      ugir<-unique(genes_in_rules)
      
      rules_identified_from_multiple_genes<-sapply(X=ugir,FUN=function(X){ #start function sapply
        
        idx<-grep(subrules.sort2[,1],pattern=X)
        
        subrules.sort.gene<-subrules.sort2[idx,]
        
        if(nrow(subrules.sort.gene)==1){ # if for the current gene there is one rule
          
          rules<-as.character(subrules.sort.gene[,1]) # gmid find the rule and mantain it
          
        } else { # if for one genes are available more rules
          
          maxsupport<-max(subrules.sort.gene[,2]) #check the maximum values of support
          minsupport<-max(subrules.sort.gene[,2]) #check the minimum values of support
          
          if(minsupport!=maxsupport){
            
            rules<-as.character(subrules.sort.gene[which.max(subrules.sort.gene[,2]),1]) #get the rules with higher support value
            
          } else {
            
          }
          nostrings_rules<-unlist(lapply(X=strsplit(as.character(subrules.sort.gene[,1]),split=","),FUN=function(X){ #get the rules with the maximum sum of properties. 
            nostrings_rules<-gsub(X,pattern="[^0-9]",replacement="")                                                   #this is useful because GMID prefer to consider properties with 1 than 0
            sum(as.numeric(nostrings_rules[!nostrings_rules==""]))                                                     #at biological level i prefer features that explain a phenoma.
          }))
          
          idxrule<-which.max(nostrings_rules)
          rules<-as.character(subrules.sort.gene[idxrule,1])
        } #close internal if statement for genes with many rules
        
      } #close function sapply
      
      ) #close sapply
      
      #repeat the same processing of before
      rules_identified_from_multiple_genes_df<- data.frame(genes=names(rules_identified_from_multiple_genes),rules=rules_identified_from_multiple_genes)
      
      list_rules<-unique(unlist(lapply(X=strsplit(as.character(rules_identified_from_multiple_genes),split=","),FUN=function(X){paste(X[2:length(X)],collapse=",")})))
      
      association_rules_genes<-sapply(list_rules,FUN=function(x){
        
        gsub(sapply(strsplit(unlist(strsplit(grep(pattern=x,x=as.character(rules_identified_from_multiple_genes_df[,2]),value=T),split="genesID=")),split=","),"[[",1),pattern="\\{",replacement="")
        
      })
      
      ARG<-list(1:length(association_rules_genes))
      
      for(l in 1:length(association_rules_genes)){
        
        current_rule<-names(association_rules_genes[l])
        number_of_genes<-length(association_rules_genes[[l]])
        data_arg_module<-data.frame(rule=rep(current_rule,number_of_genes),genes=association_rules_genes[[l]])
        ARG[[l]]<-data_arg_module
      }
      
      ARG2<-do.call(rbind,ARG)
      ARG2clean<-ARG2[ARG2[,2]!="",]
      ARG2clean[,1]<-as.character(ARG2clean[,1])
      
      #assign for each rule a number (module) that identify same rules with one id
      association_rule_module<-data.frame(rule=names(association_rules_genes),Groups_Apriori=1:length(association_rules_genes))
      association_rule_module[,1]<-as.character(association_rule_module[,1])
      
      ARM<-merge(x=association_rule_module,y=ARG2clean,by="rule")
      
    } #close statemenent if for one genes are avialable multiple rules
    
    mergeGAC_COM<-merge(dfPatientForAnalysis_GAC,ARM,by.x="genesID",by.y="genes")
    mergeGAC_COM_sort<-mergeGAC_COM[order(mergeGAC_COM$Groups_Apriori,decreasing=F),]
    
    
    ### 
    ### Analysis 2: Use k-means for clustering of gene-expression, copy-number variation, methylation data.
    ###
    
    #mergeGAC_COM_KGE<-kmeans(x=mergeGAC_COM_sort[,"GE_current_patient"],4) deprecated
    #mergeGAC_COM_KCNV<-kmeans(x=mergeGAC_COM_sort[,"CNV_current_patient"],4)
    #mergeGAC_COM_KMETH<-kmeans(x=mergeGAC_COM_sort[,"METH_current_patient"],4)
    
    #mergeGAC_COM_res_K<-cbind(mergeGAC_COM_sort,clusterGE=mergeGAC_COM_KGE$cluster,clusterCNV=mergeGAC_COM_KCNV$cluster,clusterMETH=mergeGAC_COM_KMETH$cluster)
    mergeGAC_COM_res_K<-cbind(mergeGAC_COM_sort)
    
    #obtain the list of clusters
    clusters_list_for_2<-unique(mergeGAC_COM_res_K$Groups_Apriori)
    
    mergeGAC_COM_res_K_2<-mergeGAC_COM_res_K

    return(mergeGAC_COM_res_K_2)
}