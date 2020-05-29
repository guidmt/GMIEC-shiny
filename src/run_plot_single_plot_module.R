plot_heatmap_module<-function(res_gmiec,input_GE,input_CNV,input_METH,input_MUT,input_DRUGS,subject){
  
  subject2<-as.character(subject)
  
  res_gmiec_select<-res_gmiec[res_gmiec[,1]%in%subject2,]
  print(dim(res_gmiec_select))
  genes<-unlist(strsplit(as.character(res_gmiec_select[,'genes_in_module']),split=","))
  input_GE[,-1]<-apply(X=input_GE[,-1],1,FUN=function(X){(X-mean(X))/sd(X)})
  
  input_GE_sel<-data.frame(input_GE[input_GE[,1]%in%genes,subject2])
  rownames(input_GE_sel)<-input_GE[input_GE[,1]%in%genes,1]
  colnames(input_GE_sel)<-"GE"
  
  input_CNV_sel<-input_CNV[input_CNV[,1]%in%genes,subject2]
  input_METH_sel<-input_METH[input_METH[,1]%in%genes,subject2]
  input_MUT_sel<-input_MUT[input_MUT[,2]==subject2,]
  input_MUT_sel2<-input_MUT_sel[input_MUT_sel[,1]%in%genes,]
  
  input_MUT_sel3<-rep(0,length(genes))
  
  if(dim(input_MUT_sel2)[1]==0){
    
    input_MUT_sel3<-rep(0,length(genes))
    
  } else{
      
    
    input_MUT_sel3[input_MUT_sel2[,1]%in%genes]<-1
    
  }
  
  ###count the number of drugs for genes
  VECTOR_DRUGS<-NULL
  
  for(g in genes){
    print(g)
    countDrugs<-length(input_DRUGS[input_DRUGS[,1]==g,2])
    VECTOR_DRUGS<-c(VECTOR_DRUGS,countDrugs)
  }
  print("heatmap ge")
  ht1=Heatmap(input_GE_sel,row_labels=rownames(input_GE_sel),row_names_side="left",cluster_rows = FALSE,cluster_columns = FALSE,name="GE",rect_gp = gpar(col = 'black'),rev(brewer.pal(6,"RdYlBu")),width=0.3)
  
  print("heatmap cnv")
  
  if(unique(input_CNV_sel)!=0){
  ht2=Heatmap(input_CNV_sel,cluster_rows = FALSE,cluster_columns = FALSE,name="CNV",rect_gp = gpar(col = 'black'),rev(brewer.pal(6,"RdBu")),width=0.3)
  } else{
  ht2=Heatmap(input_CNV_sel,cluster_rows = FALSE,cluster_columns = FALSE,name="CNV",width=0.3,col="gainsboro")
  }
  
  print("heatmap meth")
  
  if(unique(input_METH_sel)!=0){
  ht3=Heatmap(input_METH_sel,cluster_rows = FALSE,cluster_columns = FALSE,name="METH",rect_gp = gpar(col = 'black'),rev(brewer.pal(6,"PuOr")),width=0.3)
  }else{
  ht3=Heatmap(input_METH_sel,cluster_rows = FALSE,cluster_columns = FALSE,name="METH",rect_gp = gpar(col = 'black'),width=0.3,col="gainsboro")
  }
  
  print("heatmap mut")
  
  ht4=Heatmap(input_MUT_sel3,cluster_rows = FALSE,cluster_columns = FALSE,name="MUT",rect_gp = gpar(col = 'black'),width=0.3,col="gainsboro")

  print("heatmap drugs")
  
  if(length(unique(VECTOR_DRUGS))!=0){
  ht5=Heatmap(VECTOR_DRUGS,cluster_rows = FALSE,cluster_columns = FALSE,name="DRUGS",rect_gp = gpar(col = 'black'),rev(brewer.pal(6,"YlOrRd")),width=0.3)
  }else{
  ht5=Heatmap(VECTOR_DRUGS,cluster_rows = FALSE,cluster_columns = FALSE,name="DRUGS",rect_gp = gpar(col = 'black'),width=0.3,col="gainsboro")
  }
  
  ht1+ht2+ht3+ht4+ht5
  
  }