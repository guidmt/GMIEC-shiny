####
#### This script contain the R-code and the functions to explore the data from GMIEC analysis.
#### version 11/10/2017




#find what is the output dir to save the results from the config file.
#change with your path in which reside the configuration file.
setwd("C:/Users/guida/Desktop/stratification_tcga") 

conf_file<-read.delim(file="configuration_file2.txt",stringsAsFactors=F,header=T)
output_dir<-as.character(conf_file[conf_file[,"PARAMETERS"]=="output_dir",2])

setwd(output_dir)

#change with the file that you want analyze.
#The current version of this file allowed to manage Analysis_GMIEC_simplified_results_maxSAD.txt and Analysis_GMIEC_simplified_results_minSAD.txt
tab_input<-read.delim(file="Analysis_GMIEC_simplified_resultsmaxSAD.txt")
tab_input<-tab_input[-which(tab_input$PAM50.Subtype=="Normal-like"),] #exclude normali samples

tab_input_min<-read.delim(file="Analysis_GMIEC_simplified_resultsminSAD.txt")
tab_input_min2<-tab_input_min[-which(tab_input_min$PAM50.Subtype=="Normal-like"),] #exclude normali samples

tab_input_main<-t(read.delim(file="Analysis_GMIEC_main_results.txt"))
colnames(tab_input_main)<-as.character(tab_input_main[1,])
tab_input_main<-tab_input_main[,-1]

tab_input_main2<-tab_input_main[-which(tab_input_main[,"PAM50.Subtype"]=="Normal-like"),] #exclude normali samples



###
### Define the function number 1: Build a matrix in which the rows correspond to the patients and the columns the number of drugs
###
### tab_input: the data.frame with the data uploaded in the previously step
### clinical_data: define what is the column with the clinical data (e.g. PAM50.subtype), use the same column name of the matrix. NA values not allowed.
###                NA values will be removed.
### output_file: define the string to save the data
###
### Note: The function iterate for each clinical variable. For example if you considered PAM50 subtype, 
###       a number of output file equal to the number of cancer subtypes will be created.
###       The file will be saved with the suffix of each clinical variable. In the case of PAM50 the file created will be have      
###       the suffix ".LumA", "LumB","Basal" etc.
###
### Output Results:
### - For each clinical variable a binary heatmap with the patients and drugs
### - For each clinical variable a matrix with the results, rows and columns sorted according the clustering results 
### - For each clinical variable a two column matrix, the first column contains the list of drugs,
###   the second column the frequency of presence of a drug along the patients. The rows are sorted according with the clustering results
###
###

createPatientDrugMatrix<-function(tab_input,clinical_data,output_file){
    
      uc<-as.character(unique(tab_input[,clinical_data]))
     
      unique_clinical<-  uc[!is.na(uc)]
         
  for(un in unique_clinical){
      
    
      ##extract all compounds
      print(un)
        
      subtab_input<-tab_input[which(tab_input[,clinical_data]==un),]
      
      strdrug1<-unlist(strsplit(as.character(subtab_input[,3]),split="#"))
      strdrug2<-unique(unlist(strsplit(as.character(strdrug1),split="@")))
      
      
      #build an empty matrix with nrow=samples, ncol=drugs
      matrix_patients_drugs<-data.frame(matrix(0,nrow=nrow(subtab_input),ncol=length(strdrug2)+1))
      matrix_patients_drugs[,1]<-as.character(subtab_input$patientID)
      colnames(matrix_patients_drugs)[2:ncol(matrix_patients_drugs)]<-as.character(strdrug2)
      
      ##extract all genes
      subtab_inputgenes<-tab_input[which(tab_input[,clinical_data]==un),]
      
      strgenes1<-unique(unlist(strsplit(as.character(subtab_inputgenes[,2]),split=",")))

      
      #build an empty matrix with nrow=samples, ncol=drugs
      matrix_patients_genes<-data.frame(matrix(0,nrow=nrow(subtab_inputgenes),ncol=length(strgenes1)+1))
      matrix_patients_genes[,1]<-as.character(subtab_inputgenes$patientID)
      colnames(matrix_patients_genes)[2:ncol(matrix_patients_genes)]<-as.character(strgenes1)
      

  
  ###
  ### slow to run: this is a matrix in which the number of rows correpond with the patients and the number of columns the drugs
  ###
  
  #I'm considering the module with highest SAD along the patients in a particular cancer subtype
  # - inside each cancer subtype there are different in treatment 
  for(i in 1:length(matrix_patients_drugs[,1])){
    
    cpatient<-as.character(matrix_patients_drugs[i,1])
  
    cdrugs<-subtab_input[which(subtab_input$patientID==cpatient),3] #get the data of the drugs (drugs_in_module_With_maxSAD/drugs_in_module_With_minSAD)
    cdrug1<-unlist(strsplit(as.character(cdrugs),split="#"))
    cdrug2<-unique(unlist(strsplit(as.character(cdrug1),split="@")))
  
  resDrugsCurrentP<-sapply(X=colnames(matrix_patients_drugs[,-1]),FUN=function(X,cpat=cpatient){
      
    ok<-ifelse(length(!which(!is.na(match(cdrug2,X))))==1,1,0)
    
    })
    
  matrix_patients_drugs[i,-1]<-resDrugsCurrentP
  
  }
  
  #do clustering before
  order_column=dist(t(matrix_patients_drugs[,-1]), method="binary")
  order_column_hclust<- hclust(order_column, "average")
      
  order_row=dist((matrix_patients_drugs[,-1]), method="binary")
  order_row_hclust<- hclust(order_row, "average")
      
  mat2temp<-data.frame(patientID=matrix_patients_drugs[order_row_hclust$order,1],matrix_patients_drugs[,-1][order_row_hclust$order,order_column_hclust$order],stringsAsFactors=F)
  
  file_output1<-paste(paste(paste(output_file,".matrix.patients_drugs",sep=""),un,sep="."),".txt",sep="")
  write.table(mat2temp,file=file_output1,sep="\t",row.names=FALSE,quote=FALSE)
      
  res_for_annotation<-apply(X=matrix_patients_drugs[,-1],2,FUN=function(X){sum(X)/length(X)})
  
  file_output2<-paste(paste(paste(output_file,".heatmap.patients_drugs",sep=""),un,sep="."),".pdf",sep="")

  column_ha = HeatmapAnnotation(points = anno_points(res_for_annotation, axis = TRUE,size=unit(0.5,"mm")),annotation_height=unit(30,"mm"))
  #see the documenttion of complexHeatmap about the options to sort columns and rows
  h<- Heatmap(matrix_patients_drugs[,-1],col=c("#FFFFFF","#3794bf"),top_annotation = column_ha,show_row_names =FALSE,show_column_names=FALSE,row_dend_reorder = FALSE,column_dend_reorder = FALSE,cluster_rows=order_row_hclust,cluster_columns =order_column_hclust)
  
  order_of_the_drugs<-order_column_hclust$order
  rfad.sort<-res_for_annotation[order_of_the_drugs]
  names_compounds<-as.character(names(rfad.sort))
  length_na<-length(which(is.na(names_compounds)))
  names_compounds[is.na(names_compounds)]<-paste("fixedNA",1:length_na,sep=".")#fixed with NA values if present
    
  names(rfad.sort)<-names_compounds #fixed with NA values if present
  
  file_output3<-paste(paste(paste(output_file,".twocm.clustering.patients_drugs",sep=""),un,sep="."),".txt",sep="")
  write.table(rfad.sort,file=file_output3,sep="\t",row.names=TRUE,quote=FALSE)
 

  pdf(file=file_output2,width=20)
  print(h)
  dev.off()
  
  ####
  #### repeat the step before for genes
  ####
  for(i in 1:length(matrix_patients_genes[,1])){
    
    cpatient<-as.character(matrix_patients_genes[i,1])
    
    cgenes<-subtab_input[which(subtab_input$patientID==cpatient),2] #get the data of the genes (genes_in_module_With_maxSAD/genes_in_module_With_minSAD)
    cgene1<-unlist(strsplit(as.character(cgenes),split=","))

    resgenesCurrentP<-sapply(X=colnames(matrix_patients_genes[,-1]),FUN=function(X,cpat=cpatient){
      
      ok<-ifelse(length(!which(!is.na(match(cgene1,X))))==1,1,0)
      
    })
    
    matrix_patients_genes[i,-1]<-resgenesCurrentP
    
  }
  
  #do clustering before
  order_column=dist(t(matrix_patients_genes[,-1]), method="binary")
  order_column_hclust<- hclust(order_column, "average")
  
  order_row=dist((matrix_patients_genes[,-1]), method="binary")
  order_row_hclust<- hclust(order_row, "average")
  
  mat2temp<-data.frame(patientID=matrix_patients_genes[order_row_hclust$order,1],matrix_patients_genes[,-1][order_row_hclust$order,order_column_hclust$order],stringsAsFactors=F)
  
  file_output1<-paste(paste(paste(output_file,".matrix.patients_genes",sep=""),un,sep="."),".txt",sep="")
  write.table(mat2temp,file=file_output1,sep="\t",row.names=FALSE,quote=FALSE)
  
  res_for_annotation<-apply(X=matrix_patients_genes[,-1],2,FUN=function(X){sum(X)/length(X)})
  
  file_output2<-paste(paste(paste(output_file,".heatmap.patients_genes",sep=""),un,sep="."),".pdf",sep="")
  
  column_ha = HeatmapAnnotation(points = anno_points(res_for_annotation, axis = TRUE,size=unit(0.5,"mm")),annotation_height=unit(30,"mm"))
  #see the documenttion of complexHeatmap about the options to sort columns and rows
  h<- Heatmap(matrix_patients_genes[,-1],col=c("#FFFFFF","#cb003e"),top_annotation = column_ha,show_row_names =FALSE,show_column_names=FALSE,row_dend_reorder = FALSE,column_dend_reorder = FALSE,cluster_rows=order_row_hclust,cluster_columns =order_column_hclust)
  
  order_of_the_genes<-order_column_hclust$order
  rfad.sort<-res_for_annotation[order_of_the_genes]
  names_genes<-as.character(names(rfad.sort))
  length_na<-length(which(is.na(names_genes)))
  names_genes[is.na(names_genes)]<-paste("fixedNA",1:length_na,sep=".")#fixed with NA values if present
  
  names(rfad.sort)<-names_genes #fixed with NA values if present
  
  file_output3<-paste(paste(paste(output_file,".twocm.clustering.patients_genes",sep=""),un,sep="."),".txt",sep="")
  write.table(rfad.sort,file=file_output3,sep="\t",row.names=TRUE,quote=FALSE)
  
  
  pdf(file=file_output2,width=20)
  print(h)
  dev.off()
  
  
  
      
    }#close internal loop for clinical variable
      
}


###
### Detine the third function to explore in how many patients a gene is identify with a SADthr
###

#This function for each gene, patient. Check for a given SAD threshold if a gene is identified or not along or patient.
#The resulting output is a matrix in which the rows represnet the genes, the columns the threholds. The values of the matrix
#are the frequency of identification of a gene along all patients. 
# input:
# - tab_input2: the principal output of GMID Analysis_GMID_main_results.txt 
# - patients: a vector with the code of patients
# - list_genes: a vector with the list og genes
# - vector_sad_thr: a numeric vector with the list of thr SAD
#
#  output: 
# - A matrix in which the rows represnet the genes, the columns the threholds. 
  
  
  FreqGenesInPatientForSadThr<-function(tab_input2,patients,list_genes2,vector_sad_thr){
    
    subtab<-tab_input2[rownames(tab_input2)%in%patients,]
    subtab2<-cbind(patientIDforanalysis=rownames(subtab),subtab)
    
    results_genes_thr<-vector("list",length(vector_sad_thr)) 
    
    for(th in 1:length(vector_sad_thr)){
      
      thr<-vector_sad_thr[th]
      
      print(thr)
      genes_patient_df_current_thr<-data.frame()
      
      for(pt in samples2){
        
        idxsad<-grep(colnames(subtab2),pattern="sad",value=T)
        idxgenes<-grep(colnames(subtab2),pattern="genes_in_module",value=T)
        
        resGenesThrPatient<- sapply(X=list_genes2,function(X){
          
          #in this step i parse the list of genes for each columns using strsplit, this allowed to correctly select the genes with few characters es. C6 
          idxgenesgrep<-unlist(lapply(sapply(subtab2[subtab2[,1]==pt,idxgenes],strsplit,split=","),grep,pattern=paste(paste("^",X,sep=""),"$",sep="")))
          idxgenesgrep<-as.numeric(gsub(names(idxgenesgrep),pattern="genes_in_module",replacement=""))
          
          resGeneThr<-ifelse(subtab2[subtab2[,1]==pt,idxsad][idxgenesgrep] >=thr,1,0)
          return(resGeneThr)
        }
        )#close sapply #for each patient and all genes i obtain a binary vector 1 the genes is present at a given threshold, otherwies not.
        
        genes_patient_df_current_thr<-rbind(genes_patient_df_current_thr,resGenesThrPatient)
        
      }
      
      genes_patient_df_current_thr_t<-t(genes_patient_df_current_thr)
      
      genes_patient_df_current_thr_t_freq <-apply(X=genes_patient_df_current_thr_t,1,FUN=function(X){sum(X)/length(samples2)}) #estimate the frequency of occuring of a genes with a given thr,.
      names(genes_patient_df_current_thr_t_freq)<-list_genes2
      
      results_genes_thr[[th]]<-genes_patient_df_current_thr_t_freq
      
    }
    resGenesFrequencyforThr<-do.call(rbind,results_genes_thr)
    resGenesFrequencyforThr_t<-t(resGenesFrequencyforThr)
    colnames(resGenesFrequencyforThr_t)<-vector_sad_thr
    
    return(resGenesFrequencyforThr_t)
  }
  
#First: Extract the results from patients-drugs analysis
tab_input<-tab_input
clinical_data<-"PAM50.Subtype" #use the same string as in colnames

createPatientDrugMatrix(tab_input=tab_input,clinical_data=clinical_data,output_file="breast")
createPatientDrugMatrix(tab_input=tab_input_min2,clinical_data=clinical_data,output_file="breast.min") #min sad


#
# Control for a given set of genes and patients the frequency of identification of a gene considering different sad thr.
#


list_genes<-c("ATP8A2,PSCA,GATA2,NRTN,CXADR,CBFA2T3,RBM24,CD96,IL20RA,MOGAT2,TNFRSF17,CHI3L2,
              LAMB3,PDZK1,PNLIPRP2,CD2,COL1A2,BMP4,GPC3,KLHL13,MID1,WISP1,ATP6V0A4,SCML1,RGS1,
              MYBL2,OCA2,PNMT,LDHC,CCL8,VTCN1,STMN2,CCL28,CLCA2,ADRA2C,MMP13,TRIM2,KLRG2,MMP12,
              BFSP2,CALML5,CYP27B1,GREB1,DCN,RUNX1T1,OLFM4,DPT,SELP,RGS18,ADH1A,FMO1,THRSP,IGF1,
              CD38,GAS2,SLC13A3,PGBD5,AGT,COL4A6,COLEC11,PNOC,ABHD12B,BNC1,SLC24A5,WNT5B,PAQR5,SPINK1,
              CEACAM5,S100P,SEMA3B,TSPAN8,PKIA,C1S,TMC5,AMDHD1,CFTR,NEFH,GPR34,CD69,NPTX2,MSMB,APOD,C2orf40,
              KRT6B,CCL19,MMP7,MOXD1,ZNF683,KLRB1,ANKRD30A,CLEC4D,NGEF,CPXM2,C1GALT1,NLRP7,POU2AF1,
              GFI1,GPR65,LCK,SLAIN1,CHRNA5,PLAC8,CPA3,ALPK2,CP,QPRT,POSTN,FGFBP1,ACTL8,TRIM63,XK,
              IGSF10,CD3D,LRP8,FANCA,FABP7,KLRC2,FABP6,DDX43,HSD17B2,CD28,FABP4,CHST6,ASRGL1,
              HIST1H1B,HIST1H4L,SYT1,CX3CR1,ABCC13,IL20,KCNJ11,XDH,CCND1,FHL2,C16orf45,ATP6V1B1,
              C8orf4,CLIC6,PPEF1,KRT7,C4orf19,CD53,MMRN1,CHRDL2,MAPK8IP2,REPS2,SFTPD,C6,IYD,
              ELF5,IFNG,SDCBP2,GJB5,SPESP1,ANPEP,MPPED2,DLGAP1,EPYC,FBXO2,CRYM,IL13RA2,
              ACMSD,BCL2A1,HIST1H1A,SMARCA1,HAPLN1,TUBAL3,IL21R,SPOCD1,DPP4,CLDN14,ITK,
              CD19,LYZ,MS4A1,LGALS2,SIT1,UBASH3A,TBC1D10C,CD5,KLHL6,ACSL5,RUNX3,SLA2,SPOCK2")

list_genes2<-gsub(unlist(strsplit(unlist(strsplit(list_genes,split=",")),"\n ")),pattern=" ",replacement="") #i used this approaches for fast, the user can upload an external file
list_genes2<-list_genes2[!list_genes2==""]

#sample luminal A
samples<-c("TCGA.A2.A04N.01,TCGA.A2.A04V.01,TCGA.A2.A04Y.01,TCGA.A2.A0CP.01,TCGA.A2.A0CQ.01,TCGA.A2.A0CS.01,TCGA.A2.A0CU.01,TCGA.A2.A0CV.01,TCGA.A2.A0D3.01,TCGA.A2.A0EM.01,TCGA.A2.A0EO.01,TCGA.A2.A0ES.01,TCGA.A2.A0ET.01,TCGA.A2.A0EV.01,TCGA.A2.A0EW.01,TCGA.A2.A0EX.01,TCGA.A7.A0CD.01,TCGA.A7.A0CG.01,TCGA.A7.A0CH.01,TCGA.A7.A0DB.01,TCGA.A8.A06P.01,TCGA.A8.A06T.01,TCGA.A8.A06U.01,TCGA.A8.A06Y.01,TCGA.A8.A07E.01,TCGA.A8.A07F.01,TCGA.A8.A07G.01,TCGA.A8.A07J.01,TCGA.A8.A07P.01,TCGA.A8.A083.01,TCGA.A8.A086.01,TCGA.A8.A08A.01,TCGA.A8.A08C.01,TCGA.A8.A08T.01,TCGA.A8.A08Z.01,TCGA.A8.A090.01,TCGA.A8.A091.01,TCGA.A8.A093.01,TCGA.A8.A099.01,TCGA.A8.A09A.01,TCGA.A8.A09B.01,TCGA.A8.A09T.01,TCGA.A8.A09V.01,TCGA.A8.A0A1.01,TCGA.A8.A0A2.01,TCGA.A8.A0A4.01,TCGA.AN.A03X.01,TCGA.AN.A046.01,TCGA.AN.A04A.01,TCGA.AN.A0FD.01,TCGA.AN.A0FN.01,TCGA.AN.A0FS.01,TCGA.AN.A0FT.01,TCGA.AN.A0FW.01,TCGA.AN.A0FZ.01,TCGA.AO.A03V.01,TCGA.AO.A0J8.01,TCGA.AO.A0J9.01,TCGA.AO.A12A.01,TCGA.AO.A12H.01,TCGA.B6.A0I5.01,TCGA.B6.A0IA.01,TCGA.B6.A0IH.01,TCGA.B6.A0IN.01,TCGA.B6.A0IO.01,TCGA.B6.A0IP.01,TCGA.B6.A0RQ.01,TCGA.B6.A0WS.01,TCGA.B6.A0X0.01,TCGA.BH.A0B0.01,TCGA.BH.A0B4.01,TCGA.BH.A0BO.01,TCGA.BH.A0BP.01,TCGA.BH.A0BQ.01,TCGA.BH.A0BR.01,TCGA.BH.A0BV.01,TCGA.BH.A0C1.01,TCGA.BH.A0DE.01,TCGA.BH.A0DO.01,TCGA.BH.A0DT.01,TCGA.BH.A0DX.01,TCGA.BH.A0E7.01,TCGA.BH.A0E9.01,TCGA.BH.A0EA.01,TCGA.BH.A0EB.01,TCGA.BH.A0EI.01,TCGA.BH.A0H5.01,TCGA.BH.A0HO.01,TCGA.BH.A0W7.01,TCGA.BH.A18H.01,TCGA.BH.A18I.01,TCGA.BH.A18M.01,TCGA.BH.A18N.01,TCGA.BH.A18S.01,TCGA.BH.A18T.01,TCGA.C8.A12N.01,TCGA.C8.A12Y.01,TCGA.C8.A132.01,TCGA.C8.A133.01,TCGA.D8.A141.01,TCGA.D8.A145.01,TCGA.E2.A14Q.01,TCGA.E2.A14Z.01,TCGA.E2.A153.01,TCGA.E2.A154.01,TCGA.E2.A156.01,TCGA.E2.A15C.01,TCGA.E2.A15G.01,TCGA.E2.A15H.01,TCGA.E2.A15O.01,TCGA.E2.A15P.01,TCGA.E2.A15R.01")
samples2<-unlist(strsplit(samples,split=","))

#define the SAD thr
vector_sad_thr<-seq(0,2,by=0.11)

tab_input2=tab_input_main2
list_genes2=list_genes2
vector_sad_thr=vector_sad_thr
patients=samples2

#Run the function: slow
resFGPSthr<-FreqGenesInPatientForSadThr(tab_input2=tab_input_main2,list_genes2=list_genes2,patients=samples2,vector_sad_thr=vector_sad_thr)

#draw the heatmap
pdf("res_freqgenes_sad.pdf",width=16,height=6)
Heatmap(resFGPSthr,cluster_rows=F,cluster_columns=F,row_names_gp = gpar(fontsize = 5))
dev.off()
