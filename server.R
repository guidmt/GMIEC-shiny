options(shiny.maxRequestSize=100*1024^2)

source("./src/create_output.R", local = TRUE)
source("./src/engine_all_dataset.R", local = TRUE)
source("./src/filter_ge.R", local = TRUE)
source("./src/internal_annotation.R", local = TRUE)
source("./src/rules_for_tf.R", local = TRUE)
source("./src/rules_fornot_tf.R", local = TRUE)
source("./src/GMIEC.R", local = TRUE)
source("./src/create_report.R", local = TRUE)


function(input,output,session) { 
  
  gmiec_results<-observeEvent(input$run_gmiec,{
  
  ###
  ### Upload experiment datas
  ###
    
  ##2
  input_GE<-reactive({
  infile2<- input$ge_dataset
  read.table(file=infile2$datapath,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=F) #read empty values with 0 
  })
  
  print(paste("Load gene-expression data",dim(input_GE())))
  
  ##3
  input_CNV<-reactive({
  infile3<-input$cnv_dataset
  read.table(file=infile3$datapath,sep="\t",stringsAsFactors=F,header=T,quote=NULL,fill=T)
  })
  
  print(paste("Load cnv data",dim(input_CNV())))
  
  ##4
  input_MUT<-reactive({
    infile4<-input$mut_dataset
    read.table(file=infile4$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T)
  })
  
  print(paste("Load mut data",dim(input_MUT())))
  
  ##5
  input_METH<-reactive({
    infile5<-input$meth_dataset
    read.table(file=infile5$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T)
  })
  
  print(paste("Load meth data",dim(input_METH())))
  
  ##6
  input_clinical<-reactive({
    infile6<-input$clinical_dataset
    read.table(file=infile6$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T)
  })
  
    print(paste("Load clinical data",dim(input_clinical())))

  ##7
  input_clusters<-reactive({
    infile7<-input$clusters
  })
  
  print(paste("Number Clusters",input_clusters()))
  
  
  ###
  ### Input annotation
  ###
  
  genes_annotated_fv<-renderText({ input$genes_annotated})
  genes_annotated_TF_fv<-renderText({ input$genes_annotated_TF})
  
  #if you want use a bed file then are required an annotation file,max-gap, bed file
  if(genes_annotated_fv()==TRUE | genes_annotated_TF_fv()==TRUE){
  annotation_dataset<-reactive({
    infile9<-input$annotation_dataset
    read.table(file=infile9$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T)
  })
  
  print(paste("Load annotation data",dim(annotation_dataset())))
  
  }
  
  if(genes_annotated_fv()==TRUE | genes_annotated_TF_fv()==TRUE){
    
  distance_max<-input$distance
  print(paste("This is the distance",distance_max))
  
  }
  
  if(genes_annotated_fv()==TRUE | genes_annotated_TF_fv()==TRUE){
    
    bed_dataset<-reactive({
      infile_bed<-input$bed_file
      read.table(file=infile_bed$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T)
    })   
    print(paste("Upload bed file",dim(bed_dataset())))
  }  
  
  ###
  ### Input Drugs
  ###
  
  drugs_for_analysis<-reactive({
    infile1<-input$drugs_dataset
    read.table(file=infile1$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T)
  })
  
  print(paste("Upload drugs file",dim(drugs_for_analysis())))
  
  ###
  ### Input Type of analysis
  ###
  input_list_of_genes_test<-renderText({ input$list_of_genes})
  
  #if the user want use a custom list, upload it!  
  if(input_list_of_genes_test()==TRUE){
    
   list_genes_for_analysis<-reactive({
     infile_LG<-input$list_of_genes2
     read.table(file=infile_LG$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T)
     
   })
   
   print(paste("Upload list genes file",dim(list_genes_for_analysis())))
   
  }
   

 genes_annotated_TF2_test<-renderText({ input$genes_annotated_TF})
  
  #if the use want evalute the impact of TF
   if(genes_annotated_TF2_test()==TRUE) {
      
      # extract the gene expression data for the current TF
      input_GE_tf<-input_GE()[input_GE()[,1]== input$name_tf,]
      input_CNV_tf<-input_CNV()[input_CNV()[,1]==input$name_tf,]
      input_MUTATION_tf<-input_MUT()[input_MUT()[,1]==input$name_tf,]
      input_METH_tf<-input_METH()[input_METH()[,1]==input$name_tf,]
       }
  
  all_genes_test<-renderText({ input$all_genes})
  
  #if the use want use the entire list of genes
  if(all_genes_test()==TRUE) {
    
    list_genes_dataset<-unique(c(input_GE()[,1],input_CNV()[,1],input_MUT()[,1],input_METH()[,1]))
    all_genes_for_analysis<-list_genes_dataset
    print(length(all_genes_for_analysis))
    
  } 
  
  output_file<-input$output_file
  #
  # Genes annotated and not TF
  #  
  if(genes_annotated_fv()==TRUE){
    print("run annotation")
    
    print(dim(bed_dataset()))
    print(dim(annotation_dataset()))
    
    showNotification("Start annotation!",type="message")
    
    list_genes_for_analysis<-internal_annotation(bed_dataset(),annotation_dataset(),distance_max)
    
    print(length(list_genes_for_analysis))  
    
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()
    
    showNotification("Start annotation! wait...",type="message")
    
    list_genes_for_analysis<-internal_annotation(bed_dataset(),annotation_dataset(),distance_max)
    
    input_GE_selected<-input_GE2[input_GE2[,1]%in%list_genes_for_analysis,]
    input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%list_genes_for_analysis,]
    input_METH_selected<-input_METH2[input_METH2[,1]%in%list_genes_for_analysis,]
    input_MUTATION_selected<-input_MUTATION2[input_MUTATION2[,1]%in%list_genes_for_analysis,]
    print(dim(input_MUTATION_selected))
    
    print("run gmiec!")
    
    output_gmiec<-reactive({run_GMIEC(
      
      check_ge_for_patients=input_GE_selected,
      input_CNV_selected=input_CNV_selected,
      input_METH_selected=input_METH_selected,
      input_MUTATION_selected=input_MUTATION_selected,
      tabDrugs=drugs_for_analysis2,
      input_clinical=input_clinical2,
      parameter_discr=c("1.5;1;0.5"),
        clusters=input_clusters(),
      genes_annotated_TF_fv=FALSE
      
    )})
      
      if(!is.null(output_gmiec())){
        
        showNotification("You can download the results of analysis!",type="message")
        print("here your output")
        print(dim(output_gmiec()))  
        print(output_file)
        
        output_gmiec2<-output_gmiec()
        print(dim(output_gmiec2))
        output$downloadData <- downloadHandler(
          
          filename = paste(output_file,".txt",sep="")
          ,
          content = function(file) {
            write.table(t(output_gmiec2[-1,]),file,sep="\t",row.names=T,col.names=T,quote=F)
          }
        )
        
      }
    
  } 
  
  
    #
    # Genes annotated and  TF
    #  

    if(genes_annotated_TF_fv()==TRUE){
    print("run annotation")
    print(output_file)
    
    showNotification("Start annotation!",type="message")
      
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()

    showNotification("Start annotation! wait...",type="message")
    
    list_genes_for_analysis<-internal_annotation(bed_dataset(),annotation_dataset(),distance_max)
    
    print(length(list_genes_for_analysis))  
    
    input_GE_selected<-input_GE2[input_GE2[,1]%in%list_genes_for_analysis,]
    input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%list_genes_for_analysis,]
    input_METH_selected<-input_METH2[input_METH2[,1]%in%list_genes_for_analysis,]
    input_MUTATION_selected<-input_MUTATION2[input_MUTATION2[,1]%in%list_genes_for_analysis,]
    print(dim(input_MUTATION_selected))
    print("run gmiec!")
    
    output_gmiec<-reactive({run_GMIEC(
          
          check_ge_for_patients=input_GE_selected,
          input_CNV_selected=input_CNV_selected,
          input_METH_selected=input_METH_selected,
          input_MUTATION_selected=input_MUTATION_selected,
          tabDrugs=drugs_for_analysis2,
          input_clinical=input_clinical2,
          parameter_discr=c("1.5;1;0.5"),
          input_GE_tf=input_GE_tf,
          input_CNV_tf=input_CNV_tf,
          input_MUTATION_tf=input_MUTATION_tf,
          input_METH_tf=input_METH_tf,
            clusters=input_clusters(),
          genes_annotated_TF_fv=TRUE
      
    )})
    
    
    if(!is.null(output_gmiec())){

    showNotification("You can download the results of analysis!",type="message")
    print("here your output")
    print(dim(output_gmiec()))  
    print(output_file)
    
    output_gmiec2<-output_gmiec()
    print(dim(output_gmiec2))
    output$downloadData <- downloadHandler(
      
      filename = paste(output_file,".txt",sep="")
      ,
      content = function(file) {
        write.table(t(output_gmiec2[-1,]),file,sep="\t",row.names=T,col.names=T,quote=F)
      }
    )
  
    }
  
  }
  
  #
  # Genes list 
  #  
  
  if(input_list_of_genes_test()==TRUE) {
    
    subsetAnnoDF_unique<-list_genes_for_analysis()[,1]
    
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()
    
    input_GE_selected<-input_GE2[input_GE2[,1]%in%subsetAnnoDF_unique,]
    input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%subsetAnnoDF_unique,]
    input_METH_selected<-input_METH2[input_METH2[,1]%in%subsetAnnoDF_unique,]
    input_MUTATION_selected<-input_MUTATION2[input_MUTATION2[,1]%in%subsetAnnoDF_unique,]
    
    print("run gmiec!")
    
    output_gmiec<-reactive({run_GMIEC(
      
      check_ge_for_patients=input_GE_selected,
      input_CNV_selected=input_CNV_selected,
      input_METH_selected=input_METH_selected,
      input_MUTATION_selected=input_MUTATION_selected,
      tabDrugs=drugs_for_analysis2,
      input_clinical=input_clinical2,
      parameter_discr=c("1.5;1;0.5"),
        clusters=input_clusters(),
      genes_annotated_TF_fv=FALSE
      
    )})
      
      if(!is.null(output_gmiec())){
        
        showNotification("You can download the results of analysis!",type="message")
        print("here your output")
        print(dim(output_gmiec()))  
        print(output_file)
        
        output_gmiec2<-output_gmiec()
        print(dim(output_gmiec2))
        output$downloadData <- downloadHandler(
          
          filename = paste(output_file,".txt",sep="")
          ,
          content = function(file) {
            write.table(t(output_gmiec2[-1,]),file,sep="\t",row.names=T,col.names=T,quote=F)
          }
        )
        
      }
    
  }
  
  #
  # Genes list  all
  #  
  
  if(all_genes_test()==TRUE) {
    
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()
    
    input_GE_selected<-input_GE2[input_GE2[,1]%in%all_genes_for_analysis,]
    input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%all_genes_for_analysis,]
    input_METH_selected<-input_METH2[input_METH2[,1]%in%all_genes_for_analysis,]
    input_MUTATION_selected<-input_MUTATION2[input_MUTATION2[,1]%in%all_genes_for_analysis,]
  
    print("run gmiec!")
    
    output_gmiec<-reactive({run_GMIEC(
      
      check_ge_for_patients=input_GE_selected,
      input_CNV_selected=input_CNV_selected,
      input_METH_selected=input_METH_selected,
      input_MUTATION_selected=input_MUTATION_selected,
      tabDrugs=drugs_for_analysis2,
      input_clinical=input_clinical2,
      parameter_discr=c("1.5;1;0.5"),
        clusters=input_clusters(),
      genes_annotated_TF_fv=FALSE
      
    )})
    
    if(!is.null(output_gmiec())){
      
      showNotification("You can download the results of analysis!",type="message")
      print("here your output")
      print(dim(output_gmiec()))  
      print(output_file)
      
      output_gmiec2<-output_gmiec()
      print(dim(output_gmiec2))
      output$downloadData <- downloadHandler(
        
        filename = paste(output_file,".txt",sep="")
        ,
        content = function(file) {
          write.table(t(output_gmiec2[-1,]),file,sep="\t",row.names=T,col.names=T,quote=F)
        }
      )
      
    }
    
    }
  }#end function obeserve event
 ) #end observed event1
  
  #########################################
  # VIS-GMIEC
  #########################################
  
  observeEvent(input$run_vis,{
    
  input_for_report2<-reactive({
    infile_for_report<- input$vis_gmiec2
    read.table(file=infile_for_report$datapath,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=T) #read empty values with 0 
    })
  
  print(dim(input_for_report2()))
  
    print("Run creation report!")
    
   # resReport<-reactive({
      run_create_output(input_for_report=input_for_report2())
    #})
    
  })  
  
}#end function
