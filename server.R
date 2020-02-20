options(shiny.maxRequestSize=100*1024^2)

source("./src/create_output.R", local = FALSE)
source("./src/engine_all_dataset.R", local = FALSE)
source("./src/filter_ge.R", local = FALSE)
source("./src/internal_annotation.R", local = FALSE)
source("./src/rules_fornot_tf.R", local = FALSE)
source("./src/GMIEC.R", local = FALSE)
source("./src/gmiec_ML.R",local=FALSE)
source("./src/plot_shiny_functions.R",local=FALSE)
source("./src/simplify_output.R",local=FALSE)
source("./src/run_plot_single_plot_module.R",local=FALSE)

function(input,output,session) { 
  
gmiec_results<-observeEvent(input$run_gmiec,{
  
    GMR<-renderText({input$GMIEC_RULES})
    RFA<-renderText({input$RF_ANALYSIS})
    TD<-renderText({input$two_datasets})
    
    cb_ge_ok <- renderText({ input$cb_ge })
    cb_cnv_ok <- renderText({ input$cb_cnv })
    cb_meth_ok <- renderText({ input$cb_meth })
    cb_mutation_ok <- renderText({ input$cb_mutation })
    
    print(cb_ge_ok())
    print(cb_cnv_ok())
    print(cb_meth_ok())
    print(cb_mutation_ok())
    
    ###
    ### Upload experiment datas
    ###
    
    ### read gene expression data
      input_GE<-reactive({
      showNotification("Loading gene-expression data",type="message")
        
      if(is.null(input$ge_dataset)) {return(NULL)}
      infile2<- input$ge_dataset
      read.table(file=infile2$datapath,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=F,check.names = F) #read empty values with 0 
      
      })
     
    print(paste("Loading gene-expression data",length(input_GE())))
    

    ### read copy number variation data
      input_CNV<-reactive({
      showNotification("Loading copy number data",type="message")
      
      if(is.null(input$cnv_dataset)) {return(NULL)}
      infile3<-input$cnv_dataset
      read.table(file=infile3$datapath,sep="\t",stringsAsFactors=F,header=T,quote=NULL,fill=T,check.names = F)
      })
    
    print(paste("Loading cnv data",length(input_CNV())))
    
    ### read mutation data
  
      input_MUT<-reactive({
      showNotification("Load variants data",type="message")
        
      if(is.null(input$mut_dataset))  {return(NULL)}
      infile4<-input$mut_dataset
      read.table(file=infile4$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
      })    
    
    print(paste("Loading mut data",length(input_MUT())))
    
    ## read methylation data
    input_METH<-reactive({
      
      showNotification("Load methylation data",type="message")
      
      if(is.null(input$meth_dataset))  {return(NULL)}
      infile5<-input$meth_dataset
      read.table(file=infile5$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
    })

    print(paste("Loading meth data",length(input_METH())))
    
    
    ## read clinical data
    input_clinical<-reactive({
    showNotification("Loading clinical data",type="message")
      
    if(is.null(input$clinical_dataset))  {showNotification("The clinical file it is mandatory.",type="error")}
      infile6<-input$clinical_dataset
      read.table(file=infile6$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
    })
    
    print(paste("Loading clinical data",dim(input_clinical())))
    
    ###
    ### Input Drugs
    ###
    drugs_for_analysis<-reactive({
      infile1<-input$drugs_dataset
      read.table(file=infile1$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
    })
    
    print(paste("Upload drugs file",dim(drugs_for_analysis())))
    
    
  
    ## define number of clusters for k-modes and k-means
    input_clusters<-reactive({
    if(is.null(input$clusters))  {return(NULL)}
          infile7<-input$clusters
    })
    
    print(paste("Number Clusters",input_clusters()))
    
    
    ###
    ### Input annotation
    ###
    
    genes_annotated_fv<-renderText({ input$genes_annotated})

    #if you want use a bed file then are required an annotation file,max-gap, bed file
    if(genes_annotated_fv()==TRUE){
      annotation_dataset<-reactive({
        if(is.null(input$annotation_dataset))  {showNotification("You decided to annotate a bed file. Upload an annotation file.",type="error")}
        infile9<-input$annotation_dataset
        read.table(file=infile9$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
      })
      
      print(paste("Load annotation data",dim(annotation_dataset())))
      
    }
    
    if(genes_annotated_fv()==TRUE){
      
      distance_max<-input$distance
      print(paste("This is the distance",distance_max))
      
    }
    
    if(genes_annotated_fv()==TRUE){
      
      bed_dataset<-reactive({
        if(is.null(input$bed_file))  {showNotification("You decided to annotate a bed file. Upload a bed file.",type="error")}
        infile_bed<-input$bed_file
        read.table(file=infile_bed$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
      })   
      print(paste("Upload bed file",dim(bed_dataset())))
    }  
    
    ###
    ### Input Type of analysis
    ###
    input_list_of_genes_test<-renderText({input$list_of_genes})
    
    #if the user want use a custom list, upload it!  
    if(input_list_of_genes_test()==TRUE){
      
      list_genes_for_analysis<-reactive({
      if(is.null(input$list_of_genes2)){showNotification("You selected to run analysis with a list of genes. Upload a list of genes.",type="error")}
        infile_LG<-input$list_of_genes2
        read.table(file=infile_LG$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
        
      })
      
      print(paste("Upload list genes file",dim(list_genes_for_analysis())))
      
    }
    
    
    all_genes_test<-renderText({input$all_genes})
    
    #if the use want use the entire list of genes
    if(all_genes_test()==TRUE) {
      
    if(is.null(all_genes_test()==TRUE)){showNotification("You select to run the analysis with all genes. Set an appropriate number of clusters.",type="warning")}
      
      list_genes_dataset<-unique(c(input_GE()[,1],input_CNV()[,1],input_MUT()[,1],input_METH()[,1]))
      all_genes_for_analysis<-list_genes_dataset
      print("The number of genes in common is:")
      print(length(all_genes_for_analysis))
      
    } 
    
    output_file<-input$output_file

########## END UPLOAD FILE
##########

########## Perform some controls
    
if(!is.null(input_GE()) & !is.null(input_CNV()) & !is.null(input_METH()) & !is.null(input_MUT())){
      
      showNotification("Control matrices: ok!",type="message",duration=NULL)
      
    } else {
      
      showNotification("One of your matrix is empty",type="error")
      showNotification("If you want analyze less than 4 data-sets click on Analysis only two datasets",type="error")
      
}   
    
    

######### CREATE THE GENE LISTS
# use a list from annotation
if(genes_annotated_fv()==TRUE){
  
  showNotification("Pre-processing: selection of the genes from the data",type="message")
  
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
  
  print(length(input_GE_selected))
  print(length(input_CNV_selected))
  print(length(input_METH_selected))
  print(length(input_MUTATION_selected))
  
  }

# use a list of genes
if(input_list_of_genes_test()==TRUE){
  
  showNotification("Pre-processing: selection of the genes from the data",type="message")
  
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
  
  
  print(length(input_GE_selected))
  print(length(input_CNV_selected))
  print(length(input_METH_selected))
  print(length(input_MUTATION_selected))
  
  } 
  
# use all genes
if(all_genes_test() == TRUE){
    
    showNotification("Pre-processing: selection of the genes from the data",type="message",type='warning')
  
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()
    
    input_GE_selected<-input_GE2[input_GE2[,1]%in%all_genes_for_analysis,]
    input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%all_genes_for_analysis,]
    input_METH_selected<-input_METH2[input_METH2[,1]%in%all_genes_for_analysis,]
    input_MUTATION_selected<<-input_MUTATION2[input_MUTATION2[,1]%in%all_genes_for_analysis,] 
    
    print(length(input_GE_selected))
    print(length(input_CNV_selected))
    print(length(input_METH_selected))
    print(length(input_MUTATION_selected))
    
    } 

    
##
## Start the algorithm
##
 
      if(TD() != TRUE){
      
      ##  check
      if(dim(input_GE_selected)[1]!=0 & dim(input_CNV_selected)[1]!=0 & dim(input_METH_selected)[1]!=0 & dim(input_MUTATION_selected)[1]!=0){
        
      showNotification("Control matrices: ok!",type="message")
        
      } else {
      
      showNotification("One of your matrix is empty",type="error",duration=NULL)

      }
      
      } else {
        
        showNotification("Pre-processing of the data",type="message")
        
        if(cb_ge_ok() == FALSE & cb_cnv_ok() == FALSE & cb_meth_ok() == FALSE & cb_mutation_ok() == FALSE){showNotification("Select which omic data set you uploaded.",type="error")}
                    
        if(cb_ge_ok() == TRUE){input_GE_string="GE"} else {input_GE_string=NULL}
        if(cb_cnv_ok() == TRUE){input_CNV_string="CNV"} else{input_CNV_string=NULL}
        if(cb_meth_ok() == TRUE){input_METH_string="METH"} else{input_METH_string=NULL}
        if(cb_mutation_ok() == TRUE){input_MUT_string="MUT"} else{input_MUT_string=NULL}
        
        print(input_GE_string)
        print(input_CNV_string)
        print(input_METH_string)
        print(input_MUT_string)
        
        list_of_experiment_default<-c("GE","CNV","METH","MUT")
        list_of_experiment_user<-c(input_GE_string,input_CNV_string,input_METH_string,input_MUT_string)
        print(list_of_experiment_user)
        
        experiments_not_presents<-setdiff(list_of_experiment_default,list_of_experiment_user)
        experiments_presents<-intersect(list_of_experiment_default,list_of_experiment_user)
        
        print(experiments_not_presents)
        print(experiments_presents)
        
        dictionary<-data.frame(code_exp=list_of_experiment_default,strings_variable=c("input_GE_selected","input_CNV_selected","input_METH_selected","input_MUTATION_selected"),stringsAsFactors=F)
        var_exp_name<-dictionary[dictionary[,1]%in%experiments_presents,2]
        
        
        list_object<-list(input_GE_selected,input_CNV_selected,input_METH_selected,input_MUTATION_selected)
        names(list_object)<-c("input_GE_selected","input_CNV_selected","input_METH_selected","input_MUTATION_selected")
        
        list_object2<-list_object[var_exp_name]

        #control genes and samples availables
        GENES_AVAILABLE<-NULL
        SUBJ_AVAILABLE<-NULL
        
        #save genes and subjects available
        for(i in var_exp_name){
        
        print(i)
          
        get_i<-list_object2[[i]]
        print(dim(get_i))
        
        #if the number of columns is different from 2 are not mutation  
          if(dim(get(i))[2]!=2){
            subjects<-colnames(get(i))
            genes<-get(i)[,1]
            
          } else {
            
        #if the number of columns is equal to 2  not mutation , get the patients from the second columns 
            
            subjects<-get(i)[,2]
            genes<-get(i)[,1]
            
          }

        GENES_AVAILABLE<-c(GENES_AVAILABLE,genes)
        SUBJ_AVAILABLE<-c(SUBJ_AVAILABLE,subjects)
        
        }
        GENES_AVAILABLE2<-unique(GENES_AVAILABLE)
        SUBJ_AVAILABLE2<-unique(SUBJ_AVAILABLE)
        
        #create false data.frame for the analysis
        for(i in experiments_not_presents){
          
          if(i == "GE"){
            input_GE_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
            input_GE_selected[,1]<-GENES_AVAILABLE2
            colnames(input_GE_selected)<-SUBJ_AVAILABLE2
            print(dim(input_GE_selected))
            }
          
          if(i == "CNV"){
            input_CNV_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
            input_CNV_selected[,1]<-GENES_AVAILABLE2
            colnames(input_CNV_selected)<-SUBJ_AVAILABLE2
            print(dim(input_CNV_selected))
            
            }

          if(i=="METH"){
            input_METH_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
            input_METH_selected[,1]<-GENES_AVAILABLE2
            colnames(input_METH_selected)<-SUBJ_AVAILABLE2
            print(dim(input_METH_selected))
            
          }
          
          if(i=="MUT"){
            
            input_MUTATION_selected<-data.frame(matrix(0, ncol = 2, nrow = length(SUBJ_AVAILABLE2)-1))
            input_MUTATION_selected[,2]<-SUBJ_AVAILABLE2[-1]#the first elements is alwasy the id of genes
            input_MUTATION_selected[,1]<-'' #for mutation no genes are avaible
            colnames(input_MUTATION_selected)<-c("genes","sample_id")
          }
          
        }
      }
      
      ################################ RUN GMIEC RFK
      print("run gmiec RFK!")
      
      if(RFA() == TRUE & GMR()==TRUE){showNotification("Select only one method of analysis.",type="error")}
        
      if((RFA() == TRUE & TD() == TRUE )|(RFA() == TRUE & TD() == FALSE)){
      
      showNotification("Analysis started, wait!",type="message",duration=NULL)
        
      
      output_gmiec<-reactive({GMIEC_MLK(
        
        input_GE_selected=input_GE_selected,
        input_CNV_selected=input_CNV_selected,
        input_METH_selected=input_METH_selected,
        input_MUTATION_selected=input_MUTATION_selected,
        drugs_for_analysis2=drugs_for_analysis2,
        input_clinical=input_clinical2,
        k_user=input_clusters(),
        all_genes_for_analysis=list_genes_for_analysis()
        
      )})
      
      if(!is.null(output_gmiec())){
        
        showNotification("You can download the results of analysis!",type="warning")
        print("here your output")
        print(dim(output_gmiec()))  
      
        output$downloadData <- downloadHandler(
          
          filename = function() {
            paste('Analysis_GMIEC_RFK_main_results.', Sys.Date(), '.csv', sep='')
          }
          ,
          content = function(file) {
            write.csv(output_gmiec(),file,row.names=FALSE, na="",sep=";")
          }
        )
        
      }
      
      }
    ################################ RUN GMIEC LOGIC APPROACH
      if((GMR() == TRUE & TD() == TRUE )|(GMR()==TRUE & TD() == FALSE)){
        
        print("run gmiec LRK!")
        
        output_gmiec<-reactive({run_GMIEC(
          
          check_ge_for_patients=input_GE_selected,
          input_CNV_selected=input_CNV_selected,
          input_METH_selected=input_METH_selected,
          input_MUTATION_selected=input_MUTATION_selected,
          tabDrugs=drugs_for_analysis2,
          input_clinical=input_clinical2,
          parameter_discr=c("2;1;0.7"),
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
            
            filename = function() {
              paste('Analysis_GMIEC_main_results.', Sys.Date(), '.csv', sep='')
            }
            ,
            content = function(file) {
              write.csv(output_gmiec(),file,row.names=FALSE, na="",sep=";")
            }
          )
          
        }
        
        
      }
        

}


)

##############GMIEC-VIS
observeEvent(input$run_vis,{
  
  input_for_report2<-reactive({
    infile_for_report<- input$vis_gmiec2
    read.csv(file=infile_for_report$datapath) #read empty values with 0 
  }) 
  
  print(dim(input_for_report2()))
  
  output$list_patients<-renderUI({
      choices=input_for_report2()[,1]
      selectInput('list_patients', 'Select patient',choices, selectize=FALSE)
  })
  
  print(input$list_patients)
  print("Create heatmaps!")
  
  output$plot_heatmap_scores_drugs<-renderPlotly({
    hts<-plot_heatmap_report_gmiec(input_for_report2(),"drugs")
    hts
    }
    )
  output$plot_heatmap_scores_genes<-renderPlotly({
    hts2<-plot_heatmap_report_gmiec(input_for_report2(),"genes")
    hts2
    }
    )
  output$plot_heatmap_scores_sad<-renderPlotly({
    hts3<-plot_heatmap_report_gmiec(input_for_report2(),"sad")
    hts3
    }
    )
  
  output$number_modules<-renderUI({
    choices=paste(1:length(grep(colnames(input_for_report2()),pattern="rdg")),sep="")
    selectInput('number_modules', 'Select module',choices, selectize=FALSE)
    
  }
  )
  
  output$table_summary_scores<-renderText(plotTable(input_for_report2(),input$list_patients))
  output$table_summary_genes_module<-renderText(plot_summary_genes_drugs(input_for_report2(),input$list_patients,type="genes",module=input$number_modules))
  output$table_summary_drugs_module<-renderText(plot_summary_genes_drugs(input_for_report2(),input$list_patients,type="drugs",module=input$number_modules))
  
})

#### GMIEC PARSE RESULTS
observeEvent(input$create_simple_report,{
  
  input_for_report_simple<-reactive({
    if(is.null(input$vis_gmiec3)){showNotification("Upload a file to parse the output of GMIEC.",type="error")}
    infile_for_report<- input$vis_gmiec3
    read.csv(file=infile_for_report$datapath) #read empty values with 0 
  }) 
  
  print(dim(input_for_report_simple()))
  
  showNotification("The output will be created!",type="message")
  
  type_of_algorithm<-renderText({input$type_input_gmiec})
  type_of_output<-renderText({input$simply_res_gmiec})
  
  print(type_of_algorithm())
  print(type_of_output())
  
  export_simple_results_gmiec<-reactive({create_output_gmiec_parse(input_for_report_simple(),type_of_algorithm(),type_of_output())})
  
  print(dim(export_simple_results_gmiec))

if(!is.null(export_simple_results_gmiec())){
  
  showNotification("You can download the results!",type="message")
  print("here your output")
  print(dim(export_simple_results_gmiec()))  
  
  output_gmiec_simple<-export_simple_results_gmiec()
  
  print(dim(export_simple_results_gmiec))
  
  output$download_simple_report <- downloadHandler(
    
    filename = function() {
      paste('GMIEC_selected_module_output', Sys.Date(), '.csv', sep='')
    }
    ,
    content = function(file) {
      write.csv(output_gmiec_simple,file,row.names=FALSE, na="",sep=";")
    }
  )
  
}
  
})

##############GMIEC plot single patient
observeEvent(input$run_plot_res_gmiec,{
  showNotification("Run the analysis, wait the uploading of the data",type="message")
  
  ### read gene expression data
  input_GE_selected<-reactive({
    
    if(is.null(input$ge_dataset_res)) {return(NULL)}
    infile2<- input$ge_dataset_res
    read.table(file=infile2$datapath,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=F,check.names = F) #read empty values with 0 
    
  })
  
  print(paste("Load gene-expression data",length(input_GE_selected())))
  
  
  ### read copy number variation data
  input_CNV_selected<-reactive({
    if(is.null(input$cnv_dataset_res)) {return(NULL)}
    infile3<-input$cnv_dataset_res
    read.table(file=infile3$datapath,sep="\t",stringsAsFactors=F,header=T,quote=NULL,fill=T,check.names = F)
  })
  
  print(paste("Load cnv data",length(input_CNV_selected())))
  
  ### read mutation data
  
  input_MUTATION_selected<-reactive({
    if(is.null(input$mut_dataset_res))  {return(NULL)}
    infile4<-input$mut_dataset_res
    read.table(file=infile4$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
  })    
  
  print(paste("Load mut data",length(input_MUTATION_selected())))
  
  ## read methylation data
  input_METH_selected<-reactive({
    if(is.null(input$meth_dataset_res))  {return(NULL)}
    infile5<-input$meth_dataset_res
    read.table(file=infile5$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
  })
  
  print(paste("Load methylation data",length(input_METH_selected())))
  
  
  ## read  input drugs
  input_DRUGS<-reactive({
    if(is.null(input$genes_drugs))  {showNotification("The drug gene file it is mandatory.",type="error")}
    infile6<-input$genes_drugs
    read.table(file=infile6$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
  })
  
  print(paste("Load drugs data",length(input_DRUGS())))
  
  
  cb_ge_ok2 <- renderText({input$cb_ge2})
  cb_cnv_ok2 <- renderText({input$cb_cnv2})
  cb_meth_ok2 <- renderText({input$cb_meth2})
  cb_mutation_ok2 <- renderText({input$cb_mutation2})
  
  TD2<-renderText({input$two_datasets2})
  
  input_parse_gmiec<-reactive({
  if(is.null(input$results_gmiec_from_parse))  {showNotification("Upload the results of GMIEC.",type="error")}
    infile_for_report<- input$results_gmiec_from_parse
    read.csv(file=infile_for_report$datapath) #read empty values with 0 
  }) 
  
  if(dim(input_parse_gmiec())[2]>4){showNotification("The upload is malformed. Number of columns different from the GMIEC standard.",type="error")}
  if(dim(input_parse_gmiec())[2]>4){showNotification("The upload is malformed. Number of columns different from the GMIEC standard.",type="error")}
  
  if(length(grep(colnames(input_parse_gmiec()),pattern="patientID"))){showNotification("The uploaded file is malformed. patientID column absent, or the name of the colum is wrong",type="error")}
  if(length(grep(colnames(input_parse_gmiec()),pattern="score"))){showNotification("The uploaded file is malformed. score column absent, or the name of the column is wrong",type="error")}
  if(length(grep(colnames(input_parse_gmiec()),pattern="genes_in_module"))){showNotification("The uploaded file is malformed. genes_in_module column absent, or the name of the column is wrong",type="error")}
  if(length(grep(colnames(input_parse_gmiec()),pattern="drugs_in_module"))){showNotification("The uploaded file is malformed. drugs_in_module column absent, or the name of the column is wrong",type="error")}
  
  print(dim(input_parse_gmiec()))
  
  output$list_patients2<-renderUI({
    choices=input_parse_gmiec()[,1]
    selectInput('list_patients2', 'Select patient',choices, selectize=FALSE)
  })
  

  if(TD2() == TRUE){
  
    
    if(cb_ge_ok2() == TRUE){input_GE_string="GE"} else {input_GE_string=NULL}
    if(cb_cnv_ok2() == TRUE){input_CNV_string="CNV"} else{input_CNV_string=NULL}
    if(cb_meth_ok2() == TRUE){input_METH_string="METH"} else{input_METH_string=NULL}
    if(cb_mutation_ok2() == TRUE){input_MUT_string="MUT"} else{input_MUT_string=NULL}
    
    print(input_GE_string)
    print(input_CNV_string)
    print(input_METH_string)
    print(input_MUT_string)
    
    list_of_experiment_default<-c("GE","CNV","METH","MUT")
    list_of_experiment_user<-c(input_GE_string,input_CNV_string,input_METH_string,input_MUT_string)
    print(list_of_experiment_user)
    
    experiments_not_presents<-setdiff(list_of_experiment_default,list_of_experiment_user)
    experiments_presents<-intersect(list_of_experiment_default,list_of_experiment_user)
    
    print(experiments_not_presents)
    print(experiments_presents)
    
    dictionary<-data.frame(code_exp=list_of_experiment_default,strings_variable=c("input_GE_selected","input_CNV_selected","input_METH_selected","input_MUTATION_selected"),stringsAsFactors=F)
    var_exp_name<-dictionary[dictionary[,1]%in%experiments_presents,2]
    
    
    list_object<-list(input_GE_selected(),input_CNV_selected(),input_METH_selected(),input_MUTATION_selected())
    names(list_object)<-c("input_GE_selected","input_CNV_selected","input_METH_selected","input_MUTATION_selected")
    
    list_object2<-list_object[var_exp_name]
    
    #control genes and samples availables
    GENES_AVAILABLE<-NULL
    SUBJ_AVAILABLE<-NULL
    
    print(var_exp_name)
    #save genes and subjects available
    for(i in var_exp_name){
      
      print(i)
      
      print(names(list_object2))
      
      get_i<-list_object2[[i]]
      print("test")
      print(dim(get_i))
      
      #if the number of columns is different from 2 are not mutation  
      if(dim(get_i)[2]!=2){
        
        subjects<-colnames(get_i)
        genes<-get_i[,1]
        
      } else {
        
        #if the number of columns is equal to 2  not mutation , get the patients from the second columns 
        
        subjects<-get_i[,2]
        genes<-get_i[,1]
        
      }
      
      GENES_AVAILABLE<-c(GENES_AVAILABLE,genes)
      SUBJ_AVAILABLE<-c(SUBJ_AVAILABLE,subjects)
      
    }
    GENES_AVAILABLE2<-unique(GENES_AVAILABLE)
    SUBJ_AVAILABLE2<-unique(SUBJ_AVAILABLE)
    
    #create false data.frame for the analysis
    for(i in experiments_not_presents){
      
      if(i == "GE"){
        input_GE_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
        input_GE_selected[,1]<-GENES_AVAILABLE2
        colnames(input_GE_selected)<-SUBJ_AVAILABLE2
        print("create ge")
        print(dim(input_GE_selected))
      }
      
      if(i == "CNV"){
        input_CNV_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
        input_CNV_selected[,1]<-GENES_AVAILABLE2
        colnames(input_CNV_selected)<-SUBJ_AVAILABLE2
        print("create cnv")
        print(dim(input_CNV_selected))
        
      }
      
      if(i=="METH"){
        input_METH_selected<-data.frame(matrix(0, ncol = length(SUBJ_AVAILABLE2), nrow = length(GENES_AVAILABLE2)))
        input_METH_selected[,1]<-GENES_AVAILABLE2
        print("create meth")
        colnames(input_METH_selected)<-SUBJ_AVAILABLE2
        print(dim(input_METH_selected))
        
      }
      
      if(i=="MUT"){
        
        input_MUTATION_selected<-data.frame(matrix(0, ncol = 2, nrow = length(SUBJ_AVAILABLE2)-1))
        input_MUTATION_selected[,2]<-SUBJ_AVAILABLE2[-1]#the first elements is alwasy the id of genes
        input_MUTATION_selected[,1]<-'' #for mutation no genes are avaibl
        colnames(input_MUTATION_selected)<-c("genes","sample_id")
        print("create mut")
        print(dim(input_MUTATION_selected))
        
      }
      
    }
  }
  

  if(class(input_GE_selected)[1]!="data.frame"){input_GE_selected_plot=input_GE_selected()}else{input_GE_selected_plot=input_GE_selected}
  if(class(input_CNV_selected)[1]!="data.frame"){input_CNV_selected_plot=input_CNV_selected()}else{input_CNV_selected_plot=input_CNV_selected}
  if(class(input_METH_selected)[1]!="data.frame"){input_METH_selected_plot=input_METH_selected()}else{input_METH_selected_plot=input_METH_selected}
  if(class(input_MUTATION_selected)[1]!="data.frame"){input_MUTATION_selected_plot=input_MUTATION_selected()}else{input_MUTATION_selected_plot=input_MUTATION_selected}
  
  print(input_MUTATION_selected_plot)
  
  print(dim(input_GE_selected_plot))
  print(dim(input_CNV_selected_plot))
  print(dim(input_METH_selected_plot))
  print(dim(input_MUTATION_selected_plot))
  print(dim(input_DRUGS()))
  
  output$plot_module_patient<-renderPlot({
    showNotification("Wait the results of the analysis",type="message")
    
    plot_heatmap_module(res_gmiec=input_parse_gmiec(),
                        input_GE=input_GE_selected_plot,
                        input_CNV=input_CNV_selected_plot,
                        input_METH=input_METH_selected_plot,
                        input_MUT=input_MUTATION_selected_plot,
                        input_DRUGS=input_DRUGS(),
                        subject=input$list_patients2)
    }
  )
  
  output$table_summary_drugs_module2<-renderText(plot_summary_genes_drugs(input_parse_gmiec(),input$list_patients2,type="drugs",module=input$number_modules))

  output$download_html_report<-downloadHandler(
    
    filename = function() {
      paste('GMIEC_report', Sys.Date(), '.html', sep='')
    }
    ,
    content = function(file) {
      
      dir_src<-paste(paste(getwd(),"/src",sep=""),"template_report.Rmd",sep="/")
      
      showNotification("Wait...creation report!",type="message")
      showNotification("This process is slow!",type="message")
      
      rmarkdown::render(dir_src, params = list(
        res_gmiec = input_parse_gmiec(),
        input_GE=input_GE_selected_plot,
        input_CNV=input_CNV_selected_plot,
        input_METH=input_METH_selected_plot,
        input_MUT=input_MUTATION_selected_plot,
        input_DRUGS=input_DRUGS()))
      
      showNotification("Your report is ready!",type="message")
      
    }
  )
  
}

)
# end code server.R
}
