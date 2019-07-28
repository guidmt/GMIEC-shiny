options(shiny.maxRequestSize=100*1024^2)

source("./src/create_output.R", local = FALSE)
source("./src/engine_all_dataset.R", local = FALSE)
source("./src/filter_ge.R", local = FALSE)
source("./src/internal_annotation.R", local = FALSE)
source("./src/rules_fornot_tf.R", local = FALSE)
source("./src/GMIEC.R", local = FALSE)
source("./src/gmiec_ML.R",local=FALSE)
source("./src/plot_shiny_functions.R",local=FALSE)

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
        
      if(is.null(input$ge_dataset)) {return(NULL)}
      infile2<- input$ge_dataset
      read.table(file=infile2$datapath,sep="\t",stringsAsFactors=FALSE,header=T,quote=NULL,fill=F,check.names = F) #read empty values with 0 
      
      })
     
    print(paste("Load gene-expression data",length(input_GE())))
    

    ### read copy number variation data
      input_CNV<-reactive({
      if(is.null(input$cnv_dataset)) {return(NULL)}
      infile3<-input$cnv_dataset
      read.table(file=infile3$datapath,sep="\t",stringsAsFactors=F,header=T,quote=NULL,fill=T,check.names = F)
      })
    
    print(paste("Load cnv data",length(input_CNV())))
    
    ### read mutation data
  
      input_MUT<-reactive({
      if(is.null(input$mut_dataset))  {return(NULL)}
      infile4<-input$mut_dataset
      read.table(file=infile4$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
      })    
    
    print(paste("Load mut data",length(input_MUT())))
    
    ## read methylation data
    input_METH<-reactive({
      if(is.null(input$meth_dataset))  {return(NULL)}
      infile5<-input$meth_dataset
      read.table(file=infile5$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
    })

    print(paste("Load meth data",length(input_METH())))
    
    
    ## read clinical data
    input_clinical<-reactive({
    if(is.null(input$clinical_dataset))  {return(NULL)}
      infile6<-input$clinical_dataset
      read.table(file=infile6$datapath,sep="\t",stringsAsFactors=F,quote=NULL,header=T,fill=T,check.names = F)
    })
    
    print(paste("Load clinical data",dim(input_clinical())))
    
  
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
        infile_bed<-input$bed_file
        read.table(file=infile_bed$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
      })   
      print(paste("Upload bed file",dim(bed_dataset())))
    }  
    
    ###
    ### Input Drugs
    ###
    
    drugs_for_analysis<-reactive({
      infile1<-input$drugs_dataset
      read.table(file=infile1$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
    })
    
    print(paste("Upload drugs file",dim(drugs_for_analysis())))
    
    ###
    ### Input Type of analysis
    ###
    input_list_of_genes_test<-renderText({input$list_of_genes})
    
    #if the user want use a custom list, upload it!  
    if(input_list_of_genes_test()==TRUE){
      
      list_genes_for_analysis<-reactive({
        infile_LG<-input$list_of_genes2
        read.table(file=infile_LG$datapath,sep="\t",stringsAsFactors=F,header=T,fill=T,check.names = F)
        
      })
      
      print(paste("Upload list genes file",dim(list_genes_for_analysis())))
      
    }
    
    
    genes_annotated_TF2_test<-renderText({ input$genes_annotated_TF})
    
    #if the use want evalute the impact of TF
    if(genes_annotated_TF2_test()==TRUE & RFA()!= TRUE) {
      
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
      print("The number of genes in common is:")
      print(length(all_genes_for_analysis))
      
    } 
    
    output_file<-input$output_file

########## END UPLOAD FILE
##########


######### CREATE THE GENE LISTS
    
# use a list from annotation
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
  
  input_GE_selected<<-input_GE2[input_GE2[,1]%in%list_genes_for_analysis,]
  input_CNV_selected<<-input_CNV2[input_CNV2[,1]%in%list_genes_for_analysis,]
  input_METH_selected<<-input_METH2[input_METH2[,1]%in%list_genes_for_analysis,]
  input_MUTATION_selected<<-input_MUTATION2[input_MUTATION2[,1]%in%list_genes_for_analysis,]
  
  print(length(input_GE_selected))
  print(length(input_CNV_selected))
  print(length(input_METH_selected))
  print(length(input_MUTATION_selected))
  
  }

# use a list of genes
if(input_list_of_genes_test()==TRUE){
  
  subsetAnnoDF_unique<-list_genes_for_analysis()[,1]
  
  input_GE2<-input_GE()
  input_CNV2<-input_CNV()
  input_METH2<-input_METH()
  input_MUTATION2<-input_MUT()
  input_clinical2<-input_clinical()
  drugs_for_analysis2<-drugs_for_analysis()
  
  input_GE_selected<<-input_GE2[input_GE2[,1]%in%subsetAnnoDF_unique,]
  input_CNV_selected<<-input_CNV2[input_CNV2[,1]%in%subsetAnnoDF_unique,]
  input_METH_selected<<-input_METH2[input_METH2[,1]%in%subsetAnnoDF_unique,]
  input_MUTATION_selected<<-input_MUTATION2[input_MUTATION2[,1]%in%subsetAnnoDF_unique,]
  
  
  print(length(input_GE_selected))
  print(length(input_CNV_selected))
  print(length(input_METH_selected))
  print(length(input_MUTATION_selected))
  
  } 
  
# use all genes
if(all_genes_test() == TRUE){
    
    input_GE2<-input_GE()
    input_CNV2<-input_CNV()
    input_METH2<-input_METH()
    input_MUTATION2<-input_MUT()
    input_clinical2<-input_clinical()
    drugs_for_analysis2<-drugs_for_analysis()
    
    input_GE_selected<<-input_GE2[input_GE2[,1]%in%all_genes_for_analysis,]
    input_CNV_selected<<-input_CNV2[input_CNV2[,1]%in%all_genes_for_analysis,]
    input_METH_selected<<-input_METH2[input_METH2[,1]%in%all_genes_for_analysis,]
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
        
      input_GE2<-input_GE()
      input_CNV2<-input_CNV()
      input_METH2<-input_METH()
      input_MUTATION2<-input_MUT()
      

      } else {

                
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
      
      if((RFA() == TRUE & TD() == TRUE )|(RFA() == TRUE & TD() == FALSE)){
      
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
        
        showNotification("You can download the results of analysis!",type="message")
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
  
  output$plot_single_summary_single_patient<-renderText(plotTable(input_for_report2(),input$list_patients))
  
})

# end code server.R
}
