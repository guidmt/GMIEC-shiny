library(bcp)
library(ChIPpeakAnno)
library(plyr)
library(stats)
library(shiny)
library(shinydashboard)
library(heatmaply)
library(klaR)
library(formattable)
library(kableExtra)
library(randomForest)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)

ui <- dashboardPage(
  dashboardHeader(title = "gmiec-app"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("GMIEC-MD", tabName = "run_gmiec", icon = icon("play-circle")),
      menuItem("GMIEC-VIS", tabName = "vis_gmiec2", icon = icon("wpexplorer")),
      menuItem("GMIEC-results",tabName="gmiec_results",icon=icon("flag")),
      menuItem("Terms of use", tabName = "terms_of_use", icon = icon("user-tie"))
    )
  ),
  
  dashboardBody(
    #tab content
    tabItems(
      
      tabItem(tabName = "home",
              includeHTML("./GMIEC_www/index.html")),
      
      tabItem(tabName = "run_gmiec",
            
                fluidRow(
                  box(title="GMIEC - Input dataset",status="primary",solidHeader=TRUE,collapsible =TRUE,
                    fileInput("ge_dataset", "Upload gene-expression data",buttonLabel=icon("folder-open")),
                    fileInput("cnv_dataset", "Upload copy-number variation data",buttonLabel=icon("folder-open")),
                    fileInput("meth_dataset", "Upload methylation data",buttonLabel=icon("folder-open")),
                    fileInput("mut_dataset", "Upload mutation data",buttonLabel=icon("folder-open")),
                    fileInput("clinical_dataset", "Upload clinical data",buttonLabel=icon("folder-open"))
                         ),#end box
                  box(title="GMIEC - Input annotation",status="warning",solidHeader=TRUE,collapsible =TRUE,
                      fileInput("annotation_dataset", "Upload annotation data",buttonLabel=icon("folder")),
                      numericInput("distance","distance (bp)",value=20000)
                      ),
                  box(title="GMIEC - Input drugs",status="info",solidHeader=TRUE,collapsible =TRUE,
                      fileInput("drugs_dataset", "Upload drugs-genes data",buttonLabel=icon("folder-open"))
                  ),
                  box(title="GMIEC - Type of Analysis",status="danger",solidHeader=TRUE,collapsible =TRUE,
                      fileInput("bed_file", "Upload the bed file - for 1,2",buttonLabel=icon("folder")),
                      checkboxInput("genes_annotated","1) Use only the genes annotated",FALSE),
                      checkboxInput("all_genes","2) Use all genes",FALSE),
                      checkboxInput("list_of_genes","3) Use a list of genes",FALSE),
                      fileInput("list_of_genes2", "Upload the list of genes - only for 4",buttonLabel=icon("folder"))
                      
                  ),
                  box(title="GMIEC - Parameters Analysis",status="success",solidHeader=TRUE,collapsible =TRUE,
                      numericInput("clusters","Number Clusters k-mode/k-means",value=4),
                      checkboxInput("GMIEC_RULES","Analysis with logic rules + k-mode",FALSE),
                      checkboxInput("RF_ANALYSIS","Analysis with randomForest + k-means",TRUE),
                      checkboxInput("two_datasets","Analysis only two datasets (Select)",FALSE),
                      checkboxInput("cb_ge", label = "gene-expression", value = FALSE),
                      checkboxInput("cb_cnv", label = "copy-number", value = FALSE),
                      checkboxInput("cb_meth", label = "methylation", value = FALSE),
                      checkboxInput("cb_mutation", label = "mutation", value = FALSE)
                      
                      ),
                  
                  actionButton('run_gmiec', 'Run analysis',style = "color: white; 
                     background-color: #0066CC; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px"),
                  downloadButton('downloadData', 'Download',style = "color: white; 
                     background-color: #ec0000; 
                               position: relative; 
                               left:10%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px"))
                  
                #end fluid row
              
              ), #end tab name
      
      
      
      tabItem(tabName = "vis_gmiec2",
              tabsetPanel(
              tabPanel("Summary heatmaps GMIEC", 
              fluidRow(
              box(title="Upload results GMIEC",status="success",width=3,solidHeader=TRUE,collapsible =FALSE,  
              fileInput("vis_gmiec2", "Upload results GMIEC",buttonLabel=icon("folder")),
              actionButton('run_vis', 'Create report!',style = "color: white; 
                     background-color: #0066CC; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px")
              ),
              box(title="Heatmap scores drugs",status="warning",solidHeader=TRUE,collapsible =FALSE,
                  plotlyOutput("plot_heatmap_scores_drugs")),
              box(title="Heatmap scores genes",status="warning",solidHeader=TRUE,collapsible =FALSE,
                  plotlyOutput("plot_heatmap_scores_genes")),
              box(title="Heatmap scores combined",status="warning",solidHeader=TRUE,collapsible =FALSE,
                  plotlyOutput("plot_heatmap_scores_sad"))
              )#end fluid row
              ), # end tab panel
      tabPanel("Table summary scores, genes, drugs for patient",
               fluidRow(
               box(title="Select a patient from the list",status="success",width=3,solidHeader=TRUE,collapsible =TRUE,
               uiOutput('list_patients')),
               box(title="Select a module of the selected patient",status="success",width=3,solidHeader=TRUE,collapsible =TRUE,
               uiOutput('number_modules')),
               box(title="Table summary scores for patients",status="warning",solidHeader=TRUE,collapsible =FALSE,
               htmlOutput("table_summary_scores")),
               box(title="Table summary genes module patient",status="warning",solidHeader=TRUE,collapsible =FALSE,
               htmlOutput("table_summary_genes_module")),
               box(title="Table summary drug module patient",status="warning",solidHeader=TRUE,collapsible =FALSE,
               htmlOutput("table_summary_drugs_module"))
               ))
      
      )#end tabset panel
      ),#end tabItem
      tabItem(tabName = "gmiec_results",
              tabsetPanel(
                tabPanel("Parse GMIEC results", 
                         fluidRow(
                           box(title="Upload results GMIEC",status="success",width=3,solidHeader=TRUE,collapsible =FALSE,  
                               fileInput("vis_gmiec3", "Upload results GMIEC",buttonLabel=icon("folder")),
                               selectInput("type_input_gmiec",
                                           "Select which algorithm was used to create the GMIEC-output",
                                           c("Random forest + k-means","Logic rules + k-means")
                                           ),
                               selectInput("simply_res_gmiec", "Choose the type of output:",c("Module active",
                                                                        "Module inactive",
                                                                        "Module active with drugs",
                                                                        "Module inactive with drugs")),
                               actionButton('create_simple_report', 'Create report!',style = "color: white; 
                     background-color: #0066CC; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px"),
                               downloadButton('download_simple_report', 'Download',style = "color: white; 
                     background-color: #ec0000; 
                               position: relative; 
                               left:10%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px")
                           )
                         )
                ),
                tabPanel("Heatmap module for patient",
                         
                         box(title="GMIEC - Input dataset",status="primary",solidHeader=TRUE,collapsible =TRUE,
                             fileInput("results_gmiec_from_parse", "Upload results GMIEC",buttonLabel=icon("folder-open")),
                             fileInput("ge_dataset_res", "Upload gene-expression data",buttonLabel=icon("folder-open")),
                             fileInput("cnv_dataset_res", "Upload copy-number variation data",buttonLabel=icon("folder-open")),
                             fileInput("meth_dataset_res", "Upload methylation data",buttonLabel=icon("folder-open")),
                             fileInput("mut_dataset_res", "Upload mutation data",buttonLabel=icon("folder-open")),
                             fileInput("genes_drugs", "Upload the genes-drugs file",buttonLabel=icon("folder-open")),
                             checkboxInput("two_datasets2","Analysis only two datasets (Select)",FALSE),
                             checkboxInput("cb_ge2", label = "gene-expression", value = FALSE),
                             checkboxInput("cb_cnv2", label = "copy-number", value = FALSE),
                             checkboxInput("cb_meth2", label = "methylation", value = FALSE),
                             checkboxInput("cb_mutation2", label = "mutation", value = FALSE),
                             actionButton('run_plot_res_gmiec', 'Create report!',style = "color: white; 
                     background-color: #0066CC; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px"),
                             downloadButton('download_html_report', 'Download report pdf!',style = "color: white; 
                     background-color: #ec0000; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px")
                         ),
                         
                         box(title="Select patient",status="warning",solidHeader=TRUE,collapsible =FALSE,
                             uiOutput('list_patients2')),
                         box(title="Heatmap module patient",status="warning",solidHeader=TRUE,collapsible =FALSE,
                             plotOutput("plot_module_patient")),
                         box(title="Drugs in current module",status="warning",solidHeader=TRUE,collapsible =FALSE,
                             htmlOutput("table_summary_drugs_module2"))
                         
                         )
              )
              ),
      
      
      tabItem(tabName = "terms_of_use",
              includeHTML("./GMIEC_www/terms_of_use.html"))

    )

  )
)