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

ui <- dashboardPage(
  dashboardHeader(title = "gmiec-app"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("GMIEC-MD", tabName = "run_gmiec", icon = icon("play-circle")),
      menuItem("GMIEC-VIS", tabName = "vis_gmiec2", icon = icon("wpexplorer")),
      menuItem("Manual", tabName = "manual", icon = icon("book")),
      menuItem("FAQ", tabName = "faq", icon = icon("question")),
      menuItem("Contact", tabName = "contact", icon = icon("at")),
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
              fluidRow(
              fileInput("vis_gmiec2", "Upload results GMIEC2",buttonLabel=icon("folder"))),
              actionButton('run_vis', 'Create report!',style = "color: white; 
                     background-color: #0066CC; 
                               position: relative; 
                               left: 3%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px"),
      downloadButton('download_vis', 'Download',style = "color: white; 
                     background-color: #ec0000; 
                               position: relative; 
                               left:10%;
                               height: 35px;
                               width: 200px;
                               text-align:center;
                               text-indent: -2px;
                               border-radius: 6px;
                               border-width: 2px")),
      tabItem(tabName = "manual",
              includeHTML("./GMIEC_www/manual.html")),
      tabItem(tabName = "faq",
              includeHTML("./GMIEC_www/faq.html"))
      ,
      tabItem(tabName = "contact",
              includeHTML("./GMIEC_www/contact.html")),
      tabItem(tabName = "terms_of_use",
              includeHTML("./GMIEC_www/terms_of_use.html"))

    )

  )
)