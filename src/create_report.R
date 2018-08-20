run_create_output<-function(input_for_report){
  print("I am inside run_create_output")
  dim(input_for_report)
  print(getwd())
  withProgress(message="Start analysis!",min=0,max=1,{
  incProgress(0.15, detail = "Step1: Create Index Html")
    
  ##
  ## Create report Index
  ##
    
  rmarkdown::render(paste(getwd(),"VIS-GMIEC/index_html.Rmd",sep="/"),params= list(input_for_gmiec_report=input_for_report))
  
  ##
  ## Create report Single patient 
  ##
  
  incProgress(0.15, detail = "Step2: Create reports single patient html - Step1")
  
  number_of_patients<-2:ncol(input_for_report)
  
  for(i in number_of_patients){
    
    current_data_patient<-input_for_report[,c(1,i)]
    sampleID<-colnames(current_data_patient)[2]
    
    output_file<-paste("report",sampleID,"html",sep=".")
    gc()
    rmarkdown::render(paste(getwd(),"VIS-GMIEC/template_single_patient.Rmd",sep="/"),
                      params = list(
                        patient=current_data_patient, name_patient=sampleID
                      ),output_file=output_file,output_dir=paste(getwd(),"VIS-GMIEC/output_vis",sep="/")
    )
    gc()
    rm(current_data_patient)
    rm(sampleID)
  }
  
  ##
  ## Create report Single patient 2
  ##
  incProgress(0.15, detail = "Step3: Create reports single patient html - Step2")
  
  number_of_patients<-2:ncol(input_for_report)
  
  #this is for sad score
  for(i in number_of_patients){
    
    current_data_patient<-input_for_report[,c(1,i)]
    sampleID<-colnames(current_data_patient)[2]
    print(sampleID)
    output_file2<-paste("report",sampleID,"single.patient.html",sep=".")
    
    rmarkdown::render(paste(getwd(),"VIS-GMIEC/template_single_patient2.Rmd",sep="/"),
                      params = list(
                        patient=current_data_patient
                      ),output_file=output_file2,output_dir=paste(getwd(),"VIS-GMIEC/output_vis",sep="/"))
    
    gc()
    rm(current_data_patient)
    rm(sampleID)
    
  }
  ##
  ## End
  ##
  
  incProgress(0.15, detail = "End! Click Download Button")
  
  })
  
}