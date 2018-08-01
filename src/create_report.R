run_create_output<-function(input_for_report){

  rmarkdown::render_site("../VIS-GMIEC/index_html.Rmd")
  rmarkdown::render_site("../VIS-GMIEC/template_single_patient.Rmd")
  rmarkdown::render_site("../VIS-GMIEC/template_single_patient2.Rmd")
  
}