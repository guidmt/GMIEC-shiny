library(shiny)
setwd("C:/Users/guida/Desktop/work/GMIEC_shiny")
options(shiny.fullstacktrace=TRUE)
shinyAppDir(appDir="C:/Users/guida/Desktop/work/GMIEC_shiny")


##
## for fast debug
##

source("./src/create_output.R", local = TRUE)
source("./src/engine_all_dataset.R", local = TRUE)
source("./src/filter_ge.R", local = TRUE)
source("./src/internal_annotation.R", local = TRUE)
source("./src/rules_for_tf.R", local = TRUE)
source("./src/rules_fornot_tf.R", local = TRUE)
source("./src/GMIEC.R", local = TRUE)
setwd("C:/Users/guida/Desktop/stratification_tcga/example_data_gmiec")
input_GE2<-read.delim(file="data_expression_median.txt")
input_CNV2<-read.delim(file="data_CNA.txt")
input_METH2<-read.delim(file="data_methylation_hm27.txt")
input_MUTATION2<-read.delim(file="data_mutations_extended.txt")
colnames(input_GE2)[1]<-"genesID"
colnames(input_CNV2)[1]<-"genesID"
colnames(input_METH2)[1]<-"genesID"
colnames(input_MUTATION2)[1]<-"genesID"
tabDrugs<-read.delim(file="interactions.tsv")
colnames(tabDrugs)[1]<-"genes"
bed_dataset<-read.delim(file="BCL6_nr.bed.gz")
annotation_dataset<-read.delim(file="RefSeq.Hg19.parse.bed")

list_genes_for_analysis<-internal_annotation(bed_dataset,annotation_dataset,20000)

input_GE_selected<-input_GE2[input_GE2[,1]%in%list_genes_for_analysis,]
input_CNV_selected<-input_CNV2[input_CNV2[,1]%in%list_genes_for_analysis,]
input_METH_selected<-input_METH2[input_METH2[,1]%in%list_genes_for_analysis,]
input_MUTATION_selected<-input_MUTATION2[input_MUTATION2[,1]%in%list_genes_for_analysis,]

setwd("C:/Users/guida/Desktop/stratification_tcga/example_data_gmiec")

input_clinical<-read.delim(file="data_clinical.txt",header=T)
colnames(input_clinical)<-"SAMPLE_ID"
parameter_discr=c("1.5;1;0.5")

check_ge_for_patients<-input_GE_selected
