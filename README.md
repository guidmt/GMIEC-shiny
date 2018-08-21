### NOTE: GMIEC-shiny is again in beta version, some step of analysis are again to improve. In the next days some updates will be available:

# News:
21/08/2018:
- update server.R files for the analysis considering the rules without TF
- created script rules_fornot_tf_fd.R
- update run_GMIEC_fd.R

20/08/2018: Start implementation FD-GMIEC
- modify ui.R to implement FD GMIEC analysis
- update index.Rmd in GMIEC_www
- add function rules_for_tf_fd.R and update it. First implementation. V.beta
- add function run_GMIEC_fd.R
- update ui.R

08/09/2018: Big updates before of summer! The code of VIS-GMIEC is complete!

- index_html. Rmd: i) add links to pages created from template_single_patients.Rmd ii) add PCA analysis
- template_single_patients.Rmd: add links to output created in template_single_patients2.Rmd
- template_single_patients2.Rmd: resolved some bugs in the chunk of codes

02/08/2018

- update server.R
- update ui.R
- create function create_report.R
- debug .Rmd files of index, global single patients, reports html templates

01/08/2018

- update functions to run analysis
- update server.R for vis-gmiec, implemented upload and reading of output GMIEC
- update create_report.R and function run_create_output

19/07/2018

- removed apriori analysis and implemented kmodes from klaR 
- modify engine_all_dataset.R  
- modify create_output.R 
- add some comments 

# To do list:
20/08/2018
- debug rules_for_tf_fd.R for the analysis of TF
- implement analysis to not consider TF
- debug create_output.R

08/09/2018:
- implement FD-GMIEC
- improve code VIS-GMIEC (create output .zip with the results), add more details about the report
- improve code GMIEC
- improve documentation

02/09/2018
- update introduction, manual page, removing the description of apriori, insert k-modes
- implement VIS-GMIEC
  - i) complete the implementation
  - ii) optimize code to create the reports
  - iii) create output .zip with the results
- implement FD-GMIEC (handle multiple datasets) 

01/09/2018
- update introduction, manual page, removing the description of apriori, insert k-modes
- implement VIS-GMIEC (code ready, only to implement in gui-shiny)
- implement FD-GMIEC (handle multiple datasets)

05/07/2018

- implementation of FD-GMIEC
- update of the gui interface to upload the data (VIS-GMIEC)

### GMIEC-shiny
This is the offical repository of GMIEC-shiny
Welcome in GMIEC!! GMIEC is a framework, GMIEC, that integrates omics data (gene-expression, copy number variation, mutation, methylation) with data containing genomic coordinates of cis and/or trans epigenetic regulators (e.g. ChIP-seq). Then, at the level of the single patient, GMIEC identify groups of genes (modules) that share common genomic features. Then, using an external database, the genes in each module are associated with their own target drugs. Therefore, GMIEC allows the identification, at the level of the single patient, groups of genes that are regulated by the same epigenetic regulators (e.g. TFs), sharing common genomic features and associated with different number of drugs.

# Methods:
This shiny application contatin different two modules: GMIEC and FD-GMIEC. The input of GMIEC are: i) a file with the genomic regions of the regulators of interest; ii) omics data; ii) a file that includes the association between genes and drugs. The four principal steps of analysis are: i) the annotation of the epigenetic components (ECs) derived from a catalogue (ReMap) with their target genes; ii) a step of discretization used to mark the presence of one alteration (1) or not (0) for all genes of single patient, considering all datasets (e.g. gene-expression data, copy-number alterations data, mutations data); iii) the genes with the similar genomic features grouped using k-mode clustering.

# Functionalities of GMIEC:

Identification of groups of genes with/without genomic alterations associated with drugs using gene-expression, copy-number alteration, methylation, mutation, clinical data.
Analysis using only the annotated genes.
Analysis using only the annotated genes and TF.
Analysis using a custom list of genes.
Analysis using all genes between the datasets.

# Functionalities of FD-GMIEC:
Identification of groups of genes with/without genomic alterations associated with drugs using predefined number of datasets.
Analysis using only the annotated genes.
Analysis using only the annotated genes and TF.
Analysis using a custom list of genes.
Analysis using all genes between the datasets.

# Functionalities of VIS-GMIEC:
Exploration of the data at the level of groups of or single patient


The manual is available in [here](https://cdn.rawgit.com/guidmt/GMIEC-shiny/cf5fffbd/GMIEC_www/manual.html)
