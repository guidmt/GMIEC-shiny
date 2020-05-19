install.packages(c("shiny",
                   "shinydashboard",
                   "klaR",
                   "formattable",
                   "bcp",
                   "plyr",
                   "kableExtra",
                   "randomForest",
                   "heatmaply",
                   "NbClust"
                   ),depencies=T)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
BiocManager::install("ComplexHeatmap")