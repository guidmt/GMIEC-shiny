install.packages(c("shiny",
                   "shinydashboard",
                   "klaR",
                   "formattable",
                   "bcp",
                   "plyr",
                   "kableExtra",
                   "randomForest",
                   "heatmaply"),depencies=T)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")