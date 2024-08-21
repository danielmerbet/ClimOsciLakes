library(raster)
library(ggplot2)
library(sf)
library(RColorBrewer)

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
k <- 47 #number of total clusters
directory  <- paste0(directory_base, "data/",variable)
setwd(directory)

#summary results
summary_ccf <- read.csv("summary_ccf_global.csv")
#add pc variance explained to summary
pc_var <- read.csv("pc_data/pc_proportion_variance_global.csv")
summary_ccf$pc_variance <- NA
for (row in 1:k){
  for (pc in 1:3){
    summary_ccf$pc_variance[summary_ccf$cluster==row & summary_ccf$pc==pc] <- pc_var[row, (pc+1)]
  }
}
#pc variance greater than 10 and remove non-used indexes
summary_ccf_pc10per <- summary_ccf[summary_ccf$pc_variance>=0.10,]
summary_ccf_pc10per <- summary_ccf_pc10per[summary_ccf_pc10per$index!="npgo",] # remove index
summary_ccf_pc10per <- summary_ccf_pc10per[summary_ccf_pc10per$index!="npo",] #remove index
#At NOAA, the official ENSO indicator is the Oceanic Niño Index (ONI)
summary_ccf_pc10per$index[summary_ccf_pc10per$index=="oni"] <- "enso"

# Total number of relations:
nrow(summary_ccf_pc10per)
#bottemp: [1] 152
#surftemp: [1] 174

# Most recurrent index in relations
table(summary_ccf_pc10per$index)/sum(table(summary_ccf_pc10per$index))
#bottemp:
#enso       iod       nao       pdo 
#0.2763158 0.1776316 0.3026316 0.2434211 
#surftemp:
#enso       iod       nao       pdo 
#0.2873563 0.1781609 0.3103448 0.2241379 

