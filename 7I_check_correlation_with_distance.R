library(raster)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(geosphere)

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

centroids <- read.csv("clustering/centroids_coord.csv")

# Set index coords:
index_coord <- data.frame(index=c("enso","iod","nao","pdo"),
           lon=c(-145, 80,-25,-140),
           lat=c(0,-5, 52, 45))

# Calculate closest distance 
#The shortest distance between two points (i.e., the 'great-circle-distance' or 'as the crow flies'), according to the 'haversine method'. This method assumes a spherical earth, ignoring ellipsoidal effects.

centroids$enso <- NA; centroids$iod <- NA; 
centroids$nao <- NA; centroids$pdo <- NA

centroids$enso_cor <- NA; centroids$iod_cor <- NA; 
centroids$nao_cor <- NA; centroids$pdo_cor <- NA

oscillations <- c("enso", "iod", "nao", "pdo")
for (r in 1:nrow(centroids)){
  for (i in oscillations){
    if (i %in% summary_ccf_pc10per[summary_ccf_pc10per$cluster==centroids$DN[r],]$index){
      #centroids[c(i)][r,] <- "Y"
      
      # Define coordinates (longitude, latitude)
      coord1 <- c(index_coord$lon[index_coord$index==i], 
                  index_coord$lat[index_coord$index==i])
      coord2 <- c(centroids$X[r], centroids$Y[r])
      
      # Calculate distance in meters
      distance <- distHaversine(coord1, coord2)
      centroids[i][r,] <- distance/1000
      
      centroids[paste0(i,"_cor")][r,] <- mean(summary_ccf_pc10per[summary_ccf_pc10per$cluster==centroids$DN[r] & summary_ccf_pc10per$index==i,]$cor)
      
    }
  }
}


cor.test(centroids$enso, centroids$enso_cor, na.rm=T)
cor.test(centroids$iod, centroids$iod_cor, na.rm=T)
cor.test(centroids$nao, centroids$nao_cor, na.rm=T)
cor.test(centroids$pdo, centroids$pdo_cor, na.rm=T)

range(centroids[oscillations],na.rm=T)
