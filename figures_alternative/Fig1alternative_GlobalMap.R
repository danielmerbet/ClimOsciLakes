
# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
k <- 47 #number of total clusters
directory  <- paste0(directory_base, "data/",variable)
setwd(directory)

#temporal while running data again
#rl <- length(list.files(path="ccf/"))
summary_ccf <- read.csv("summary_ccf_global.csv")
#summary_ccf <- summary_ccf[1:(nrow(summary_ccf)/6),]
#summary_ccf <- summary_ccf[,c("index","cluster","pc","cor")]
#temporal while running data again

#summary for significant correlation between PCs of clusters and climate oscillation (index)
summ_04 <- summary_ccf[summary_ccf$cor>=0.4,]

#add lag and and maximum correlation
lag_cor_04 <- read.csv("pca_04/ccf_04_global.csv")
lag_cor_04$lag <- -5:5 
c_cor <- c(); c_lag <- c()
for (r in 1:(ncol(lag_cor_04)-1)){
  pos_max <- which(abs(lag_cor_04[,r])==max(abs(lag_cor_04[,r])))
  c_cor <- c(c_cor, lag_cor_04[pos_max,r])
  c_lag <- c(c_lag, lag_cor_04[pos_max,(nrow(summ_04)+1)])
}
summ_04$cor <- c_cor
summ_04$lag <- c_lag

#PC proportion  of variance
pc_var <- read.csv("pc_data/pc_proportion_variance_global.csv")
pc_var$cluster <- 1:nrow(pc_var)
prop_var <- c()
for (r in 1:nrow(summ_04)){
  summ_04$cluster[r]
  r_row <- which(pc_var$cluster==summ_04$cluster[r])
  r_col <- (1+summ_04$pc[r])
  
  prop_var <- c(prop_var , pc_var[r_row, r_col])
}
summ_04$prop_var_pc <- prop_var

clustering <- raster(paste0(directory,"/clustering/global.tiff"))

#add pc variance explained to summary
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

centroid_coords <-read.csv("clustering/centroids_coord.csv")
summary_ccf_pc10per$lon <- NA
summary_ccf_pc10per$lat <- NA
for (k_temp in 1:k){
  if (k_temp!=16){ #outside the tropical cluster number 16, it has 3 regions in different longitudes
    summary_ccf_pc10per$lon[summary_ccf_pc10per$cluster==k_temp] <- centroid_coords$X[centroid_coords$DN==k_temp]
    summary_ccf_pc10per$lat[summary_ccf_pc10per$cluster==k_temp] <- centroid_coords$Y[centroid_coords$DN==k_temp]
  }else{
    summary_ccf_pc10per$lon[summary_ccf_pc10per$cluster==k_temp] <- centroid_coords$X[centroid_coords$DN==k_temp][1]
    summary_ccf_pc10per$lat[summary_ccf_pc10per$cluster==k_temp] <- centroid_coords$Y[centroid_coords$DN==k_temp][1]
  }
  
}

add_summ <- summary_ccf_pc10per[summary_ccf_pc10per$cluster==16,]
#add second coords for tropical cluster
add_summ_temp <- add_summ
add_summ_temp$lon <- centroid_coords$X[centroid_coords$DN==16][2]
add_summ_temp$lat <- centroid_coords$Y[centroid_coords$DN==16][2]
summary_ccf_pc10per <- rbind(summary_ccf_pc10per, add_summ_temp)
#add third coords for tropical cluster
add_summ_temp <- add_summ
add_summ_temp$lon <- centroid_coords$X[centroid_coords$DN==16][3]
add_summ_temp$lat <- centroid_coords$Y[centroid_coords$DN==16][3]
summary_ccf_pc10per <- rbind(summary_ccf_pc10per, add_summ_temp)

#Plot map relating the climate oscillation with each cluster
library(ggplot2)
library(ggmap)
library(ncdf4)
library(raster)
library(dplyr)

summary_ccf_pc10per$lon_index <- NA
summary_ccf_pc10per$lat_index <- NA

#NAO location
summary_ccf_pc10per[summary_ccf_pc10per$index=="nao",]$lon_index <- -25
summary_ccf_pc10per[summary_ccf_pc10per$index=="nao",]$lat_index <- 52

#ENSO location for  clusters in America
#summary_ccf_pc10per$lon_index[(summary_ccf_pc10per$cluster %in% c(1,2,4:22)) & (summary_ccf_pc10per$index=="oni")] <- -145
#summary_ccf_pc10per$lat_index[(summary_ccf_pc10per$cluster %in% c(1,2,4:22)) & (summary_ccf_pc10per$index=="oni")] <- 0
#ENSO location for  clusters in Africa and Asia
#summary_ccf_pc10per$lon_index[(summary_ccf_pc10per$cluster %in% c(3,23:47)) & (summary_ccf_pc10per$index=="oni")] <- 160
#summary_ccf_pc10per$lat_index[(summary_ccf_pc10per$cluster %in% c(3,23:47)) & (summary_ccf_pc10per$index=="oni")] <- 0
summary_ccf_pc10per[summary_ccf_pc10per$index=="oni",]$lon_index <- -145
summary_ccf_pc10per[summary_ccf_pc10per$index=="oni",]$lat_index <- 0

#PDO location
summary_ccf_pc10per[summary_ccf_pc10per$index=="pdo",]$lon_index <- -140
summary_ccf_pc10per[summary_ccf_pc10per$index=="pdo",]$lat_index <- 45

#IOD location
summary_ccf_pc10per[summary_ccf_pc10per$index=="iod",]$lon_index <- 80
summary_ccf_pc10per[summary_ccf_pc10per$index=="iod",]$lat_index <- -5

#Load clustering map to set the pc_variance
clustering <- raster(paste0(directory,"/clustering/global.tiff"))

index_names <- c("ENSO", "NAO", "PDO", "IOD")
c <- 0
for (i in c("oni", "nao", "pdo", "iod")){
  c <- c+1
  data_index <- summary_ccf_pc10per[summary_ccf_pc10per$index==i, ] 
  
  arrow_data <- data.frame(
    start_lon = data_index$lon_index, # Longitudes of start points
    start_lat = data_index$lat_index,  # Latitudes of start points
    end_lon = data_index$lon,  # Longitudes of end points
    end_lat = data_index$lat,   # Latitudes of end points
    color = data_index$index, #Color of the arrow for each index  
    width = (data_index$cor) #relate the size of the arrow with the correlation with the indexes
  )
  
  # Define index as text for the map
  #text_data <- data.frame(
  #  lon = c(-25, -145, 160, -140, 80),  # Longitudes where you want the text
  #  lat = c(52, 0, 0, 45, -5),      # Latitudes where you want the text
  #  label = c("NAO", "ENSO", "ENSO", "PDO", "IOD")  # Text labels
  #)
  
  
  #for the repeated correlation in cluster 16
  add_after <- data_index[(nrow(data_index)-1):nrow(data_index),]
  data_index <- data_index[1:(nrow(data_index)-2),]
  
  # Aggregate by cluster using sum function
  aggregated_data <- data_index %>%
    group_by(cluster) %>%
    summarize(
      pc_variance = sum(pc_variance)
    )
  
  aggregated_data <- rbind(aggregated_data, add_after[c("cluster", "pc_variance")])
  
  
  #Set map
  pcvariance_map <- clustering
  for (k_temp in 1:k){
    if (k_temp %in% aggregated_data$cluster){
      pcvariance_map[clustering[]==k_temp] <- aggregated_data$pc_variance[aggregated_data$cluster==k_temp][1]
    }else{
      pcvariance_map[clustering[]==k_temp] <- 0
    }
    
  }
  
  # Transform the raster data to a data frame
  r_df <- as.data.frame(pcvariance_map, xy = TRUE)
  combinations <- length(unique(r_df$global))-1
  
  blue_colors <- c("gray",colorRampPalette(c("lightblue", "darkblue"))(combinations))
  
  # Plot the map with average surface temperature and arrows (WITH WIDTH)
  #ggplot(r_df, aes(x = x, y = y)) +
  ggplot()+
    geom_raster(data = r_df, aes(x = x, y = y, fill = global))+
    scale_fill_gradientn(colours = blue_colors, na.value = NA)+
    geom_curve(data = arrow_data, aes(x = start_lon, y = start_lat, xend = end_lon, yend = end_lat, color=color, size=width),
               curvature = 0.2, arrow = arrow(type = "closed", length = unit(1, "mm")),
               linetype = "longdash") +
    scale_color_manual(values = c("black"),
                       labels = c(i=index_names[c])) +
    scale_size_continuous(range = c(min(arrow_data$width), max(arrow_data$width))) +
    #geom_text(data = text_data, aes(x = lon, y = lat, label = label), 
    #          size = 5, color = "black", fontface = "bold") +
    coord_fixed(ratio = 1.3) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude",
         fill = "Variance Explained", size = "Correlation", color="")+  
    theme(
      #legend.spacing = unit(0.1, "mm"),           # Adjust the space between legend items
      #legend.key.size = unit(0.6, "mm"),           # Adjust the size of the legend keys
      #legend.text = element_text(size = 5),        # Adjust the text size in the legend
      #legend.title = element_text(size = 10),       # Adjust the title size in the legend
      panel.grid.major = element_blank(),           # Remove major grid lines
      panel.grid.minor = element_blank()            # Remove minor grid lines
    )+
    coord_cartesian(ylim = c(-60, 84))
  
  p1
  
  #sabe plot  
  ggsave(filename = paste0("plot/fig1_",i,".png"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm")
  ggsave(filename = paste0("plot/fig1_",i,".pdf"), plot = p1, width = 20, height = 10, units = "cm")  
  ggsave(filename = paste0("plot/fig1_",i,".tiff"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm", compression = "lzw")
  
}

