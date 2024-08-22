library(foreign)

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
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

# Load the NetCDF file
nc_file <- paste0(directory_base, "nc_laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc")
temperature_data <- brick(nc_file, varname = "surftemp")

# Calculate the average temperature over time
avg_temperature <- calc(temperature_data, mean, na.rm = TRUE)

# Convert the averaged temperature data to a data frame
temperature_df <- as.data.frame(rasterToPoints(avg_temperature))
colnames(temperature_df) <- c("lon", "lat", "temperature")

summary_ccf_pc10per$lon_index <- NA
summary_ccf_pc10per$lat_index <- NA

#NAO location
summary_ccf_pc10per[summary_ccf_pc10per$index=="nao",]$lon_index <- -25
summary_ccf_pc10per[summary_ccf_pc10per$index=="nao",]$lat_index <- 52

#ENSO location for  clusters in America
summary_ccf_pc10per$lon_index[(summary_ccf_pc10per$cluster %in% c(1,2,4:22)) & (summary_ccf_pc10per$index=="oni")] <- -145
summary_ccf_pc10per$lat_index[(summary_ccf_pc10per$cluster %in% c(1,2,4:22)) & (summary_ccf_pc10per$index=="oni")] <- 0
#ENSO location for  clusters in Africa and Asia
summary_ccf_pc10per$lon_index[(summary_ccf_pc10per$cluster %in% c(3,23:47)) & (summary_ccf_pc10per$index=="oni")] <- 160
summary_ccf_pc10per$lat_index[(summary_ccf_pc10per$cluster %in% c(3,23:47)) & (summary_ccf_pc10per$index=="oni")] <- 0

#PDO location
summary_ccf_pc10per[summary_ccf_pc10per$index=="pdo",]$lon_index <- -140
summary_ccf_pc10per[summary_ccf_pc10per$index=="pdo",]$lat_index <- 45

#IOD location
summary_ccf_pc10per[summary_ccf_pc10per$index=="iod",]$lon_index <- 80
summary_ccf_pc10per[summary_ccf_pc10per$index=="iod",]$lat_index <- -5

summary_ccf_pc10per$color_arrow <- NA

summary_ccf_pc10per$color_arrow[summary_ccf_pc10per$index=="oni"] <- "darkred"
summary_ccf_pc10per$color_arrow[summary_ccf_pc10per$index=="nao"] <- "black"
summary_ccf_pc10per$color_arrow[summary_ccf_pc10per$index=="pdo"] <- "blue4"
summary_ccf_pc10per$color_arrow[summary_ccf_pc10per$index=="iod"] <- "purple4"

arrow_data <- data.frame(
  start_lon = summary_ccf_pc10per$lon_index, # Longitudes of start points
  start_lat = summary_ccf_pc10per$lat_index,  # Latitudes of start points
  end_lon = summary_ccf_pc10per$lon,  # Longitudes of end points
  end_lat = summary_ccf_pc10per$lat,   # Latitudes of end points
  color = summary_ccf_pc10per$index, #Color of the arrow for each index  
  width = (summary_ccf_pc10per$cor) #relate the size of the arrow with the correlation with the indexes
)

# Define index as text for the map
text_data <- data.frame(
  lon = c(-25, -145, 160, -140, 80),  # Longitudes where you want the text
  lat = c(52, 0, 0, 45, -5),      # Latitudes where you want the text
  label = c("NAO", "ENSO", "ENSO", "PDO", "IOD")  # Text labels
)

# Plot the map with average surface temperature and arrows (WITH WIDTH)
ggplot() +
  geom_raster(data = temperature_df, aes(x = lon, y = lat, fill = temperature-273.15)) +  # Plot the temperature data
  #scale_fill_gradientn(colours = rev(terrain.colors(10))) +
  #scale_fill_viridis_c(option = "cividis")+
  #scale_fill_gradientn(colours = heat.colors(10))+
  #scale_fill_gradientn(colours = topo.colors(10))+
  #scale_fill_gradientn(colours = rainbow(10))+
  scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")))+
  geom_curve(data = arrow_data, aes(x = start_lon, y = start_lat, xend = end_lon, yend = end_lat, color=color, size=width),
             curvature = 0.2, arrow = arrow(type = "closed", length = unit(1, "mm"))) +
  scale_color_manual(values = c("oni" = "darkred", 
                                "nao" = "black", 
                                "pdo" = "blue4",
                                "iod" = "purple4"),
                     labels = c("oni"="ENSO","nao"="NAO","pdo"="PDO","iod"="IOD")) +
  scale_size_continuous(range = c(min(arrow_data$width), max(arrow_data$width))) +
  geom_text(data = text_data, aes(x = lon, y = lat, label = label), 
            size = 5, color = "black", fontface = "bold") +
  coord_fixed(ratio = 1.3) +
  theme_minimal() +
  labs(fill = "Degree C", color = "Climate index", size = "Correlation")

# Plot the map with average surface temperature and arrows (NO WIDTH)
ggplot() +
  geom_raster(data = temperature_df, aes(x = lon, y = lat, fill = temperature-273.15)) +  # Plot the temperature data
  #scale_fill_gradientn(colours = rev(terrain.colors(10))) +
  #scale_fill_viridis_c(option = "cividis")+
  #scale_fill_gradientn(colours = heat.colors(10))+
  #scale_fill_gradientn(colours = topo.colors(10))+
  #scale_fill_gradientn(colours = rainbow(10))+
  scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")))+
  geom_curve(data = arrow_data, aes(x = start_lon, y = start_lat, xend = end_lon, yend = end_lat, color=color, size=width),
             curvature = 0.2, arrow = arrow(type = "closed", length = unit(1, "mm")),
            size = 0.1) +
  scale_color_manual(values = c("oni" = "darkred", 
                                "nao" = "black", 
                                "pdo" = "darkslategray",
                                "iod" = "purple4"),
                     labels = c("oni"="ENSO","nao"="NAO","pdo"="PDO","iod"="IOD")) +
  coord_fixed(ratio = 1.3) +
  theme_minimal() +
  labs(fill = "Degree C", color = "Climate index")

#Plot a gray map
ggplot() +
  borders("world", colour = "gray85", fill = "gray80") +  # Plot the world map
  theme_minimal() +
  geom_curve(data = arrow_data, aes(x = start_lon, y = start_lat, xend = end_lon, yend = end_lat),
             curvature = 0.2, arrow = arrow(type = "closed", length = unit(2, "mm")),
             color = "blue", size = 0.1) +
  coord_fixed(ratio = 1.3)
