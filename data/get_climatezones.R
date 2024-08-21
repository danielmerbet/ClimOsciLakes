# Load required libraries
library(raster)
library(sf)    # Alternatively to rgdal and rgeos
library(dplyr) # For data manipulation

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
k <- 47 #number of total clusters
directory  <- paste0(directory_base, "data/",variable)
setwd(directory)

# Load the raster
raster_data <- raster("../../clim_data/koppen_geiger_0p00833333_1991-2020.tif")

# Load the shapefile
shapefile_data <- st_read("clustering/clustering.shp")  # Using sf package

# Climate zones according to Beck et al. (2018):
table_beck <- data.frame(id=1:30, zone=c(rep("tropical",3),
                                         rep("dry",4),
                                         rep("temperate",9),
                                         rep("continental",12),
                                         rep("polar",2)))
# Create list to save
save_list <- list(); save_zones <- list()

for (dn in 1:k){
  
  # Filter the shapefile based on the cluster attribute
  filtered_shape <- shapefile_data %>% filter(DN == dn)
  
  # Extract values where the shapefile intersects the raster
  extracted_values <- extract(raster_data, filtered_shape)
  extracted_values <- unlist(extracted_values)
  
  # Display and save the extracted values
  table_values <- table(extracted_values)
  print(dn)
  print(table_values)
  
  # Select for more than 1% of pixels
  percent_values <- (table_values/sum(table_values))*100
  percent_values <- percent_values[percent_values>5]
  
  save_list[[dn]] <- table_values
  zones_values <- unique(table_beck$zone[table_beck$id %in% as.numeric(names(percent_values))])
  print(zones_values)
  save_zones[[dn]] <- zones_values
}

save(save_zones, file = "clustering/climate_zones.RData")





 