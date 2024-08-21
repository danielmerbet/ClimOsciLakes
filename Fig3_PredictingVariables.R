library(raster)
library(ggplot2)
library(sf)
library(dplyr)

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
directory  <- paste0(directory_base, "data/",variable)
k <- 47
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


#Load clustering map
clustering <- raster(paste0(directory,"/clustering/global.tiff"))

#extract indexes per cluster
list_index <- list()
for (k_temp in 1:k){
  list_index[[k_temp]] <- (sort(unique(summary_ccf_pc10per$index[summary_ccf_pc10per$cluster==k_temp])))
}
sorted_strings <- sapply(list_index, function(x) paste(sort(x), collapse = "-"))
table(sorted_strings)
sorted_strings[sorted_strings==""] <- NA
sorted_strings <- toupper(sorted_strings)
sorted_strings[is.na(sorted_strings)] <- "Unrelated"

# Load XGBoots results
xgboost_data <- read.csv("XGBoost/save_data.csv")

list_var <- list()
for (k_temp in 1:k){
  list_var[[k_temp]] <- ((unique(xgboost_data$var[xgboost_data$cluster==k_temp])))
}
paste_strings <- sapply(list_var, function(x) paste((x), collapse = "-"))
table(paste_strings)
paste_strings[paste_strings==""] <- NA
#sorted_strings <- toupper(sorted_strings)
paste_strings[is.na(paste_strings)] <- "Unrelated"

#Set index map
vars_map_symbol <- clustering
for (k_temp in 1:k){
  vars_map_symbol[clustering[]==k_temp] <- paste_strings[k_temp]
}

# Transform the raster data to a data frame
r_df <- as.data.frame(vars_map_symbol, xy = TRUE)
combinations <- length(unique(r_df$global))-2

if (variable=="surftemp"){
  var_plot <- "a.                                                   Surface Temperature"
}else{
  var_plot <- "b.                                                    Bottom Temperature"
}

# Read in the shapefile
centroid_data <- st_read("clustering/centroids.shp") 

# Extract the coordinates and combine them with the DN field
coordinates <- st_coordinates(centroid_data)  # This extracts the coordinates
shape_df <- centroid_data %>%
  as.data.frame() %>%                       # Convert to a regular data frame
  mutate(x = coordinates[, 1],              # Add x coordinate
         y = coordinates[, 2])    

# Select the required columns including the DN field
shape_df <- shape_df %>%
  select(x, y, DN)  # Replace 'DN' with the actual name of the column if necessary

meteo_vars <- c("tas", "hurs", "pr", "ps", "rsds", "rlds", "sfcwind", "tasmax", "tasmin")

for (m in meteo_vars){
  shape_df[m] <- NA
}

c <- 0
for (k_temp in shape_df$DN){
  c <- c+1
  for (var in meteo_vars){
    data_temp_row <- xgboost_data[xgboost_data$cluster==k_temp & xgboost_data$var==var,]
    if (nrow(data_temp_row)>0){
      shape_df[var][c,] <- sum(data_temp_row$predict*100)
    }
  }
}

load("clustering/order_to_apply.RData")
load("clustering/climate_zones.RData")
climate_zones <- sapply(save_zones, function(x) paste(x, collapse = "-"))
climate_zones <- data.frame(climate_zone=climate_zones, cluster=1:47)
climate_zones <- climate_zones[order_to_apply,]
climate_zones$cluster_reorder <- 1:47

# Copy to reorder
clustering_reorder <- clustering
shape_df$DN_reorder <- NA
summary_ccf_pc10per$cluster_reorder <- NA

c<-0
for (l in order_to_apply){
  c <- c+1
  clustering_reorder[clustering[]==l] <- c
  summary_ccf_pc10per$cluster_reorder[summary_ccf_pc10per$cluster==l] <- c
  shape_df$DN_reorder[shape_df$DN==l] <- c
}

shapefile <- st_read("clustering/climate_zones_vectorized.shp")

shapefile$majority <- NA
#tropical climates
shapefile$majority[shapefile$X_majority %in% 1:3] <- "Tropical" 
#arid climates
shapefile$majority[shapefile$X_majority %in% 4:7] <- "Arid"
#temperate climates
shapefile$majority[shapefile$X_majority %in% 8:16] <- "Temperate"
#cold climates
shapefile$majority[shapefile$X_majority %in% 17:28] <- "Cold" 
#polar climates
shapefile$majority[shapefile$X_majority %in% 29:30] <- "Polar"

offset_point <- 0
# Create the global map
p <- ggplot() +
  geom_sf(data = shapefile, aes(fill = majority), color=NA) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  #theme(
  #  panel.grid.major = element_blank(),  # Remove major grid lines
  #  panel.grid.minor = element_blank()   # Remove minor grid lines
  #) +
  labs(title = var_plot,
       x = "Longitude",
       y = "Latitude") +
  geom_text(data = shape_df, aes(x = x, y = y, label = DN_reorder), 
            size = 3, vjust = -0.6, hjust = 1.0, color = "black") +
  # Add point for tas non-NA values
  geom_point(data = shape_df[!is.na(shape_df$tas) | !is.na(shape_df$tasmin) | !is.na(shape_df$tasmax), ], 
             aes(x = x, y = y), shape = 0, size = 2, color = "black") +
  # Add point for pr non-NA values
  geom_point(data = shape_df[!is.na(shape_df$pr), ], 
             aes(x = x- offset_point*2, y = y), shape = 1, size = 2, color = "black") +
  # Add point for hurs non-NA values
  geom_point(data = shape_df[!is.na(shape_df$sfcwind), ], 
             aes(x = x- offset_point, y = y), shape = 2, size = 2, color = "black")+
  # Add point for pr non-NA values
  geom_point(data = shape_df[!is.na(shape_df$rsds), ], 
             aes(x = x+offset_point, y = y), shape = 3, size = 2, color = "black") +
  # Add point for hurs non-NA values
  geom_point(data = shape_df[!is.na(shape_df$rlds), ], 
             aes(x = x+offset_point*2, y = y), shape = 4, size = 2, color = "black")

p
#save plot  
ggsave(filename = "plot/plot3/fig3.pdf", plot = p, width = 20, height = 15, units = "cm")  




