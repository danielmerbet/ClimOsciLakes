library(raster)
library(ggplot2)
library(sf)
library(RColorBrewer)

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
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


# Base world map
world_map <- map_data("world")

#my_colors <- c("#5E4FA2", "#7570B3","#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF",
#               "#FEE08B", "#FDAE61", "#F46D43", "salmon", "#D53E4F","#9E0142")

# Function to create a bar plot for a single location
create_bar_grob <- function(tas,hurs,pr,ps,rsds,rlds,sfcwind,tasmax,tasmin) {
  # Create a data frame for the bar plot, removing zero values
  plot_data <- data.frame(value = c(tas,hurs,pr,ps,rsds,rlds,sfcwind,tasmax,tasmin),
                          category = c("tm","rh","pr","ps","sr",
                                       "lr","ws","tx","tn"))
  plot_data <- plot_data[plot_data$value != 0, ]
  
  # Create the bar plot
  bar_plot <- ggplot(plot_data) +
    geom_bar(aes(x = category, y = value, fill = category), stat = "identity", width = 0.5) +
    #geom_text(aes(x = category, y = value, label = category), vjust = -0.5, size=1.5) +
    scale_fill_manual(values = c("tm" = "salmon", "rh" = "#F46D43", "pr" = "#3288BD",
                                 "ps"= "#66C2A5", "sr" ="#E6F598", "lr"= "#FFFFBF",
                                 "ws" = "#FEE08B", "tx" ="#9E0142" , "tn"="#D53E4F"
    )) +
    theme_void() +
    theme(legend.position = "none")
  
  ggplotGrob(bar_plot)
}

shapefile <- st_read("clustering/climate_zones_vectorizes.shp")

# Create the global map
p <- ggplot() +
  geom_sf(data = shapefile, aes(fill = X_majority)) + #
  coord_sf(expand = FALSE) +
  #geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray90", color = "gray60") +
  #coord_fixed(ratio = 1.3) +
  theme_minimal() +
  labs(title = var_plot,
       x = "Longitude",
       y = "Latitude")+
  geom_text(data = shape_df, aes(x = x, y = y, label = DN_reorder), 
            size = 3, vjust = -0.6, hjust= 0.8,color = "black")

# Add bar plots to the map
for (i in 1:nrow(shape_df)) {
  bar_grob <- create_bar_grob(shape_df$tas[i], shape_df$hurs[i], shape_df$pr[i],
                              shape_df$ps[i], shape_df$rsds[i],shape_df$rlds[i],
                              shape_df$sfcwind[i], shape_df$tasmax[i], shape_df$tasmin[i])
  x <- shape_df$x[i]
  y <- shape_df$y[i]
  p <- p + annotation_custom(bar_grob, xmin = x, xmax = x + 10, ymin = y, ymax = y + 10)
}

# Print the plot
#print(p)


# Create a legend
legend_data <- data.frame(
  category = c("tm","rh","pr","ps","sr",
               "lr","ws","tx","tn"),
  value = c(1, 1, 1,1,1,1,1,1,1)
)

legend_plot <- ggplot(legend_data) +
  geom_bar(aes(x = category, y = value, fill = category), stat = "identity") +
  scale_fill_manual(values = c("tm" = "salmon", "rh" = "#F46D43", "pr" = "#3288BD",
                               "ps"= "#66C2A5", "sr" ="#E6F598", "lr"= "#FFFFBF",
                               "ws" = "#FEE08B", "tx" ="#9E0142" , "tn"="#D53E4F")) +
  theme_void() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "", nrow = 1))

# Extract legend as a grob
legend_grob <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]

# Add legend to the map
p1 <- p + annotation_custom(legend_grob, xmin = -180, xmax = 180, ymin = -90, ymax = -80)

#save plot  
#ggsave(filename = "plot/fig1.png", plot = p1, width = 20, height = 10, dpi = 300, units = "cm")
ggsave(filename = "plot/plot3/fig3.pdf", plot = p1, width = 20, height = 15, units = "cm")  
#ggsave(filename = "plot/fig1.tiff", plot = p1, width = 20, height = 10, dpi = 300, units = "cm", compression = "lzw")




