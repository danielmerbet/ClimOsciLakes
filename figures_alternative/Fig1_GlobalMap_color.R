library(raster)
library(ggplot2)
library(sf)
library(dplyr)

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
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


#Load clustering map
clustering <- raster(paste0(directory,"/clustering/global.tiff"))

# Load centroids data
centroid_data <- st_read("clustering/centroids.shp") 

# Extract the coordinates and combine them with the DN field
coordinates <- st_coordinates(centroid_data)  # This extracts the coordinates
shape_df <- centroid_data %>%
  as.data.frame() %>%                       # Convert to a regular data frame
  mutate(x = coordinates[, 1],              # Add x coordinate
         y = coordinates[, 2])    

#load climate zone data and reorder accordingly
load("clustering/climate_zones.RData")
climate_zones <- sapply(save_zones, function(x) paste(x, collapse = "-"))
custom_order <- c("tropical", "tropical-dry","tropical-dry-temperate",
                  "tropical-dry-temperate-continental",
                  "tropical-dry-temperate-polar","tropical-temperate",
                  "dry", "dry-temperate", "dry-temperate-continental",
                  "dry-temperate-continental-polar","dry-temperate-polar",
                  "dry-continental", "dry-polar",
                  "temperate","temperate-continental",
                  "continental","continental-polar", "polar")
get_rank <- function(x) {
  # Find the index of the first match in custom_order
  match_index <- match(x, custom_order)
  # Handle cases where the exact match isn't found (e.g., if there's a typo)
  if (is.na(match_index)) {
    # Loop through the custom order to find the first pattern that matches
    for (i in seq_along(custom_order)) {
      if (grepl(custom_order[i], x)) {
        return(i)
      }
    }
  }
  return(match_index)
}

# Sort the vector using the custom order
order_to_apply <- order(sapply(climate_zones, get_rank))
#sorted_vector <- climate_zones[order_to_apply]
save(order_to_apply, file = "clustering/order_to_apply.RData")

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

#extract indexes per cluster
list_index <- list()
for (k_temp in 1:k){
  list_index[[k_temp]] <- (sort(unique(summary_ccf_pc10per$index[summary_ccf_pc10per$cluster_reorder==k_temp])))
}
sorted_strings <- sapply(list_index, function(x) paste(sort(x), collapse = "-"))
table(sorted_strings)
sorted_strings[sorted_strings==""] <- NA
sorted_strings <- toupper(sorted_strings)
sorted_strings[is.na(sorted_strings)] <- "Unrelated"

#Set index map
indexes_map_symbol <- clustering_reorder
for (k_temp in 1:k){
  indexes_map_symbol[clustering_reorder[]==k_temp] <- sorted_strings[k_temp]
}

# Transform the raster data to a data frame
r_df <- as.data.frame(indexes_map_symbol, xy = TRUE)
combinations <- length(unique(r_df$global))-2

# Set title for plot
if (variable=="surftemp"){
  var_plot <- "a.                                                   Surface Temperature"
}else{
  var_plot <- "b.                                                    Bottom Temperature"
}

factor_plot <- c("ENSO-NAO-PDO", "ENSO-IOD-NAO","ENSO-IOD-NAO-PDO",
                 "NAO", "ENSO-PDO", "ENSO-NAO", "NAO-PDO",     
                 "ENSO-IOD-PDO", "IOD-PDO", "PDO",            
                 "IOD-NAO", "ENSO-IOD", "IOD-NAO-PDO", "ENSO", "Unrelated")

my_colors <- c(colorRampPalette(c("#053061","lightblue",
                                  "lightcoral","#67001F"
                                  ))(length(factor_plot)-1), "gray")

sel_colors <- which(sort(factor_plot) %in% unique(r_df$global))
my_colors <- my_colors[sel_colors]

#Final plot
p1 <- ggplot(r_df, aes(x = x, y = y)) +
  #geom_point(aes(shape = factor(global), color = factor(global)), size = 0.1) +
  geom_point(aes(color = factor(global)), size = 0.1, alpha=0.6) +
  #scale_color_brewer(palette = "RdBu")+
  scale_color_manual(values=my_colors, na.translate = FALSE)+
  #for legend:
  #scale_color_manual(values=rep("white",length(my_colors)), na.translate = FALSE)+
  theme_minimal() +
  #coord_sf(expand = FALSE) +
  coord_fixed() +
  labs(x = 'Longitude', y = 'Latitude', 
       title = var_plot, color = 'Significant relation') +
  #ggtitle('Global Lake Response to Climate Oscillations')+
  guides(
    #shape = guide_legend(override.aes = list(size = 4)),
    color = guide_legend(override.aes = list(size = 4))
  )+
  theme(
    legend.spacing = unit(0.1, "mm"),           # Adjust the space between legend items
    legend.key.size = unit(0.6, "mm"),           # Adjust the size of the legend keys
    legend.text = element_text(size = 8),        # Adjust the text size in the legend
    legend.title = element_text(size = 10),       # Adjust the title size in the legend
    #panel.grid.major = element_blank(),           # Remove major grid lines
    panel.grid.minor = element_blank()            # Remove minor grid lines
  )+
  #coord_cartesian(ylim = c(-60, 84))+
  #scale_y_continuous(breaks = seq(-40, 80, by = 20), limits = c(-60, 84)) +
  #scale_x_continuous(breaks = seq(-120, 120, by = 60)) +
  scale_y_continuous(
    breaks = seq(-40, 80, by = 20),
    limits = c(-60, 84),
    labels = function(lat) {
      sapply(lat, function(y) {
        if (y > 0) {
          paste0(y, "\u00B0N")
        } else if (y < 0) {
          paste0(-y, "\u00B0S")
        } else {
          "0\u00B0"
        }
      })
    }
  ) +
  scale_x_continuous(
    breaks = seq(-120, 120, by = 60),
    labels = function(lon) {
      sapply(lon, function(x) {
        if (x > 0) {
          paste0(x, "\u00B0E")
        } else if (x < 0) {
          paste0(-x, "\u00B0W")
        } else {
          "0\u00B0"
        }
      })
    }
  ) +
  geom_text(data = shape_df, aes(x = x, y = y, label = DN_reorder), 
            size = 3, vjust = 0, color = "black") 
p1

#save plot  
ggsave(filename = "plot/plot1/fig1.pdf", plot = p1, width = 20, height = 15, units = "cm", dpi=300)  
