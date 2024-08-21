library(colorRamps)
library(raster)
library(dplyr)
library(ggplot2)
library(paletteer)

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

# Read oscillations table
oscillations <- read.table(paste0(directory_base, "clim_index/oscillations.txt"), header = T)
colnames(oscillations) <- c("year", "enso", "nao", "pdo", "iod", "npgo", "npo")

###################################################
#ANALYSIS USING PRINCIPAL COMPONENTS

#data_temp <- summary_ccf_pc10per
#row<-1
#direction_index <- "positive"
#i <- "enso"

all_data <- c()
for (i in c("enso", "nao", "pdo", "iod")){
  for (direction_index in c("positive", "negative")){
    
    slope_save <- c(); intercept_save <- c()
    for (row in 1:nrow(summary_ccf_pc10per)){
      cluster_temp <- summary_ccf_pc10per$cluster[row]
      index_temp <- summary_ccf_pc10per$index[row]
      pc_temp <- paste0("PC", summary_ccf_pc10per$pc[row])
      pca_data <- read.csv(paste0("pc_data/global", cluster_temp,".csv"))
      
      pca_data_index <- data.frame(oscillations[index_temp], 
                                   pca_data[pc_temp])
      colnames(pca_data_index) <- c("index", "pc")
      
      if (direction_index=="positive"){
        pca_data_index <- pca_data_index[pca_data_index$index>0,]
      }else{
        pca_data_index <- pca_data_index[pca_data_index$index<0,]
      }
      
      
      #plot
      if (i==index_temp){
        pdf(file=paste0("plot/plot2/regression/",i,"_",
                        direction_index,"_c",cluster_temp,
                        "_i",index_temp,"_pc",pc_temp,".pdf"))
        plot(pc ~ index, data=pca_data_index)
        abline(lm(pc ~ index, data=pca_data_index), col = "blue", lwd = 2)
        dev.off()
      }
      
      # Fit the linear regression model
      model <- lm(pc ~ index, data=pca_data_index)
      #add regression to plot
      #abline(model, col="red")
      # Extract the slope and intercept
      coefficients <- coef(model)
      slope_save <- c(slope_save, coefficients[2])
      intercept_save <- c(intercept_save, coefficients[1])
    }
    summary_ccf_pc10per$slope <- slope_save
    summary_ccf_pc10per$intercept <- intercept_save
    
    #Load clustering map
    clustering <- raster(paste0(directory,"/clustering/global.tiff"))
    
    data_temp <- summary_ccf_pc10per[summary_ccf_pc10per$index==i,]
    
    # Aggregate by cluster using sum function
    aggregated_data <- data_temp %>%
      group_by(cluster) %>%
      summarize(
        pc_variance_sum = sum(pc_variance),
        slope_mean = mean(slope),
        intercept_mean = mean(intercept)
      )
    
    aggregated_data$direction <- direction_index
    aggregated_data$index <- i
    
    all_data <- rbind(all_data, aggregated_data)
    
    #print(range(aggregated_data$slope_mean))
    #print((aggregated_data$slope_mean))
    
    
    #Set map
    slope_map_symbol <- clustering
    for (k_temp in 1:k){
      if (k_temp %in% aggregated_data$cluster){
        slope_map_symbol[clustering[]==k_temp] <- aggregated_data$slope_mean[aggregated_data$cluster==k_temp][1]
        #pcvariance_map[clustering[]==k_temp] <- aggregated_data$pc_variance[aggregated_data$cluster==k_temp][1]
      }else{
        #slope_map_symbol[clustering[]==k_temp] <- round(sort(unique(aggregated_data$slope_mean))[1]-1)
        slope_map_symbol[clustering[]==k_temp] <- NA
      }
      
    }
    
    #limits = c(sort(unique(r_df$global))[2],range(r_df$global, na.rm=T)[2]))
    # Transform the raster data to a data frame
    r_df <- as.data.frame(slope_map_symbol, xy = TRUE)
    combinations <- length(unique(r_df$global))-2
    
    #plot_colors <- c("gray",colorRampPalette(c("darkblue","yellow", "darkred"))(combinations))
    #plot_colors <- c("gray", blue2red(10))
    
    
    p1 <- ggplot() +
      geom_raster(data = r_df, aes(x = x, y = y, fill = global)) +
      #scale_fill_viridis_c(option = "cividis", na.value = NA)+
      #scale_fill_gradientn(colours = c("pink",viridisLite::cividis(combinations)), na.value = NA)+
      scale_fill_gradientn(colours = c(rev(paletteer_c("grDevices::Spectral", combinations))), na.value = NA)+ #, limits=c(-24.25,42.23)
      #scale_fill_gradientn(colours = c(rev(paletteer_c("ggthemes::Classic Red-Blue", 30))), na.value = NA)+
      
      #scale_fill_gradientn(colours = plot_colors, na.value = NA)+
      #scale_size_continuous(range = c(sort(unique(r_df$global))[2], range(r_df$global, na.rm=T)[2])) +
      #scale_fill_gradientn(colours = rev(terrain.colors(combinations)), na.value = NA )+
      theme_minimal() +
      coord_fixed() +
      labs(x = 'Longitude', y = 'Latitude', fill = 'Slope') +
      #ggtitle('Global Lake Response to Climate Oscillations')+
      guides(
        #shape = guide_legend(override.aes = list(size = 4)),
        color = guide_legend(override.aes = list(size = 4))
      )+
      theme(
        #legend.spacing = unit(0.1, "mm"),           # Adjust the space between legend items
        #legend.key.size = unit(0.6, "mm"),           # Adjust the size of the legend keys
        #legend.text = element_text(size = 8),        # Adjust the text size in the legend
        #legend.title = element_text(size = 10),       # Adjust the title size in the legend
        panel.grid.major = element_blank(),           # Remove major grid lines
        panel.grid.minor = element_blank()            # Remove minor grid lines
      )+
      coord_cartesian(ylim = c(-60, 84))
    #p1
    
    
    #ggsave(filename = paste0("plot/fig2_",i,"_",direction_index, ".png"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm")
    #ggsave(filename = paste0("plot/plot2/alternative/fig2_",i,"_",direction_index,".pdf"), plot = p1, width = 20, height = 10, units = "cm")  
    #ggsave(filename = paste0("plot/fig2_",i,"_",direction_index,".tiff"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm", compression = "lzw")
    
    
    
    
    
  }
}
#################33

for (i in c("enso", "nao", "pdo", "iod")){
  
  df <- all_data[all_data$index==i,]
  
  # Plot the data
  p2 <- ggplot(df, aes(x = factor(cluster), y = slope_mean, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("positive" = "darkred", "negative" = "blue4")) +
    labs(x = "Cluster", y = "Slope Mean", title = toupper(i)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_minimal()
  
  ggsave(filename = paste0("plot/plot2/fig2_",i,".pdf"), plot = p2, width = 20, height = 10, units = "cm")  
  
}
  



#clustering with quanlitative (direction) and quantitave (slope) variables
# Normalize the quantitative data
df_normalized <- df %>%
  mutate(across(c(pc_variance_sum, slope_mean, intercept_mean), scale))

# Compute the distance matrix
dist_matrix <- dist(df_normalized)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Convert hclust object to dendrogram object
dend <- as.dendrogram(hc)

# Create a ggplot dendrogram
dend_data <- dendro_data(dend)

# Plot the dendrogram
ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(dend_data), aes(x = x, y = y, label = label, hjust = 1), size = 3) +
  theme_minimal() +
  labs(title = "Dendrogram of Clusters (Including Direction)", x = "Clusters", y = "Height") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

df_normalized$direction<- as.factor(df_normalized$direction)
df_normalized$slope_mean <- as.numeric(df_normalized$slope_mean)

# Compute Gower distance
gower_dist <- daisy(df_normalized %>% select(-index, -pc_variance_sum, -intercept_mean), metric = "gower")

# Perform hierarchical clustering
hc <- hclust(gower_dist, method = "ward.D2")

# Convert hclust object to dendrogram object
dend <- as.dendrogram(hc)

# Create a ggplot dendrogram
dend_data <- dendro_data(dend)

# Plot the dendrogram
ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(dend_data), aes(x = x, y = y, label = label, hjust = 1), size = 3) +
  theme_minimal() +
  labs(title = "Dendrogram of Clusters (Including Direction)", x = "Clusters", y = "Height") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))







#plot NA
#Set map
slope_map_symbol <- clustering
for (k_temp in 1:k){
  if (k_temp %in% aggregated_data$cluster){
    slope_map_symbol[clustering[]==k_temp] <- NA
    #pcvariance_map[clustering[]==k_temp] <- aggregated_data$pc_variance[aggregated_data$cluster==k_temp][1]
  }else{
    #slope_map_symbol[clustering[]==k_temp] <- round(sort(unique(aggregated_data$slope_mean))[1]-1)
    slope_map_symbol[clustering[]==k_temp] <- 20
  }
  
}

#limits = c(sort(unique(r_df$global))[2],range(r_df$global, na.rm=T)[2]))
# Transform the raster data to a data frame
r_df <- as.data.frame(slope_map_symbol, xy = TRUE)
combinations <- length(unique(r_df$global))-2

#plot_colors <- c("gray",colorRampPalette(c("darkblue","yellow", "darkred"))(combinations))
#plot_colors <- c("gray", blue2red(10))


p1 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = global)) +
  #scale_fill_viridis_c(option = "cividis", na.value = NA)+
  #scale_fill_gradientn(colours = c("pink",viridisLite::cividis(combinations)), na.value = NA)+
  scale_fill_gradientn(colours = c("gray"), na.value = NA)+
  
  #scale_fill_gradientn(colours = plot_colors, na.value = NA)+
  #scale_size_continuous(range = c(sort(unique(r_df$global))[2], range(r_df$global, na.rm=T)[2])) +
  #scale_fill_gradientn(colours = rev(terrain.colors(combinations)), na.value = NA )+
  theme_minimal() +
  coord_fixed() +
  labs(x = 'Longitude', y = 'Latitude', fill = 'Slope') +
  #ggtitle('Global Lake Response to Climate Oscillations')+
  guides(
    #shape = guide_legend(override.aes = list(size = 4)),
    color = guide_legend(override.aes = list(size = 4))
  )+
  theme(
    #legend.spacing = unit(0.1, "mm"),           # Adjust the space between legend items
    #legend.key.size = unit(0.6, "mm"),           # Adjust the size of the legend keys
    #legend.text = element_text(size = 8),        # Adjust the text size in the legend
    #legend.title = element_text(size = 10),       # Adjust the title size in the legend
    panel.grid.major = element_blank(),           # Remove major grid lines
    panel.grid.minor = element_blank()            # Remove minor grid lines
  )+
  coord_cartesian(ylim = c(-60, 84))
p1

ggsave(filename = paste0("plot/fig2_",i,"_",direction_index,"_NA.png"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm")
ggsave(filename = paste0("plot/fig2_",i,"_",direction_index,"_NA.pdf"), plot = p1, width = 20, height = 10, units = "cm")  
ggsave(filename = paste0("plot/fig2_",i,"_",direction_index,"_NA.tiff"), plot = p1, width = 20, height = 10, dpi = 300, units = "cm", compression = "lzw")














###################################################
#ANALYSIS USING WATER TEMP DATA
#Load temperature data
load("HiClimR_output/global.RData")

#ENSO
pos_max <- which(oscillations$oni==max(oscillations$oni))
oscillations$year[pos_max]
pos_min <- which(oscillations$oni==min(oscillations$oni))
oscillations$year[pos_min]

#MAXIMUN VALUE
#Convert to raster data
temp_data<- y$data[,pos_max]
lon_lat_names <- strsplit(names(y$data[,1]), ",")
# Convert the list to a data frame
df_lon_lat <- do.call(rbind, lon_lat_names)
df_lon_lat <- data.frame(df_lon_lat, stringsAsFactors = FALSE)
names(df_lon_lat) <- c("lon", "lat")

temp_data <- data.frame(lon=as.numeric(df_lon_lat$lon),
                        lat=as.numeric(df_lon_lat$lat),
                        temp=as.numeric(y$data[,1]),
                        cluster=y$region[!is.na(y$region)])
temp_data_r <- temp_data
# Create a SpatialPointsDataFrame
coordinates(temp_data_r) <- ~lon+lat
# Define the resolution of the raster (based on the data provided)
res <- 0.5
# Create an empty raster
temp_r <- raster(extent(temp_data_r), resolution=res)
# Assign values to the raster
temp_r <- rasterize(temp_data_r, r, field="temp")



r_df <- as.data.frame(clustering, xy = TRUE)
combinations <- length(unique(r_df$global))-1

blue_colors <- c("gray",colorRampPalette(c("lightblue", "darkblue"))(combinations))


p1 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = global))
  
ggplot()


     