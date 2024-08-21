library(colorRamps)
library(raster)
library(dplyr)
library(ggplot2)
library(paletteer)
library(foreign)

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
    
  }
}
#################33
# Order by climate zones according to Beck et al. (2018):
climate_zones <- read.dbf("clustering/climate_zones.dbf")
num_tropical <- sort(climate_zones$DN[climate_zones$X_majority %in% 1:3])
num_dry <- sort(climate_zones$DN[climate_zones$X_majority %in% 4:7])
num_temperate <- sort(climate_zones$DN[climate_zones$X_majority %in% 8:16])
num_continental <- sort(climate_zones$DN[climate_zones$X_majority %in% 17:28])
num_polar <- sort(climate_zones$DN[climate_zones$X_majority %in% 29:30])
order_cz <- c(num_tropical, num_dry, num_temperate, num_continental, num_polar)

for (i in c("enso", "nao", "pdo", "iod")){
  
  df <- all_data[all_data$index==i,]
  
  uni_cluster <- unique(df$cluster)
  l_tropical <- length(which(uni_cluster %in% num_tropical))
  l_dry <- length(which(uni_cluster %in% num_dry))
  l_temperate <- length(which(uni_cluster %in% num_temperate))
  l_continental <- length(which(uni_cluster %in% num_continental))
  l_polar <- length(which(uni_cluster %in% num_polar))
  
  # Specify the order of clusters
  df$cluster <- factor(df$cluster, levels = order_cz)
  
  p2 <- ggplot(df, aes(x = cluster, y = slope_mean, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("positive" = "black", "negative" = "gray")) +
    labs(x = "Cluster", fill="Index phase", y = "Mean Decrease/Increase (%)", title = toupper(i)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
    theme(legend.position = "none")+
    annotate("rect", xmin = 0.5, xmax = l_tropical+0.5, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = l_tropical+0.5,
             xmax = sum(l_tropical+l_dry)+0.5, 
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "goldenrod") +
    annotate("rect", xmin = sum(l_tropical+l_dry)+0.5, 
             xmax = sum(l_tropical+l_dry+l_temperate)+0.5,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "green") +
    annotate("rect", xmin = sum(l_tropical+l_dry+l_temperate)+0.5, 
             xmax = sum(l_tropical+l_dry+l_temperate+l_continental)+0.5,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = sum(l_tropical+l_dry+l_temperate+l_continental)+0.5, 
             xmax = sum(l_tropical+l_dry+l_temperate+l_continental+l_polar)+0.5,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "violet") +
    annotate("text", x = l_tropical/2+0.5, y = max(df$slope_mean) + 5, label = "Tropical", size = 4, color = "black") +
    annotate("text", x = sum(l_tropical+l_dry/2)+0.5, y = max(df$slope_mean) + 5, label = "Dry", size = 4, color = "black")+ 
    annotate("text", x = sum(l_tropical+l_dry+l_temperate/2)+0.5, y = max(df$slope_mean) + 5, label = "Temperate", size = 4, color = "black")+
    annotate("text", x = sum(l_tropical+l_dry+l_temperate+l_continental/2)+0.5, y = max(df$slope_mean) + 5, label = "Continental", size = 4, color = "black") +
    annotate("text", x = sum(l_tropical+l_dry+l_temperate+l_continental+l_polar/2)+0.5, y = max(df$slope_mean) + 5, label = "Polar", size = 4, color = "black")
  
  # Plot the data
  #p2 <- ggplot(df, aes(x = factor(cluster), y = slope_mean, fill = direction)) +
  #  geom_bar(stat = "identity", position = "dodge") +
  #  scale_fill_manual(values = c("positive" = "darkred", "negative" = "blue4")) +
  #  labs(x = "Cluster", y = "Slope Mean", title = toupper(i)) +
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #  theme_minimal()
  
  ggsave(filename = paste0("plot/plot2/fig2_",i,".pdf"), plot = p2, width = 20, height = 10, units = "cm")  
  
}
