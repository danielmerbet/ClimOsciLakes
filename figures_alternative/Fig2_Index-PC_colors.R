library(colorRamps)
library(raster)
library(dplyr)
library(ggplot2)
library(paletteer)
library(foreign)
library(reshape2)

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

#Add order an climate zone data
load("clustering/order_to_apply.RData")
load("clustering/climate_zones.RData")
climate_zones <- sapply(save_zones, function(x) paste(x, collapse = "-"))
climate_zones <- data.frame(climate_zone=climate_zones, cluster=1:47)
climate_zones <- climate_zones[order_to_apply,]
climate_zones$cluster_reorder <- 1:47

# Merge regression data with climate data and add new order
all_data <- merge(all_data, climate_zones, by=c("cluster"))

for (i in c("enso", "nao", "pdo", "iod")){
  
  df <- all_data[all_data$index==i,]
  
  df_long <- df%>%
    group_by(cluster_reorder) %>%
    slice_max(abs(slope_mean), with_ties = FALSE)
  
  # Create a new dataframe with separate rows for each match
  df_long <- df_long %>%
    mutate(tropical = grepl("tropical", climate_zone),
           dry = grepl("dry", climate_zone),
           temperate = grepl("temperate", climate_zone),
           continental = grepl("continental", climate_zone),
           polar = grepl("polar", climate_zone)) %>%
    melt(id.vars = c("cluster","cluster_reorder", "slope_mean", "direction", "climate_zone"),
         measure.vars = c("tropical", "dry", "temperate", "continental", "polar"),
         variable.name = "keyword",
         value.name = "present") %>%
    filter(present == TRUE)
  
  # Map colors to keywords
  df_long$circle_color <- with(df_long, 
                               ifelse(keyword == "tropical", "blue",
                                      ifelse(keyword == "dry", "red",
                                             ifelse(keyword == "temperate", "yellow",
                                                    ifelse(keyword == "continental", "green",
                                                           ifelse(keyword == "polar", "violet", NA))))))
  
  if (variable=="bottemp"){
    if (i=="enso"){
      offset_value <- 6
    }else if(i=="iod"){
      offset_value <- 10
    }else if(i=="pdo"){
      offset_value <- 4
    }else{
      offset_value <- 5
    }
  }else{
    if (i=="enso"){
      offset_value <- 5
    }else if(i=="iod"){
      offset_value <- 10
    }else if(i=="pdo"){
      offset_value <- 4
    }else{
      offset_value <- 8
    }
  }

  
  # Step 3: Create a stacking offset for each circle
  df_long <- df_long %>%
    group_by(cluster_reorder) %>%
    mutate(offset = row_number() * offset_value) 
  
  df_long[df_long$slope_mean<0,]$offset <-  (-1)*df_long[df_long$slope_mean<0,]$offset
  
  # Specify the order of clusters
  df$cluster_reorder <- factor(df$cluster_reorder, levels = sort(unique(df$cluster_reorder)))
  df_long$cluster_reorder <- factor(df_long$cluster_reorder, levels = sort(unique(df_long$cluster_reorder)))
  

  
  #df_long <- df_long[df_long$direction=="negative",]
  

  
  p2 <- ggplot(df, aes(x = cluster_reorder, y = slope_mean, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "green") +  # Add horizontal line at 10%
    geom_hline(yintercept = -10, linetype = "dashed", color = "green") +  # Add horizontal line at 10%
    scale_fill_manual(values = c("positive" = "black", "negative" = "gray")) +
    labs(x = "Cluster", fill = "Index phase", y = "Inverse/direct relation based on mean slope (%)", title = toupper(i)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    theme(legend.position = "none")+
    
    # Add geom_point to add the circles
    geom_point(data = df_long, 
               aes(x = cluster_reorder, y = slope_mean+offset),
               color = df_long$circle_color, 
               fill = df_long$circle_color,
               size = 3, 
               shape = 21, 
               stroke = 2) +
    scale_color_identity()
  
  ggsave(filename = paste0("plot/plot2/fig2_",i,".pdf"), plot = p2, width = 20, height = 10, units = "cm")  
  
}

# Save regressions with new order explained in section 2
#all_data[1,]
#system("plot/plot2/regression/enso_negative_c1_ienso_pcPC1.pdf")


