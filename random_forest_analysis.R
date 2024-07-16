# Load necessary libraries
library(randomForest)
library(tidyverse)
library(ncdf4)
library(astsa)

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
directory  <- paste0(directory_base, variable, "/")
print(paste0("############## START:", variable, "##############"))
setwd(directory)
dir.create("randomforest")

# Function to train and evaluate a random forest model for a given PC
train_random_forest <- function(pc) {
  formula <- as.formula(paste(pc, "~ . -year"))
  rf_model <- randomForest(formula, data = df, importance = TRUE)
  
  # Plot importance of variables
  importance <- importance(rf_model)
  var_importance <- data.frame(Variables = rownames(importance), Importance = importance[,1])
  var_importance <- var_importance[order(var_importance$Importance, decreasing = TRUE), ]
  
  print(rf_model)
  varImpPlot(rf_model, main = paste("Variable Importance for", pc))
  
  return(list(model = rf_model, importance = var_importance))
}

year_range_total <- 1901:2021 #total year range of the data
year_range <- 1950:2021 #subset of years to use (to fit with climate indexes)
meteo_vars <- c("tas", "hurs", "pr", "ps", "rsds", "rlds", "sfcwind", "tasmax", "tasmin")

# List of all significant ccf
for (ccf in list.files(path = "ccf_04/")){
  ccf_sign_temp <- strsplit(ccf, "_")
  
  region <- ccf_sign_temp[[1]][1]
  ncin_region <- nc_open(paste0("clustering/", region, ".nc"))
  var_nc_region <- ncvar_get(ncin_region, "region")
  
  #clusters_region <- sort(unique(as.vector(var_nc_region)))
  cluster <- as.numeric(strsplit(ccf_sign_temp[[1]][2], "cluster")[[1]][2])
  #for (cluster in clusters_region){
  
  # Find positions where var_nc_region is equal to the number of cluster
  positions <- which(var_nc_region == cluster, arr.ind = TRUE)

  values_mean <- list()
  values_mean[["year"]] <- year_range_total
  
  for (meteo_var in meteo_vars){
    
    # Open NetCDF file and extract data
    ncin <- nc_open(paste0(directory_base, "clim_data/20crv3-era5_obsclim_", meteo_var,"_global_yearly_1901_2021.nc"))
    var_nc <- ncvar_get(ncin, meteo_var)
    
    # Loop through the positions and extract values from var_nc
    values_list <- list()
    for (i in 1:nrow(positions)) {
      row <- positions[i, 1]
      col <- positions[i, 2]
      values_list[[i]] <- var_nc[row, col, ]
    }
    
    # Convert the list to a matrix (or array) if needed
    values_matrix <- do.call(rbind, values_list)
    values_mean[[meteo_var]] <- apply(values_matrix, 2, mean)
    
  }
  
  meteorological_data <- data.frame(values_mean)
  pc_data <- read.csv(paste0("pc_data/", region, cluster, ".csv"))
  pc_data$year <- year_range
  pc_value <- strsplit(ccf_sign_temp[[1]][3], "pc")[[1]][2]
  pc_data[c(paste0("PC", pc_value), "year")]
  
  # Merge PC scores with meteorological data
  # Assuming both data frames have a common 'year' column
  df <- merge(pc_data, meteorological_data, by = "year")
  #Detrend the data
  df_det <- c(); df_det$year <- year_range; df_det <- data.frame(df_det)
  for (i in names(df)[2:ncol(df)]){
    df_det[i] <- as.vector(scale(detrend(df[i], lowess = T)))
  }
  
  # Train Random Forest model for each principal component
  set.seed(123)  # For reproducibility
  
  # Train models for the first three PCs
  results_pc1 <- train_random_forest("PC1")
  results_pc2 <- train_random_forest("PC2")
  results_pc3 <- train_random_forest("PC3")
  
  # Display the variable importance
  print("Variable Importance for PC1:")
  print(results_pc1$importance)
  print("Variable Importance for PC2:")
  print(results_pc2$importance)
  print("Variable Importance for PC3:")
  print(results_pc3$importance)
  
  importance_rf <- data.frame(pc1=results_pc1$importance, pc2=results_pc2$importance, pc3=results_pc3$importance)
  
  write.csv(importance_rf, file=paste0("randomforest/", region, cluster,"_", "PC", pc_value, ".csv"), row.names = F, quote = F)
  
  #}
}

print(paste0("############## END:", variable, "##############"))
