# Load necessary libraries

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"

print(paste0("############## START:", variable, "##############"))

year_range_total <- 1901:2021 #total year range of the data
year_range <- 1950:2021 #subset of years to use (to fit with climate indexes)
meteo_vars <- c("tas", "hurs", "pr", "ps", "rsds", "rlds", "sfcwind", "tasmax", "tasmin")
region <- "global"


print(paste0("############## START: k", k, "##############"))

directory  <- paste0(directory_base, "data/",variable)
setwd(directory)
dir.create("XGBoost/")
dir.create("XGBoost/merge_table")

summary_file <- read.csv("summary_ccf_global.csv")

summary_file <- summary_file %>% distinct(cluster, pc, .keep_all = TRUE)

k <- 47

# List of all significant ccf
for (combi in 1:nrow(summary_file)){
  
  ncin_region <- nc_open(paste0("clustering/", region, ".nc"))
  var_nc_region <- ncvar_get(ncin_region, "region")
  
  #clusters_region <- sort(unique(as.vector(var_nc_region)))
  cluster <- summary_file$cluster[combi]
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
  #pc_value <- strsplit(ccf_sign_temp[[1]][3], "pc")[[1]][2]
  pc_value <- summary_file$pc[combi]
  pc_data[c(paste0("PC", pc_value), "year")]
  
  # Merge PC scores with meteorological data
  # Assuming both data frames have a common 'year' column
  df <- merge(pc_data, meteorological_data, by = "year")
  #Detrend the data
  df_det <- c(); df_det$year <- year_range; df_det <- data.frame(df_det)
  for (i in names(df)[2:ncol(df)]){
    df_det[i] <- as.vector(scale(detrend(df[i], lowess = T)))
  }
  
  write.csv(df, file=paste0("XGBoost/merge_table/", region, cluster,"_", "PC", pc_value, ".csv"))
}



print(paste0("############## END:", variable, "##############"))
