library(xgboost)

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
directory  <- paste0(directory_base, "data/",variable)
setwd(directory)
dir.create("XGBoost/result")

print(paste0("############## START:", variable, "##############"))

meteo_vars <- c("tas", "hurs", "pr", "ps", "rsds", "rlds", "sfcwind", "tasmax", "tasmin")

summary_file <- read.csv("summary_ccf_global.csv")
summary_file <- summary_file %>% distinct(cluster, pc, .keep_all = TRUE)

predictive_perc <- 0.25

var_save <- c(); predict_save <- c(); cluster_save <- c(); pc_save <- c()
# Loop for all cases in summary
for (combi in 1:nrow(summary_file)){
  
  #Load data
  cluster <- summary_file$cluster[combi]
  pc_value <- summary_file$pc[combi]
  pc <- paste0("PC", pc_value)
  data_temp <- read.csv(paste0("XGBoost/merge_table/global",cluster, 
                          "_PC", pc_value,".csv"))
  
  # Prepare data
  x <- as.matrix(data_temp[meteo_vars])  # Exclude the target variable (mpg)
  y <- unlist(data_temp[pc])
  
  # Convert data to DMatrix format
  dtrain <- xgb.DMatrix(data = x, label = y)
  
  # Set parameters for XGBoost
  params <- list(
    objective = "reg:squarederror",
    max_depth = 3,
    eta = 0.1,
    nthread = 2
  )
  
  # Train the XGBoost model
  set.seed(42)
  xgb_model <- xgboost(
    data = dtrain,
    params = params,
    nrounds = 100,
    verbose = 0
  )
  
  # Get feature importance scores
  importance_matrix <- xgb.importance(feature_names = colnames(x), model = xgb_model)
  
  # Print the importance matrix
  #print(importance_matrix)
  
  # Plot feature importance
  #xgb.plot.importance(importance_matrix)
  
  # Save data
  write.csv(importance_matrix, file=paste0("XGBoost/result/global", cluster,"_", "PC", pc_value, ".csv"))
  
  if (max(importance_matrix$Gain)>predictive_perc){
    #print(importance_matrix$Feature[importance_matrix$Gain>predictive_perc])
    #print(importance_matrix$Gain[importance_matrix$Gain>predictive_perc])
    var_save <- c(var_save, importance_matrix$Feature[importance_matrix$Gain>predictive_perc][1]) 
    predict_save <- c(predict_save, importance_matrix$Gain[importance_matrix$Gain>predictive_perc][1])
    cluster_save <- c(cluster_save, cluster)
    pc_save <- c(pc_save, pc_value)
  }
}

save_all <- data.frame(cluster=cluster_save, pc=pc_save, var=var_save, predict=predict_save)

write.csv(save_all, file="XGBoost/save_data.csv", row.names = F)