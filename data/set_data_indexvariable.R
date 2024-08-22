
# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
directory  <- paste0(directory_base, "data/",variable)
setwd(directory)
k <- 47

load("clustering/order_to_apply.RData")
sum_xgboost <- read.csv("XGBoost/save_data.csv")
clim_index <- read.csv("../../clim_index/oscillations.txt", sep= " ")
colnames(clim_index) <- c("year", "enso",  "nao",  "pdo",  "iod",  "npgo", "npo")

tel <- read.csv(paste0("../../supplementary/teleconnections_", variable,".csv"))
n <- length(clim_index$year)
conf_interval <- 1.96 / sqrt(n)
cor_save <- c(); cluster_save <- c(); index_save <- c(); var_save <- c();recluster_save <- c();
for (r in 1:nrow(sum_xgboost)){
  xgb_temp <- sum_xgboost[r,]
  pc_data <- read.csv(paste0("XGBoost/merge_table/","global",xgb_temp$cluster,"_PC",xgb_temp$pc,".csv"))
  #pc_data <- pc_data[paste0("PC", xgb_temp$pc)]
  for (index in unique(tel[tel$cluster==xgb_temp$cluster,]$index)){
    ccf_temp <- ccf(clim_index[index], pc_data[xgb_temp$var])
    if(max(abs(ccf_temp$acf))>conf_interval){
      cor_save <- c(cor_save, max(abs(ccf_temp$acf)))
      cluster_save <- c(cluster_save, xgb_temp$cluster)
      index_save <- c(index_save, index)
      var_save <- c(var_save, xgb_temp$var)
    }else{
      cor_save <- c(cor_save, NA)
      cluster_save <- c(cluster_save, xgb_temp$cluster)
      index_save <- c(index_save, index)
      var_save <- c(var_save, xgb_temp$var)
    }
    
  }
  
}

index_var <- data.frame(cluster_save,index_save,var_save,cor_save)
index_var$cluster_reorder <- NA; c <- 0
for (l in order_to_apply){
  c <- c+1
  index_var$cluster_reorder[index_var$cluster_save==l] <- c
}

write.csv(index_var, file = "XGBoost/index_variable.csv", row.names = F, quote = F)

