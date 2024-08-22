
# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"
directory  <- paste0(directory_base, "data/",variable)
k <- 47
setwd(directory)

index_var <- read.csv("XGBoost/index_variable.csv")
stable1 <- read.csv(paste0("../../supplementary/STable1_",variable,".csv"))

# Function to convert a single string to a numeric vector
convert_ranges_to_vector <- function(range_string) {
  parts <- strsplit(range_string, ",")[[1]]
  numeric_vector <- unlist(lapply(parts, function(x) eval(parse(text = gsub("-", ":", x)))))
  return(numeric_vector)
}

# Apply the function to each row in the column
clusters_num <- lapply(stable1$Clusters, convert_ranges_to_vector)

stable1_new <- data.frame()
for (r in 1:nrow(stable1)){
  temp_df <- data.frame(cluster_reorder=unlist(clusters_num[r]),
                        zone=stable1[r,]$Climate.zones,
                        total=stable1[r,]$Total)
  stable1_new <- rbind(stable1_new, temp_df)
}

all_data <- merge(stable1_new, index_var, by="cluster_reorder")

write.csv(all_data, file=paste0("../../supplementary/STable1_",variable,"_variable.csv"))

