#Select the best k possible

library(raster)

#SURFTEMP

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"

#Extract de cluster with the minimum amount of lakes (pixels) for all k values
min_values <- c(); max_values <- c()
for (k in 5:100){
  directory  <- paste0(directory_base, "global/",variable, "/k", k)
  clustering <- raster(paste0(directory,"/clustering/global.tiff"))
  table(clustering[])
  min_values <- c(min_values, min(table(clustering[])))
  max_values <- c(max_values, max(table(clustering[])))
}


nonunique_cluster_pc <- c()
unique_cluster_pc <- c()
unique_PC1 <- c()
unique_PC1_04 <- c()
unique_PC1_noNPGO <- c()
unique_PC1_04_noNPGO <- c()
for (k in 5:100){
  
  directory  <- paste0(directory_base, "global/",variable, "/k", k)
  setwd(directory)
  
  #temporal while running data again
  #rl <- length(list.files(path="ccf/"))
  summary_ccf <- read.csv("summary_ccf_global.csv")
  #summary_ccf <- summary_ccf[1:(nrow(summary_ccf)/6),]
  #summary_ccf <- summary_ccf[,c("index","cluster","pc","cor")]
  #temporal while running data again
  
  #summary for significant correlation>04 between PCs of clusters and climate oscillation (index)
  summ_04 <- summary_ccf[summary_ccf$cor>=0.4,]
  #nonunique_cluster_pc <- c(nonunique_cluster_pc, nrow(summ_04))
  #unique_cluster_pc <- c(unique_cluster_pc, nrow(unique(summ_04[c("cluster", "pc")])))
  
  
  #summary significant correlation with PC1 (strongest PC)
  summary_PC1 <- summary_ccf[summary_ccf$pc==1,]
  
  #summary for all
  unique_PC1 <- c(unique_PC1, nrow(unique(summary_PC1[c("index","cluster", "pc")])))
  
  #summary for correlacion>0.4
  summary_PC1_04 <- summary_PC1[summary_PC1$cor>=0.4,]
  unique_PC1_04 <- c(unique_PC1_04, nrow(unique(summary_PC1_04[c("index","cluster","pc")])))
  
  #remove NPGO and NPO
  summary_PC1_noNPGO <- summary_PC1[summary_PC1$index!="npgo",]
  summary_PC1_noNPGO <- summary_PC1_noNPGO[summary_PC1_noNPGO$index!="npo",]
  summary_PC1_04_noNPGO <- summary_PC1_04[summary_PC1_04$index!="npgo",]
  summary_PC1_04_noNPGO <- summary_PC1_04_noNPGO[summary_PC1_04_noNPGO$index!="npo",]
  unique_PC1_noNPGO <- c(unique_PC1_noNPGO, nrow(unique(summary_PC1_noNPGO[c("index","cluster", "pc")])))
  unique_PC1_04_noNPGO <- c(unique_PC1_04_noNPGO, nrow(unique(summary_PC1_04_noNPGO[c("index","cluster", "pc")])))
  
}



#Select best k
best_k_surf <- data.frame(min=min_values, #dif_04=(nonunique_cluster_pc- unique_cluster_pc), 
                     unique_PC1,unique_PC1_04,cluster=5:100,
                     unique_PC1_04/unique_PC1, 
                     unique_PC1_noNPGO,
                     unique_PC1_04_noNPGO,
                     unique_PC1_04_noNPGO/unique_PC1_noNPGO)

#BOTTEMP


# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"

#Extract de cluster with the minimum amount of lakes (pixels) for all k values
min_values <- c()
for (k in 5:100){
  directory  <- paste0(directory_base, "global/",variable, "/k", k)
  clustering <- raster(paste0(directory,"/clustering/global.tiff"))
  table(clustering[])
  min_values <- c(min_values, min(table(clustering[])))
}


nonunique_cluster_pc <- c()
unique_cluster_pc <- c()
unique_PC1 <- c()
unique_PC1_04 <- c()
unique_PC1_noNPGO <- c()
unique_PC1_04_noNPGO <- c()
for (k in 5:100){
  
  directory  <- paste0(directory_base, "global/",variable, "/k", k)
  setwd(directory)
  
  #temporal while running data again
  #rl <- length(list.files(path="ccf/"))
  summary_ccf <- read.csv("summary_ccf_global.csv")
  #summary_ccf <- summary_ccf[1:(nrow(summary_ccf)/6),]
  #summary_ccf <- summary_ccf[,c("index","cluster","pc","cor")]
  #temporal while running data again
  
  #summary for significant correlation>04 between PCs of clusters and climate oscillation (index)
  summ_04 <- summary_ccf[summary_ccf$cor>=0.4,]
  #nonunique_cluster_pc <- c(nonunique_cluster_pc, nrow(summ_04))
  #unique_cluster_pc <- c(unique_cluster_pc, nrow(unique(summ_04[c("cluster", "pc")])))
  
  
  #summary significant correlation with PC1 (strongest PC)
  summary_PC1 <- summary_ccf[summary_ccf$pc==1,]
  
  #summary for all
  unique_PC1 <- c(unique_PC1, nrow(unique(summary_PC1[c("index","cluster", "pc")])))
  
  #summary for correlacion>0.4
  summary_PC1_04 <- summary_PC1[summary_PC1$cor>=0.4,]
  unique_PC1_04 <- c(unique_PC1_04, nrow(unique(summary_PC1_04[c("index","cluster","pc")])))
  
  #remove NPGO and NPO
  summary_PC1_noNPGO <- summary_PC1[summary_PC1$index!="npgo",]
  summary_PC1_noNPGO <- summary_PC1_noNPGO[summary_PC1_noNPGO$index!="npo",]
  summary_PC1_04_noNPGO <- summary_PC1_04[summary_PC1_04$index!="npgo",]
  summary_PC1_04_noNPGO <- summary_PC1_04_noNPGO[summary_PC1_04_noNPGO$index!="npo",]
  unique_PC1_noNPGO <- c(unique_PC1_noNPGO, nrow(unique(summary_PC1_noNPGO[c("index","cluster", "pc")])))
  unique_PC1_04_noNPGO <- c(unique_PC1_04_noNPGO, nrow(unique(summary_PC1_04_noNPGO[c("index","cluster", "pc")])))
  
}



#Select best k
best_k_bot <- data.frame(min=min_values, #dif_04=(nonunique_cluster_pc- unique_cluster_pc), 
                     unique_PC1,unique_PC1_04,cluster=5:100,
                     unique_PC1_04/unique_PC1, 
                     unique_PC1_noNPGO,
                     unique_PC1_04_noNPGO,
                     unique_PC1_04_noNPGO/unique_PC1_noNPGO)


#PLOT
# Plot the first dataset
setwd(directory_base)
pdf("figures/sfig1.pdf")
plot(best_k_surf$cluster,best_k_surf$unique_PC1_04_noNPGO.unique_PC1_noNPGO, xlab= "k", ylab="%04/total")
points(best_k_bot$cluster, best_k_bot$unique_PC1_04_noNPGO.unique_PC1_noNPGO, col="red")
# Add the second dataset
par(new = TRUE)
plot(best_k_surf$cluster, best_k_surf$min, type="l",  axes = FALSE, ylab=" ", xlab=" ")
points(best_k_bot$cluster, best_k_bot$min, type="l", col="red")
axis(side = 4)
mtext("Minimum lakes in cluster", side = 4, line = 3)
abline(v=47, col="green")
dev.off()

# The minimum number in cluster
print(best_k_surf$min[43])
print(best_k_bot$min[43])

