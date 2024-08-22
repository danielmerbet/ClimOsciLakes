library(ncdf4); library(HiClimR); library(raster);library(astsa);library(tidyverse) #library(pcaMethods)

# Set working directory and create necessary folders
variable <- "surftemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/regional/"
directory  <- paste0(directory_base, variable, "/")
print(paste0("############## START:", variable, "##############"))
setwd(directory)
dirs <- c("clustering", "permutation", 
          "pca_04_explained", "pca_04", "ccf_04", "loadings_04",
          "pca_explained", "pca", "ccf", "loadings",
          "acf", "HiClimR_output", "pc_data")
sapply(dirs, dir.create, showWarnings = FALSE)
dir.create("clustering/shp")

# Open NetCDF file and extract data
ncin <- nc_open(paste0(directory_base, "nc_laketemp/gotm_20crv3-era5_obsclim_histsoc_default_", variable,"_global_annual_1901_2021.nc"))
var_nc <- ncvar_get(ncin, variable)
year_range_total <- 1901:2021 #total year range of the data
year_range <- 1950:2021 #subset of years to use (to fit with climate indexes)

# Function to perform PCA and permutation test
pca_eigenperm<- function(data, nperm = 1000){
  pca_out<- prcomp(data, scale. = T)
  eigenperm<- data.frame(matrix(NA, nperm, ncol(data)))
  n<- ncol(data)
  data_i<- data.frame(matrix(NA, nrow(data), ncol(data)))
  for (j in 1: nperm){
    for (i in 1:n){
      data_i[,i]<- sample(data[,i], replace = F)
    }
    pca.perm<- prcomp(data_i, scale. = T)
    eigenperm[j,]<- pca.perm$sdev^2
  }
  colnames(eigenperm)<- colnames(pca_out$rotation)
  eigenperm
  
}


# Extract latitude and longitude
lat <- ncvar_get(ncin, "lat")
lon <- ncvar_get(ncin, "lon")

# Prepare data matrix
mdata <- apply(var_nc,c(3,1),c) #it must be c(3,1) to have a proper estructure in the final matrix
mdata <- apply(mdata,2,c) #points, years
colnames(mdata) <- year_range_total
rownames(mdata) <- paste(sort(rep(lon, length(lat))), rep(lat, length(lon)), sep=",")
# Subset data from 1950 (climate indexes start in this year)
mdata <- mdata[,(50:ncol(mdata))]
#the structure of de data should be lon, lat, e.g., (1st row) -19.75,-39.75 (2nd row) -19.75,-38.75
#see example in TestCases$x[1:10,]

# Set clustering parameters
k <- c(4,4,6,7,3,7)
colors_plot <- c(1:7)

# Read oscillations table
oscillations <- read.table(paste0(directory_base, "clim_index/oscillations.txt"), header = T)

# Calculate confidence interval of the pca
n <- length(oscillations$year)
conf_interval <- 1.96 / sqrt(n)

# Process each region
system(paste0("rm ", directory, "clustering/*.nc"))
c<-0; names_ccf <- c(); ccf_temp_all <- c(); names_ccf_04 <- c(); ccf_temp_all_04 <- c()
pc_importance <- c(); row_names <- c()
index_list <- c(); region_list <- c(); cluster_list <- c(); cor04_list <- c(); pc_list <- c()
for (r in c("Africa","Latinamerica","Asia","Europe","Oceania", "Northern America ")){
  ###1. Clustering starts
  print(paste0("clustering starts for: ", r))
  c<-c+1
  
  if (r=="Latinamerica"){
    y <- HiClimR(mdata, lon = lon, lat = lat, lonStep = 1, latStep = 1, geogMask = TRUE,
                 region = c("South America", "Central America"), meanThresh = NULL, varThresh = NULL, detrend = TRUE,
                 standardize = TRUE, nPC = NULL, method = "ward", hybrid = FALSE, kH = NULL,
                 members = NULL, nSplit = 1, upperTri = TRUE, verbose = TRUE,
                 validClimR = TRUE, k = k[c], minSize = 1, alpha = 0.01,
                 plot = FALSE, colPalette = NULL, hang = -1, labels = FALSE)
  }else if(r=="Northern America "){
    y <- HiClimR(mdata, lon = lon, lat = lat, lonStep = 1, latStep = 1, geogMask = TRUE,
                 region = r, meanThresh = NULL, varThresh = NULL, detrend = TRUE,
                 standardize = TRUE, nPC = NULL, method = "ward", hybrid = FALSE, kH = NULL,
                 members = NULL, nSplit = 1, upperTri = TRUE, verbose = TRUE,
                 validClimR = TRUE, k = k[c], minSize = 1, alpha = 0.01,
                 plot = FALSE, colPalette = NULL, hang = -1, labels = FALSE)
  }else{
    y <- HiClimR(mdata, lon = lon, lat = lat, lonStep = 1, latStep = 1, geogMask = TRUE,
                 continent = r, meanThresh = NULL, varThresh = NULL, detrend = TRUE,
                 standardize = TRUE, nPC = NULL, method = "ward", hybrid = FALSE, kH = NULL,
                 members = NULL, nSplit = 1, upperTri = TRUE, verbose = TRUE,
                 validClimR = TRUE, k = k[c], minSize = 1, alpha = 0.01,
                 plot = FALSE, colPalette = NULL, hang = -1, labels = FALSE)
  }

  # Save HiClimR results to avoid doing it again
  save(y, file=paste0("HiClimR_output/",r,".RData"))
  
  # Save dendogram plot
  pdf(paste0("clustering/",r, "_dendo.pdf"))
  plot(y, hang = -1, labels = FALSE)
  dev.off()

  # Plot the resulting clustering on a map
  coords_df <- data.frame(lon = y$coords[, 1], lat = y$coords[, 2], region = as.factor(y$region))
  ggplot(coords_df, aes(x = lon, y = lat, fill = region)) +
    geom_tile() +
    scale_fill_manual(values = c(colors_plot[1:k[c]])) +
    theme_minimal() +
    labs(title = paste("Region:", r), x = "Longitude", y = "Latitude") +
    coord_fixed()
  ggsave(paste0("clustering/", r, "_map.pdf"),  width = 20, height = 13, units = "cm")
  
  # Export region map and mean timeseries into NetCDF-4 file and raster (tiff)
  y.nc <- HiClimR2nc(y=y, ncfile=paste0("clustering/",r,".nc"), timeunit="years", dataunit="mm")
  y.nctotiff <- brick(paste0("clustering/",r,".nc"), varname="region")
  writeRaster(y.nctotiff, paste0("clustering/",r,".tiff"), overwrite=TRUE)  

  ###1. Clustering finished
 
  # Loop among cluster to implement PCA and permutation analysis 
  for (cluster in as.numeric(names(table(y$region)))){
    
    cluster_mdata <- mdata[(y$region==cluster & (!is.na(y$region))),]
    
    #  Scaling and detrending (lowess method) all the data
    cluster_mdata <- apply(cluster_mdata, MARGIN=1, FUN = function(x){scale(x)})
    cluster_mdata <- apply(cluster_mdata, MARGIN=2, FUN = function(x){detrend(x, lowess = T)})

    ###2. Principal components starts
    print(paste0("PCA for cluster: ", cluster, " region: ",r))
    pca_prcomp <- prcomp(cluster_mdata, center=TRUE, scale=TRUE)
    
    ###3. Permutation test starts
    pca_prcomp_perm<- pca_eigenperm(t(cluster_mdata))
    pca_rand95<- 
      data.frame(Random_Eigenvalues = sapply(pca_prcomp_perm, quantile, 0.95)) %>%
      #95% percentile of randome eigenvalues
      mutate(PC = colnames(pca_prcomp$rotation)) %>%
      #add PC IDs as discrete var
      cbind(Eigenvalues = pca_prcomp$sdev^2)
    #combine rand95 with real eigenvals
    
    ## Only the first 3 PCs will be analysed because they represent most of the variance
    pca_rand95_long<-
      gather(pca_rand95[1:3, ], key = Variable, value = Value, -PC)
    
    # Save permutation test
    ggplot(pca_rand95_long, aes(PC, Value, fill = Variable)) +
      geom_bar(stat = "identity", position = position_dodge())+
      labs(y="Eigenvalue", x="", fill= "") +
      theme_classic()
    ggsave(paste0("permutation/", r, "_cluster",cluster,".pdf"),  width = 20, height = 13, units = "cm")

    ###3. Permutation test finished

    summ <- summary(pca_prcomp)
    #print(summ)
    
    # Save importance (Proportion of Variance) of each pc for each cluster
    pc_importance <- rbind(pc_importance, summ$importance[2,])
    row_names <- c(row_names, paste0(r,cluster))
    
    # Save the first 3 PCs of each cluster
    write.csv(pca_prcomp$x[,1:3], file=paste0("pc_data/", r, cluster,".csv"), quote = F, row.names = F)
    
    variance_explained_3pcs <- summ$importance[2,][1:3]
    raster_variance_explained_pc <- y.nctotiff

    other_clusters <- as.numeric(names(table(y$region)))
    other_clusters <- other_clusters[-cluster]

    raster_variance_explained_pc[raster_variance_explained_pc %in% other_clusters] <- NA
    
    # Analysis made only for the first 3 PCs because they represent most of the variance
    for (pc in 1:3){
      #print("Inside a pc loop")
      if (pca_rand95_long$Value[3+pc]>pca_rand95_long$Value[pc]){
        
        # Significant PCA after permutation test
        scores_plot <- data.frame(year=year_range, score=pca_prcomp$x[,pc])
        df_plot <- merge(scores_plot, oscillations)
        
        # Save al acf for every cluster
        pdf(paste0("acf/",r, "_cluster",cluster,"_pc", pc,".pdf"))
        acf(as.numeric(df_plot$score))
        dev.off()

        for(ind in names(df_plot)[3:8]){
          
          # Cross-correlation between climate indexes and PC
          ccf_temp <- ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5, plot = F)
          
          # Select only the significant ccf values
          if (max(abs(ccf_temp$acf))>conf_interval){
            #print(paste0("Inside a significant ccf value, greater than ", conf_interval))
            
            index_list <- c(index_list, ind)
            region_list <- c(region_list, r)
            cluster_list <- c(cluster_list, cluster)
            pc_list <- c(pc_list, pc)
            if (max(abs(ccf_temp$acf))>0.4){
              cor04_list <- c(cor04_list, "cor04")
            }else{
              cor04_list <- c(cor04_list, as.character(max(abs(ccf_temp$acf))))
            }

            pdf(paste0("ccf/",r, "_cluster",cluster,"_pc", pc, "_", ind, ".pdf"))
            ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5)
            dev.off()
            
            raster_variance_explained_pc[raster_variance_explained_pc==cluster] <- variance_explained_3pcs[pc]
            writeRaster(raster_variance_explained_pc, paste0("pca_explained/",r,cluster,"_pc", pc,"_",ind,".tiff"), overwrite=TRUE)
            
            # Save the ccf plot between de pc and the climate indexes
            pdf(paste0("pca/",r,"_cluster",cluster,"_ccf",pc,"_",ind,".pdf"))
            ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5, ylim=c(-0.5,0.5),ylab=paste0("Correlation PC", pc, " - ", ind), main=paste0("CCF between PC and ", ind," oscillation"))
            dev.off()
            names_ccf_04 <- c(names_ccf_04, paste0(r, "_cluster",cluster,"_ccf", pc, "_", ind))
            ccf_temp_all_04 <- rbind(ccf_temp_all_04, ccf_temp$acf)
            
            # Extract the loadings of each pc for each cluster
            split_data <- strsplit(names(pca_prcomp$rotation[, 1]), ",")
            loadings_cluster <- data.frame(matrix(unlist(split_data), ncol=2, byrow=TRUE))
            names(loadings_cluster) <- c("lon", "lat")
            loadings_cluster$lat <- as.numeric(loadings_cluster$lat)
            loadings_cluster$lon <- as.numeric(loadings_cluster$lon)
            loadings_cluster$loading <- as.numeric(pca_prcomp$rotation[, pc])
            
            #plot the loadings
            ggplot(loadings_cluster, aes(x = lon, y = lat, fill = loading)) +
              geom_tile() +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
              theme_minimal() +
              labs(title = paste("PC", pc, "Loadings for Cluster", cluster, "in", r), x = "Longitude", y = "Latitude") +
              coord_fixed()
            ggsave(paste0("loadings/", r, "_cluster", cluster, "_pc", pc, "_", ind, ".pdf"), width = 20, height = 13, units = "cm")
          }
          
          if (max(abs(ccf_temp$acf))>0.4){
            print("Inside a ccf value greater than 0.4")
            
            pdf(paste0("ccf_04/",r, "_cluster",cluster,"_pc", pc, "_", ind, ".pdf"))
            ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5)
            dev.off()

            raster_variance_explained_pc[raster_variance_explained_pc==cluster] <- variance_explained_3pcs[pc]
            writeRaster(raster_variance_explained_pc, paste0("pca_04_explained/",r,cluster,"_pc", pc,"_",ind,".tiff"), overwrite=TRUE)
            
            # Save the ccf plot between de pc and the climate indexes
            pdf(paste0("pca_04/",r,"_cluster",cluster,"_ccf",pc,"_",ind,".pdf"))
            ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5, ylim=c(-0.5,0.5),ylab=paste0("Correlation PC", pc, " - ", ind), main=paste0("CCF between PC and ", ind," oscillation"))
            dev.off()
            names_ccf_04 <- c(names_ccf_04, paste0(r, "_cluster",cluster,"_ccf", pc, "_", ind))
            ccf_temp_all_04 <- rbind(ccf_temp_all_04, ccf_temp$acf)
            
            # Extract the loadings of each pc for each cluster
            split_data <- strsplit(names(pca_prcomp$rotation[, 1]), ",")
            loadings_cluster <- data.frame(matrix(unlist(split_data), ncol=2, byrow=TRUE))
            names(loadings_cluster) <- c("lon", "lat")
            loadings_cluster$lat <- as.numeric(loadings_cluster$lat)
            loadings_cluster$lon <- as.numeric(loadings_cluster$lon)
            loadings_cluster$loading <- as.numeric(pca_prcomp$rotation[, pc])
            
            #plot the loadings
            ggplot(loadings_cluster, aes(x = lon, y = lat, fill = loading)) +
              geom_tile() +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
              theme_minimal() +
              labs(title = paste("PC", pc, "Loadings for Cluster", cluster, "in", r), x = "Longitude", y = "Latitude") +
              coord_fixed()
            ggsave(paste0("loadings_04/", r, "_cluster", cluster, "_pc", pc, "_", ind, ".pdf"), width = 20, height = 13, units = "cm")
          }

        }
        
      }else{
        print(paste0("NON SIGNIFICANT PCA ADDED", r, "_cluster",cluster,"_ccf", pc))        
      }
      
      
    }
  }
}

save(names_ccf_04, ccf_temp_all_04, file="pca_04/ccf_04.RData")
ccf_temp_all_04 <- t(ccf_temp_all_04)
colnames(ccf_temp_all_04) <- names_ccf_04
write.csv(ccf_temp_all_04, file="pca_04/ccf_04.csv", quote=F, row.names=F)

row.names(pc_importance) <- row_names
write.csv(pc_importance, file = "pc_data/pc_proportion_variance.csv", quote = F)

summary_ccf <- data.frame(index=index_list, region=region_list, cluster=cluster_list, pc=pc_list, cor=cor04_list)
write.csv(summary_ccf, file = "summary_ccf.csv", quote = F, row.names = F)

print(paste0("############## END:", variable, "##############"))
