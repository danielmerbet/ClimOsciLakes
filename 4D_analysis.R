library(ncdf4); library(HiClimR); library(raster);library(astsa);library(tidyverse) #library(pcaMethods)

#open file data
setwd("/data/brussel/vo/000/bvo00012/vsc10623/OsciLakes")
dir.create("clustering")
dir.create("permutation")
dir.create("pca_04_explained")
dir.create("pca_04")
ncin <- nc_open("nc_laketemp/gotm_ensmean_obsclim_histsoc_default_surftemp_global_daily_1901_2021.nc")
var_nc <- ncvar_get(ncin, "surftemp")
year_range_total <- 1901:2021 #total year range of the data
year_range <- 1950:2021 #subset of years to use

#functions to use
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
#

lat <- ncvar_get(ncin, "lat")
lon <- ncvar_get(ncin, "lon")

mdata <- apply(var_nc,c(3,1),c) #it must be c(3,1) to have a proper estructure in the final matrix
mdata <- apply(mdata,2,c) #points, years
colnames(mdata) <- year_range_total
rownames(mdata) <- paste(sort(rep(lon, length(lat))), rep(lat, length(lon)), sep=",")

#subset data from 1950 (climate indexes start in this year)
mdata <- mdata[,(50:ncol(mdata))]

#the structure of de data should be lon, lat, e.g., (1st row) -19.75,-39.75 (2nd row) -19.75,-38.75
#see example in TestCases$x[1:10,]
#k<-c(4,3,4,5,3,4)
k <- c(4,4,6,7,3,7)
colors_plot <- c(1:7)

#read oscillations table
oscillations <- read.table("clim_index/oscillations.txt", header = T)

## Colors for gridded data
#colPalette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
c<-0
names_ccf <- c(); ccf_temp_all <- c()
names_ccf_04 <- c(); ccf_temp_all_04 <- c()
#for (r in geogMask()$continent[-3]){
for (r in c("Africa","Americas","Asia","Europe","Oceania", "Northern America ")){
  ###1. clustering starts
  print(paste0("clustering starts for: ", r))
  c<-c+1
  
  if (r=="Americas"){
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

  
  #save plot
  
  pdf(paste0("clustering/",r, "_dendo.pdf"))
  plot(y, hang = -1, labels = FALSE)
  dev.off()
  #plot(y$coords[, 1], y$coords[, 2], col = colPalette(max(4, na.rm = TRUE)),
  #     pch = 15, cex = 1)
  pdf(paste0("clustering/",r, ".pdf"))
  RegionsMap <- matrix(y$region, nrow = length(unique(y$coords[, 1])), byrow = TRUE)
  image(unique(y$coords[, 1]), unique(y$coords[, 2]), RegionsMap, col = colors_plot[1:k[c]])
  legend("bottomleft", legend=colors_plot[1:k[c]], col=colors_plot[1:k[c]], pch=rep(15,k[c]), cex=0.5, ncol=2)
  dev.off()
  
  ## Export region map and mean timeseries into NetCDF-4 file and raster (tiff)
  y.nc <- HiClimR2nc(y=y, ncfile=paste0("clustering/",r,".nc"), timeunit="years", dataunit="mm")
  y.nctotiff <- brick(paste0("clustering/",r,".nc"), varname="region")
  writeRaster(y.nctotiff, paste0("clustering/",r,".tiff"))  

  ###clustering finished
 
  #loop among cluter to implement PCA and permutation analysis 
  for (cluster in as.numeric(names(table(y$region)))){
    
    cluster_mdata <- mdata[(y$region==cluster & (!is.na(y$region))),]
    
    #scaling and detrending (lowess method) all the data
    cluster_mdata <- apply(cluster_mdata, MARGIN=1, FUN = function(x){scale(x)})
    cluster_mdata <- apply(cluster_mdata, MARGIN=2, FUN = function(x){detrend(x, lowess = T)})

    ###2. Principal components starts
    print(paste0("PCA for cluster: ", cluster, r))
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
    
    ## only the first 3 PCs
    pca_rand95_long<-
      gather(pca_rand95[1:3, ], key = Variable, value = Value, -PC)
    
    #pdf(paste0("permutation/", r, "_cluster",cluster,".pdf"))
    ggplot(pca_rand95_long, aes(PC, Value, fill = Variable)) +
      geom_bar(stat = "identity", position = position_dodge())+
      labs(y="Eigenvalue", x="", fill= "") +
      theme_classic()
    ggsave(paste0("permutation/", r, "_cluster",cluster,".pdf"),  width = 20, height = 13, units = "cm")
    #dev.off()
    ###3. Permutation test finished
    

    summ <- summary(pca_prcomp)
    #print(summ)
    
    variance_explained_3pcs <- summ$importance[2,][1:3]
    raster_variance_explained_pc <- y.nctotiff

    other_clusters <- as.numeric(names(table(y$region)))
    other_clusters <- other_clusters[-cluster]

    raster_variance_explained_pc[raster_variance_explained_pc %in% other_clusters] <- NA
    
    #we will only work with the first 3 PCs
    for (pc in 1:3){
      print("entro a pc loop")
      if (pca_rand95_long$Value[3+pc]>pca_rand95_long$Value[pc]){
        
        #Significant PCA after permutation test
        scores_plot <- data.frame(year=year_range, score=pca_prcomp$x[,pc])
        df_plot <- merge(scores_plot, oscillations)
        
        pdf(paste0("acf/",r, "_cluster",cluster,"_pc", pc,".pdf"))
        acf(as.numeric(df_plot$score))
        dev.off()

        for(ind in names(df_plot)[3:8]){
          #cross-correlation between climate indexes and PC
          ccf_temp <- ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5)
          if (max(abs(ccf_temp$acf))>0.3){

            raster_variance_explained_pc[raster_variance_explained_pc==cluster] <- variance_explained_3pcs[pc]
            #writeRaster(raster_variance_explained_pc, paste0("pca_explained/",r,cluster,"_pc", pc,"_",ind,".tiff"))

            #pdf(paste0("pca/",r, "_cluster",cluster,"_ccf", pc,"_",ind,".pdf"))
            #ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5, ylim=c(-0.5,0.5),ylab=paste0("Correlation PC", pc, " - ", ind), main=paste0("CCF between PC ", pc," and ", ind," oscillation"))
            #dev.off()
            names_ccf <- c(names_ccf, paste0(r, "_cluster", cluster,"_ccf", pc, "_", ind))
            ccf_temp_all <- rbind(ccf_temp_all, ccf_temp$acf)
          }
          if (max(abs(ccf_temp$acf))>0.4){
           print("entro a 0.4")

            raster_variance_explained_pc[raster_variance_explained_pc==cluster] <- variance_explained_3pcs[pc]
            writeRaster(raster_variance_explained_pc, paste0("pca_04_explained/",r,cluster,"_pc", pc,"_",ind,".tiff"))

            pdf(paste0("pca_04/",r,"_cluster",cluster,"_ccf",pc,"_",ind,".pdf"))
            ccf(as.numeric(df_plot$score),df_plot[ind], lag.max = 5, ylim=c(-0.5,0.5),ylab=paste0("Correlation PC", pc, " - ", ind), main=paste0("CCF between PC and ", ind," oscillation"))
            dev.off()
            names_ccf_04 <- c(names_ccf_04, paste0(r, "_cluster",cluster,"_ccf", pc, "_", ind))
            ccf_temp_all_04 <- rbind(ccf_temp_all_04, ccf_temp$acf)
          }
        


        }
        
      }else{
        print(paste0("NON SIGNIFICANT PCA ADDED", r, "_cluster",cluster,"_ccf", pc))        
      }
      
      
    }
    
  }

  
}

#save(names_ccf, ccf_temp_all, file="pca/ccf.RData")
save(names_ccf_04, ccf_temp_all_04, file="pca_04/ccf_04.RData")

#ccf_temp_all <- t(ccf_temp_all)
ccf_temp_all_04 <- t(ccf_temp_all_04)

colnames(ccf_temp_all_04) <- names_ccf_04
#colnames(ccf_temp_all) <- names_ccf

#write.csv(ccf_temp_all, file="pca/ccf.csv", quote=F, row.names=F)
write.csv(ccf_temp_all_04, file="pca_04/ccf_04.csv", quote=F, row.names=F)
