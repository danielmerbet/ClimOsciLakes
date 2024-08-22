library(WaveletComp)

# Set working directory and create necessary folders
variable <- "bottemp" #surftemp or bottemp
directory_base <- "/home/dmercado/Documents/ClimOsciLakes/"

print(paste0("############## START:", variable, "##############"))

# Read oscillations table
oscillations <- read.table(paste0(directory_base, "clim_index/oscillations.txt"), header = T)

for (k in 5:50){
  
  print(paste0("############## START: k", k, "##############"))
  
  directory  <- paste0(directory_base, "global/",variable, "/k", k)
  setwd(directory)
  dir.create("wavelet")
  
  # read summary of PC analysis
  summary_pca <- read.csv("summary_ccf_global.csv")
  
  #Subset only the significant ccf between PC and climate index
  ccf_sign <- list.files(path = "ccf_04/")
  wavelet_results <- list()
  for (relation in 1:length(ccf_sign)){
    ccf_sign_temp <- strsplit(ccf_sign, "_")[[relation]]
    r <- paste0(ccf_sign_temp[1], as.numeric(strsplit(ccf_sign_temp[2], "cluster")[[1]][2]))
    
    # Read pca analysis table
    #list.files(path = "pc_data/")
    #r <- "Africa1.csv"
    pca_comp <- read.csv(paste0("pc_data/", r, ".csv"))
    year_range <- 1950:2021 #subset of years to use (to fit with climate indexes)
    pca_comp <- data.frame(year=year_range, pca_comp)
    
    # Merge pca with oscillation data
    pca_osci <- merge(pca_comp, oscillations, by.x = "year", by.y = "year", all = TRUE)
    
    #PC
    pc <- paste0("PC", as.numeric(strsplit(ccf_sign_temp[3], "pc")[[1]][2]))
    
    #index
    index <- strsplit(ccf_sign_temp[4], ".pdf")[[1]]
    
    # Wavelet Decomposition
    #for (pc in colnames(pca_comp)[-1]) {
    #  for (index in colnames(oscillations)[-1]) {  # Assuming the first column is 'year'
    combined_data <- data.frame(time = as.numeric(pca_osci$year), pc = as.numeric(pca_osci[[pc]]), index = as.numeric(pca_osci[[index]]))
    #combined_data <- combined_data[complete.cases(combined_data), ]  # Remove rows with NA values
    
    # Perform wavelet analysis
    #wavelet_pc <- analyze.wavelet(combined_data, "pc", loess.span = 0.75)
    #wavelet_index <- analyze.wavelet(combined_data, "index", loess.span = 0.75)
    
    # Cross-wavelet power
    cross_wavelet <- analyze.coherency(combined_data, my.pair = c("pc", "index"), loess.span = 0.75)
    
    pdf(paste0("wavelet/", r, "_", pc, "_", index, ".pdf"))
    wt.image(cross_wavelet,  
             main = paste("Cross-Wavelet Power Spectrum of ", r, "for ", pc, "and", index),
             timelab = "time")
    dev.off()
    
    #wavelet_results[[paste0(pc, "_", index)]] <- list(wavelet_pc = wavelet_pc, wavelet_index = wavelet_index, cross_wavelet = cross_wavelet)
    wavelet_results[[paste0(pc, "_", index)]] <- list(cross_wavelet = cross_wavelet)
    #  }
    #}
    
  }
  
}

print(paste0("############## END:", variable, "##############"))

