a<-read.csv("~/Documents/ClimOsciLakes/bottemp/pc_data/pc_proportion_variance.csv")
a
pdf("test.pdf")
plot(a[1,2:73])
dev.off()

for (i in 1:31){print(cumsum(as.numeric(a[i,2:5])))}

a<-read.csv("~/Documents/ClimOsciLakes/surftemp/pc_data/pc_proportion_variance.csv")
a
pdf("test.pdf")
plot(a[1,2:73])
dev.off()

for (i in 1:31){print(cumsum(as.numeric(a[i,2:5])))}

b <- read.csv("~/Documents/ClimOsciLakes/bottemp/summary_ccf.csv")
