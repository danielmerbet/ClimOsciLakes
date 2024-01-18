oscillations <- data.frame(year=1901:2019)
#plot with ONI (ENSO) oscillation
oni <- read.table("/home/ry4902/Documents/OscillationsLakes/index/ONI", 
                  header = T)
oni <- aggregate(ANOM ~ YR, data=oni, mean)
oni <- data.frame(year=oni$YR, oni=oni$ANOM)
oscillations <- merge(oscillations,oni, by="year")
#plot with NAO oscillation
nao <- read.table("/home/ry4902/Documents/OscillationsLakes/index/NAO")
nao <- data.frame(year=nao[,1], nao=rowMeans(nao[,2:13]))
oscillations <- merge(oscillations,nao, by="year")
#plot with PDO oscillation
pdo <- read.table("/home/ry4902/Documents/OscillationsLakes/index/PDO", header=T)
pdo <- data.frame(year=pdo[,1], pdo=rowMeans(pdo[,2:13]))
oscillations <- merge(oscillations,pdo, by="year")
#plot with IOD oscillation
iod<- read.table("/home/ry4902/Documents/OscillationsLakes/index/IOD")
iod <- data.frame(year=iod[,1], iod=rowMeans(iod[,2:13]))
oscillations <- merge(oscillations,iod, by="year")

#plot with NPGO oscillation
npgo<- read.table("/home/ry4902/Documents/OscillationsLakes/index/NPGO", header=T)
npgo<- aggregate(NPGO ~ YEAR, data=npgo, FUN=mean)
npgo <- data.frame(year=npgo[,1], npgo=npgo[,2])
oscillations <- merge(oscillations,npgo, by="year")

#plot with NPO oscillation
npo<- read.table("/home/ry4902/Documents/OscillationsLakes/index/NPO")
npo <- data.frame(year=npo[,1], npo=npo[,2])
oscillations <- merge(oscillations,npo, by="year")

write.table(oscillations, file="/home/ry4902/Documents/OscillationsLakes/index/oscillations.txt", quote = F, row.names = F)

