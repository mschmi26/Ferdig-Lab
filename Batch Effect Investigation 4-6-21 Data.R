setwd("/Users/kbuttons/Documents/Research/P01/Analysis/MKK2835xNHP1337 Drug Data/")

Individual_Drug_Replicates <- read.csv("Notre Dame individual replicates 032821_fixed.csv",header=TRUE,as.is=TRUE)
BatchInfo <- read.csv("Batch Info.csv",as.is=TRUE,header=TRUE)

DrugNames <- unique(Individual_Drug_Replicates[,"Sample.ID"])
DrugNames <- cbind(DrugNames,DrugNames)
for(i in 1:length(DrugNames[,1])){
  DrugRows <- which(Individual_Drug_Replicates[,"Sample.ID"]==DrugNames[i,1])
  if(Individual_Drug_Replicates[DrugRows[1],"Sample.Name"]!=""){
    DrugNames[i,2] <- Individual_Drug_Replicates[which(Individual_Drug_Replicates[,"Sample.ID"]==DrugNames[i,1])[1],"Sample.Name"]
  }
}

### Transform to IC50 and AUC matrices with rows as individual reps and collumns as drugs
library(reshape2)


DrugRepWAC50 <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="AC50..uM.")
DrugRepWAUC <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="AUC")
DrugRepWCC <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="CC.v2")

### Screen any data from curves that are not CURVE_CLASS -1.1, -1.2 or -2.1
DrugRepWAC50filt <- DrugRepWAC50
DrugRepWAUCfilt <- DrugRepWAUC
DrugRepWAC50filt[,c(-1,-2)] <- NA
DrugRepWAUCfilt[,c(-1,-2)] <- NA
for(i in 3:length(DrugRepWAC50filt)){
  DrugRepWAC50filt[which(DrugRepWCC[,i]=="-1.1"),i] <- DrugRepWAC50[which(DrugRepWCC[,i]=="-1.1"),i]
  DrugRepWAC50filt[which(DrugRepWCC[,i]=="-1.2"),i] <- DrugRepWAC50[which(DrugRepWCC[,i]=="-1.2"),i]
  DrugRepWAC50filt[which(DrugRepWCC[,i]=="-2.1"),i] <- DrugRepWAC50[which(DrugRepWCC[,i]=="-2.1"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-1.1"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-1.1"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-1.2"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-1.2"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-2.1"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-2.1"),i]
}

### Add batch info
Batch <- rep(NA,length(DrugRepWAC50filt[,1]))
ProgID <- rep(NA,length(DrugRepWAC50filt[,1]))
for(i in 1:length(BatchInfo[,1])){
  Batch[which(DrugRepWAC50filt[,"Protocol.Name"]==BatchInfo[i,"PROTOCOL_NAME"])] <- BatchInfo[i,"Batch"]
  ProgID[which(DrugRepWAC50filt[,"Protocol.Name"]==BatchInfo[i,"PROTOCOL_NAME"])] <- BatchInfo[i,"ID"]
}

DrugRepWAC50filt <- cbind(Batch,ProgID,DrugRepWAC50)
DrugRepWAUCfilt <- cbind(Batch,ProgID,DrugRepWAUCfilt)

### Remove extra reps 3-8, they are NA for most drugs
DrugRepWAC503R <- DrugRepWAC50filt[which(DrugRepWAC50filt[,"Sample.Replicate"] %in% c(100,101,102)),]
DrugRepWAUC3R <- DrugRepWAUCfilt[which(DrugRepWAUCfilt[,"Sample.Replicate"] %in% c(100,101,102)),]
DrugRepWAC503RFC <- DrugRepWAC503R

na.count <- rep(NA,length(DrugRepWAC503R[1,]))
for(i in 6:length(na.count)){
  na.count[i] <- length(which(is.na(DrugRepWAC503R[,i])=="TRUE"))
}

### Scale AC50 matrix by Drug Mean
parentReps <- c("Malaria-NMAC-MxN-MKK2835_1","Malaria-NMAC-MxN-NHP1337_1","Malaria-NMAC-MxN-MKK2835_2","Malaria-NMAC-MxN-NHP1337_2","Malaria-NMAC-MxN-MKK2835_3","Malaria-NMAC-MxN-NHP1337_3","Malaria-NMAC-MxN-MKK2835_4","Malaria-NMAC-MxN-NHP1337_4")
ParentMean <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% parentReps),c(-1,-2,-3,-4)],na.rm=TRUE)
DrugMean <- colMeans(DrugRepWAC503R[,c(-1,-2,-3,-4)],na.rm=TRUE)

DrugRepWAC503RFC <- DrugRepWAC503R

### Get fold change in relation to Drug Mean
for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFC[notNA,i] <- DrugRepWAC503R[notNA,i]/DrugMean[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFC[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFC[,"Protocol.Name"],"_",DrugRepWAC503RFC[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFC[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==4)] <- "black"

library(RColorBrewer)
library(gplots)
library(vegan)
library(grDevices)

## Make vector of colors for values below threshold
rc1 <- colorpanel(80,hsv(0.175,1,0.95), "gray90")
## Make vector of colors for values above threshold
rc2 <- colorpanel(80,"gray90", hsv(0.62,1,0.9))
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51.
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256)
rb1 <- seq(-6, 0, length.out=80)
rb2 <- seq(0, 6, length.out=80)[-1]
rampbreaks <- c(-20,rb1,rb2,20)

png("Chemogenomic Heat Map logFC filt Drug Mean NS.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Get fold change in relation to Drug Mean
DrugRepWAC503RFCP <- DrugRepWAC503R
for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFCP[notNA,i] <- DrugRepWAC503R[notNA,i]/ParentMean[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFCP[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFCP[,"Protocol.Name"],"_",DrugRepWAC503RFCP[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFCP[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFCP[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFCP[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFCP[,"Batch"]==4)] <- "black"


png("Chemogenomic Heat Map logFC filt Parent Mean NS.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Get fold change in relation to Mean per Batch
DrugRepWAC503RFCB <- DrugRepWAC503R

DrugMean1 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Batch"] == 1),c(-1,-2,-3,-4)],na.rm=TRUE)
DrugMean2 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Batch"] == 2),c(-1,-2,-3,-4)],na.rm=TRUE)
DrugMean3 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Batch"] == 3),c(-1,-2,-3,-4)],na.rm=TRUE)
DrugMean4 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Batch"] == 4),c(-1,-2,-3,-4)],na.rm=TRUE)

### Plot each batch drug mean against total mean
plot(log2(DrugMean),log2(DrugMean1))
points(log2(DrugMean),log2(DrugMean2),col="red")
points(log2(DrugMean),log2(DrugMean3),col="blue")
points(log2(DrugMean),log2(DrugMean4),col="green")

Batch1 <- which(DrugRepWAC503R[,"Batch"]==1)
Batch2 <- which(DrugRepWAC503R[,"Batch"]==2)
Batch3 <- which(DrugRepWAC503R[,"Batch"]==3)
Batch4 <- which(DrugRepWAC503R[,"Batch"]==4)

for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFCB[intersect(notNA,Batch1),i] <- DrugRepWAC503R[intersect(notNA,Batch1),i]/DrugMean1[i-4]
  DrugRepWAC503RFCB[intersect(notNA,Batch2),i] <- DrugRepWAC503R[intersect(notNA,Batch2),i]/DrugMean2[i-4]
  DrugRepWAC503RFCB[intersect(notNA,Batch3),i] <- DrugRepWAC503R[intersect(notNA,Batch3),i]/DrugMean3[i-4]
  DrugRepWAC503RFCB[intersect(notNA,Batch4),i] <- DrugRepWAC503R[intersect(notNA,Batch4),i]/DrugMean4[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFCB[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFCB[,"Protocol.Name"],"_",DrugRepWAC503RFCB[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFCB[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFCB[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFCB[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFCB[,"Batch"]==4)] <- "black"

png("Chemogenomic Heat Map logFC filt Drug Mean By Batch NS.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Get fold change in relation to Parent Mean per Batch
DrugRepWAC503RFCBP <- DrugRepWAC503R
ParentMean1 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_1","Malaria-NMAC-MxN-NHP1337_1")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean2 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_2","Malaria-NMAC-MxN-NHP1337_2")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean3 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_3","Malaria-NMAC-MxN-NHP1337_3")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean4 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_4","Malaria-NMAC-MxN-NHP1337_4")),c(-1,-2,-3,-4)],na.rm=TRUE)

plot(log(ParentMean),log(ParentMean1),ylab="log2 each Parent Replicates",xlab="log2 Mean Parent Replicates")
points(log(ParentMean),log(ParentMean2),col="red")
points(log(ParentMean),log(ParentMean3),col="blue")
points(log(ParentMean),log(ParentMean4),col="green")

Batch1 <- which(DrugRepWAC503R[,"Batch"]==1)
Batch2 <- which(DrugRepWAC503R[,"Batch"]==2)
Batch3 <- which(DrugRepWAC503R[,"Batch"]==3)
Batch4 <- which(DrugRepWAC503R[,"Batch"]==4)

for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFCBP[intersect(notNA,Batch1),i] <- DrugRepWAC503R[intersect(notNA,Batch1),i]/ParentMean1[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch2),i] <- DrugRepWAC503R[intersect(notNA,Batch2),i]/ParentMean2[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch3),i] <- DrugRepWAC503R[intersect(notNA,Batch3),i]/ParentMean3[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch4),i] <- DrugRepWAC503R[intersect(notNA,Batch4),i]/ParentMean4[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFCBP[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFCBP[,"Protocol.Name"],"_",DrugRepWAC503RFCBP[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==4)] <- "black"

png("Chemogenomic Heat Map log FC filt Parent Mean By Batch.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Get fold change in relation to NHP1337 Mean per Batch
DrugRepWAC503RFCBP <- DrugRepWAC503R
ParentMean1 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-NHP1337_1")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean2 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-NHP1337_2")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean3 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-NHP1337_3")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean4 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-NHP1337_4")),c(-1,-2,-3,-4)],na.rm=TRUE)

plot(log(ParentMean),log(ParentMean1),ylab="log2 each Parent Replicates",xlab="log2 Mean Parent Replicates")
points(log(ParentMean),log(ParentMean2),col="red")
points(log(ParentMean),log(ParentMean3),col="blue")
points(log(ParentMean),log(ParentMean4),col="green")

Batch1 <- which(DrugRepWAC503R[,"Batch"]==1)
Batch2 <- which(DrugRepWAC503R[,"Batch"]==2)
Batch3 <- which(DrugRepWAC503R[,"Batch"]==3)
Batch4 <- which(DrugRepWAC503R[,"Batch"]==4)

for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFCBP[intersect(notNA,Batch1),i] <- DrugRepWAC503R[intersect(notNA,Batch1),i]/ParentMean1[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch2),i] <- DrugRepWAC503R[intersect(notNA,Batch2),i]/ParentMean2[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch3),i] <- DrugRepWAC503R[intersect(notNA,Batch3),i]/ParentMean3[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch4),i] <- DrugRepWAC503R[intersect(notNA,Batch4),i]/ParentMean4[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFCBP[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFCBP[,"Protocol.Name"],"_",DrugRepWAC503RFCBP[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==4)] <- "black"

png("Chemogenomic Heat Map log FC filt NHP1337 Mean By Batch.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Get fold change in relation to MKK2835 Mean per Batch
DrugRepWAC503RFCBP <- DrugRepWAC503R
ParentMean1 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_1")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean2 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_2")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean3 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_3")),c(-1,-2,-3,-4)],na.rm=TRUE)
ParentMean4 <- colMeans(DrugRepWAC503R[which(DrugRepWAC503R[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_4")),c(-1,-2,-3,-4)],na.rm=TRUE)

plot(log(ParentMean),log(ParentMean1),ylab="log2 each Parent Replicates",xlab="log2 Mean Parent Replicates")
points(log(ParentMean),log(ParentMean2),col="red")
points(log(ParentMean),log(ParentMean3),col="blue")
points(log(ParentMean),log(ParentMean4),col="green")

Batch1 <- which(DrugRepWAC503R[,"Batch"]==1)
Batch2 <- which(DrugRepWAC503R[,"Batch"]==2)
Batch3 <- which(DrugRepWAC503R[,"Batch"]==3)
Batch4 <- which(DrugRepWAC503R[,"Batch"]==4)

for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWAC503R[,i])=="FALSE")
  DrugRepWAC503RFCBP[intersect(notNA,Batch1),i] <- DrugRepWAC503R[intersect(notNA,Batch1),i]/ParentMean1[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch2),i] <- DrugRepWAC503R[intersect(notNA,Batch2),i]/ParentMean2[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch3),i] <- DrugRepWAC503R[intersect(notNA,Batch3),i]/ParentMean3[i-4]
  DrugRepWAC503RFCBP[intersect(notNA,Batch4),i] <- DrugRepWAC503R[intersect(notNA,Batch4),i]/ParentMean4[i-4]
}

### Make matrix for HeatMap and color code row labels by batch
mat=data.matrix((DrugRepWAC503RFCBP[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFCBP[,"Protocol.Name"],"_",DrugRepWAC503RFCBP[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFCBP[,"Batch"]==4)] <- "black"

png("Chemogenomic Heat Map log FC filt MKK2835 Mean By Batch.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Take average across reps for each drug
DrugWAC50FCB <- cbind(BatchInfo,matrix(NA,ncol=(length(DrugRepWAC503RFCB[1,])-4),nrow=length(BatchInfo[,1])))
NaNCount <- rep(NA, length(DrugWAC50FCB[1,]))
colnames(DrugWAC50FCB)[c(-1,-2,-3)] <- colnames(DrugRepWAC503RFCB)[c(-1,-2,-3,-4)]
for(i in 1:length(DrugWAC50FCB[,1])){
  DrugWAC50FCB[i,4:length(DrugWAC50FCB[1,])] <- colMeans(DrugRepWAC503RFCB[which(DrugRepWAC503RFCB[,"Protocol.Name"]==DrugWAC50FCB[i,"PROTOCOL_NAME"]),c(-1,-2,-3,-4)],na.rm=TRUE)
}
for(i in 4:length(DrugWAC50FCB[1,])){
  NaNCount[i] <- length(which(is.nan(DrugWAC50FCB[,i])=="TRUE"))
}

### Make matrix for HeatMap and color code row labels by batch

mat=data.matrix((DrugWAC50FCB[,-c(1,2,3,which(NaNCount>(length(DrugWAC50FCB[,1])*.1)))]))
rownames(mat)<-DrugWAC50FCB[,"PROTOCOL_NAME"]
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugWAC50FCB[,"Batch"]==1)] <- "blue"
rowCol[which(DrugWAC50FCB[,"Batch"]==2)] <- "red"
rowCol[which(DrugWAC50FCB[,"Batch"]==3)] <- "green"
rowCol[which(DrugWAC50FCB[,"Batch"]==4)] <- "black"


png("Chemogenomic Heat Map log FC filt Drug Mean By Batch Mean NS.png", height=24, width=15.5, units="in", res=220)
par(oma=c(4,0,0,15))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=1.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 6","4 to 5.9","2 to 3.9","No Change","-2 to -3.9","-4 to -5.9","<= -6","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()


### Take average across reps for each drug
DrugWAC50FC <- cbind(BatchInfo,matrix(NA,ncol=(length(DrugRepWAC503RFC[1,])-4),nrow=length(BatchInfo[,1])))
NaNCount <- rep(NA, length(DrugWAC50FCB[1,]))
colnames(DrugWAC50FCB)[c(-1,-2,-3)] <- colnames(DrugRepWAC503RFC)[c(-1,-2,-3,-4)]
for(i in 1:length(DrugWAC50FCB[,1])){
  DrugWAC50FC[i,4:length(DrugWAC50FCB[1,])] <- colMeans(DrugRepWAC503RFC[which(DrugRepWAC503RFC[,"Protocol.Name"]==DrugWAC50FC[i,"PROTOCOL_NAME"]),c(-1,-2,-3,-4)],na.rm=TRUE)
}
for(i in 4:length(DrugWAC50FC[1,])){
  NaNCount[i] <- length(which(is.nan(DrugWAC50FC[,i])=="TRUE"))
}

### Make matrix for HeatMap and color code row labels by batch

mat=data.matrix((DrugWAC50FC[,-c(1,2,3,which(NaNCount>(length(DrugWAC50FC[,1])*.1)))]))
rownames(mat)<-DrugWAC50FC[,"PROTOCOL_NAME"]
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugWAC50FC[,"Batch"]==1)] <- "blue"
rowCol[which(DrugWAC50FC[,"Batch"]==2)] <- "red"
rowCol[which(DrugWAC50FC[,"Batch"]==3)] <- "green"
rowCol[which(DrugWAC50FC[,"Batch"]==4)] <- "black"

### for each drug fill NAs with average
#mat[which(is.nan(mat)=="TRUE")] <- NA
#for(i in 1:length(mat[1,])){
# mat[which(is.na(mat[,i])=="TRUE"),i] <- mean(mat[,i],na.rm=TRUE)
#}

png("Chemogenomic Heat Map log mean FC filt to Mean no BC NS.png", height=24, width=15.5, units="in", res=220)
par(oma=c(4,0,0,15))

heatmap.2(log2(mat), trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=1.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 6","4 to 5.9","2 to 3.9","No Change","-2 to -3.9","-4 to -5.9","<= -6","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()



### ComBat batch corrected data
library("sva")
batch <- DrugRepWAC503RFC[,"Batch"]
mod <- model.matrix(~as.factor(ProgID),data=DrugRepWAC503RFC)

mat=data.matrix((DrugRepWAC503RFC[,-c(1,2,3,4,which(na.count>20))]))
rownames(mat)<-paste0(DrugRepWAC503RFC[,"Protocol.Name"],"_",DrugRepWAC503RFC[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503RFC[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503RFC[,"Batch"]==4)] <- "black"

combat_mat <- ComBat(dat=t(mat),batch=batch,mod=mod,par.prior=TRUE,prior.plots=FALSE)


### Take Mean Fold Change to MKK2835 Mean of batch corrected data
DrugRepWAC503RBC <- cbind(DrugRepWAC503R[,c(1,2,3,4,5)],t(combat_mat))

### Scale AC50 matrix by MKK2835 Mean 
MKK2835Mean <- colMeans(DrugRepWAC503RBC[which(DrugRepWAC503RBC[,"Protocol.Name"] %in% c("Malaria-NMAC-MxN-MKK2835_1","Malaria-NMAC-MxN-MKK2835_2","Malaria-NMAC-MxN-MKK2835_3","Malaria-NMAC-MxN-MKK2835_4")),c(-1,-2,-3,-4)],na.rm=TRUE)
DrugMean <- colMeans(DrugRepWAC503RBC[,c(-1,-2,-3,-4)],na.rm=TRUE)

DrugRepWAC503RBCFC <- DrugRepWAC503RBC
### Get fold change in relation to MKK2835 Mean
for(i in 5:length(DrugRepWAC503RBCFC[1,])){
  DrugRepWAC503RBCFC[,i] <- DrugRepWAC503RBC[,i]/MKK2835Mean[i-4]
}

mat=data.matrix((DrugRepWAC503RBC[,-c(1,2,3,4)]))
rownames(mat)<-paste0(DrugRepWAC503R[,"Protocol.Name"],"_",DrugRepWAC503R[,"Sample.Replicate"])
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugRepWAC503R[,"Batch"]==1)] <- "blue"
rowCol[which(DrugRepWAC503R[,"Batch"]==2)] <- "red"
rowCol[which(DrugRepWAC503R[,"Batch"]==3)] <- "green"
rowCol[which(DrugRepWAC503R[,"Batch"]==4)] <- "black"

library(RColorBrewer)
library(gplots)
library(vegan)
library(grDevices)

## Make vector of colors for values below threshold
rc1 <- colorpanel(80,hsv(0.175,1,0.95), "gray90")    
## Make vector of colors for values above threshold
rc2 <- colorpanel(80,"gray90", hsv(0.62,1,0.9))
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(1/8, 1, length.out=80)
rb2 <- seq(1, 8, length.out=80)[-1]
rampbreaks <- c(0,rb1,rb2,270843)



png("Chemogenomic Heat Map FC filt Drug Mean Combat.png", height=24, width=15.5, units="in", res=220)
par(oma=c(6,0,0,10))

heatmap.2(mat, trace="none", density="none", 
          scale="none",col = rampcols, breaks = rampbreaks, 
          cexRow=0.5,cexCol=0.5,na.rm=T,na.color="grey",key=FALSE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
legend(0,1,c("Sensitive",">= 8","4 to 7.9","2 to 3.9","No Change","-2 to -3.9","-4 to -7.9","<= -8","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()


write.csv(DrugRepWAC503RBC,"MKK2835xNHP1337_DrugMeanFoldChangetoBatchMean.csv")
