library(RColorBrewer)
library(gplots)

library(tidyverse)
setwd('/Users/macke/Documents/RSA/MKK2835xNHP1337/normalized')

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

### Transform to IC50 and AUC matrices with rows as individual reps and columns as drugs
library(reshape2)


#DrugRepWIC50mean <- dcast(Individual_Drug_Replicates, Protocol.Name + Assay.Plate.ID ~ Sample.ID, fun.aggregate=mean, value.var="IC50..uM.")
DrugRepWIC50 <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="IC50..uM.")
DrugRepWAUC <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="AUC")
DrugRepWCC <- dcast(Individual_Drug_Replicates, Protocol.Name + Sample.Replicate ~ Sample.ID, value.var="CC.v2")

### Screen any data from curves that are not CURVE_CLASS -1.1, -1.2 or -2.1
DrugRepWIC50filt <- DrugRepWIC50
DrugRepWAUCfilt <- DrugRepWAUC
DrugRepWIC50filt[,c(-1,-2)] <- NA
DrugRepWAUCfilt[,c(-1,-2)] <- NA
for(i in 3:length(DrugRepWIC50filt)){
  DrugRepWIC50filt[which(DrugRepWCC[,i]=="-1.1"),i] <- DrugRepWIC50[which(DrugRepWCC[,i]=="-1.1"),i]
  DrugRepWIC50filt[which(DrugRepWCC[,i]=="-1.2"),i] <- DrugRepWIC50[which(DrugRepWCC[,i]=="-1.2"),i]
  DrugRepWIC50filt[which(DrugRepWCC[,i]=="-2.1"),i] <- DrugRepWIC50[which(DrugRepWCC[,i]=="-2.1"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-1.1"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-1.1"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-1.2"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-1.2"),i]
  DrugRepWAUCfilt[which(DrugRepWCC[,i]=="-2.1"),i] <- DrugRepWAUC[which(DrugRepWCC[,i]=="-2.1"),i]
}

### Add batch info
Batch <- rep(NA,length(DrugRepWIC50filt[,1]))
ProgID <- rep(NA,length(DrugRepWIC50filt[,1]))
for(i in 1:length(BatchInfo[,1])){
  Batch[which(DrugRepWIC50filt[,"Protocol.Name"]==BatchInfo[i,"PROTOCOL_NAME"])] <- BatchInfo[i,"Batch"]
  ProgID[which(DrugRepWIC50filt[,"Protocol.Name"]==BatchInfo[i,"PROTOCOL_NAME"])] <- BatchInfo[i,"ï..ID"]
}

DrugRepWIC50filt <- cbind(Batch,ProgID,DrugRepWIC50)
DrugRepWAUCfilt <- cbind(Batch,ProgID,DrugRepWAUCfilt)


na.count <- rep(NA,length(DrugRepWIC50[1,]))
for(i in 6:length(na.count)){
  na.count[i] <- length(which(is.na(DrugRepWIC50[,i])=="TRUE"))
}

### Make replicate matrix to fill in residuals
DrugRepWIC50RC <- DrugRepWIC50filt

### Get residuals from LM including intercept and drug mean and fill in drug matrix
for(i in 5:length(na.count)){
  notNA <- which(is.na(DrugRepWIC50filt[,i])=="FALSE")
  if(length(notNA)>length(DrugRepWIC50filt[,i])*0.1){
    LMDrug <- lm(DrugRepWIC50filt[,i]~as.factor(DrugRepWIC50filt[,1]))
    DrugRepWIC50RC[names(residuals(LMDrug)),i] <- residuals(LMDrug)
  } else {
    DrugRepWIC50RC[,i] <- NA
  }
}

na.count <- rep(NA,length(DrugRepWIC50RC[1,]))
for(i in 5:length(na.count)){
  na.count[i] <- length(which(is.na(DrugRepWIC50RC[,i])=="TRUE"))
}


### Take average across reps for each drug 
DrugWIC50 <- cbind(BatchInfo,matrix(NA,ncol=(length(DrugRepWIC50[1,])-2),nrow=length(BatchInfo[,1])))
NaNCount <- rep(NA, length(DrugWIC50[1,]))
colnames(DrugWIC50)[c(-1,-2,-3)] <- colnames(DrugRepWIC50)[c(-1,-2)]
for(i in 1:length(DrugWIC50[,1])){
  DrugWIC50[i,c(-1,-2,-3)] <- colMeans(DrugRepWIC50[which(DrugRepWIC50[,"Protocol.Name"]==DrugWIC50[i,"PROTOCOL_NAME"]),c(-1,-2)],na.rm=TRUE)
}
for(i in 4:length(DrugWIC50[1,])){
  NaNCount[i] <- length(which(is.nan(DrugWIC50[,i])=="TRUE"))
}



### Make matrix for HeatMap and color code row labels by batch

mat=data.matrix((DrugWIC50[,-c(1,2,3,which(NaNCount>(length(DrugWIC50[,1])*.1)))]))
rownames(mat)<-DrugWIC50[,"PROTOCOL_NAME"]
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugWIC50[,"Batch"]==1)] <- "blue"
rowCol[which(DrugWIC50[,"Batch"]==2)] <- "red"
rowCol[which(DrugWIC50[,"Batch"]==3)] <- "green"
rowCol[which(DrugWIC50[,"Batch"]==4)] <- "black"
colfunc <- colorRampPalette(c("blue", "white", "yellow"))

png("Chemogenomic Heat Map Mean IC50 NS check.png", height=24, width=15.5, units="in", res=220)
par(oma=c(4,0,0,15))

heatmap.2(mat, trace="none", density="none", 
          scale="column",col = colfunc(40), 
          cexRow=1.5,cexCol=0.5,na.rm=T,na.color="grey",key=TRUE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)
#legend(0,1,c("Sensitive",">= 6","4 to 5.9","2 to 3.9","No Change","-2 to -3.9","-4 to -5.9","<= -6","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

### Take average across reps for each drug for IC50 corrected by batch
DrugWIC50B <- cbind(BatchInfo,matrix(NA,ncol=(length(DrugRepWIC50RC[1,])-4),nrow=length(BatchInfo[,1])))
NaNCount <- rep(NA, length(DrugWIC50B[1,]))
colnames(DrugWIC50B)[c(-1,-2,-3)] <- colnames(DrugRepWIC50RC)[c(-1,-2,-3,-4)]
for(i in 1:length(DrugWIC50B[,1])){
  DrugWIC50B[i,c(-1,-2,-3)] <- colMeans(DrugRepWIC50RC[which(DrugRepWIC50RC[,"Protocol.Name"]==DrugWIC50B[i,"PROTOCOL_NAME"]),c(-1,-2,-3,-4)],na.rm=TRUE)
}
for(i in 4:length(DrugWIC50B[1,])){
  NaNCount[i] <- length(which(is.nan(DrugWIC50B[,i])=="TRUE"))
}

good = as.vector(as.factor(names_gt$name))
DrugWIC50B= drugB[ , -which(names(drugB) %in% good)]

DrugWIC50B = subset(drugB,select = -good)
### Make matrix for HeatMap and color code row labels by batch

mat=data.matrix((DrugWIC50B[,-c(1,2,3,which(NaNCount>(length(DrugWIC50B[,1])*.1)))]))
rownames(mat)<-DrugWIC50B[,"PROTOCOL_NAME"]
rowCol <- rep(NA,length(rownames(mat)))
rowCol[which(DrugWIC50B[,"Batch"]==1)] <- "blue"
rowCol[which(DrugWIC50B[,"Batch"]==2)] <- "red"
rowCol[which(DrugWIC50B[,"Batch"]==3)] <- "green"
rowCol[which(DrugWIC50B[,"Batch"]==4)] <- "black"
colfunc <- colorRampPalette(c("blue", "white", "yellow"))

png("Chemogenomic Heat Map Batch Mean NS check.png", height=24, width=15.5, units="in", res=220)
par(oma=c(4,0,0,15))

heatmap.2(mat, trace="none", density="none", 
          scale="column",col = colfunc(40), 
          cexRow=1.5,cexCol=0.5,na.rm=T,na.color="grey",key=TRUE,
          hclust=function(x) hclust(x,method="average"),
          distfun=function(x) as.dist((1-cor(t(x),method="spearman",use="pairwise.complete.obs"))/2),colRow=rowCol)

#legend(0,1,c("Sensitive",">= 6","4 to 5.9","2 to 3.9","No Change","-2 to -3.9","-4 to -5.9","<= -6","Resistant"),col=c(hsv(0.175,1,0.95),hsv(0.175,1,0.95),hsv(0.175,0.6,0.95),hsv(0.175,0.3,0.95),"gray90",hsv(0.62,0.3,0.9),hsv(0.62,0.6,0.9),hsv(0.62,1,0.9),hsv(0.62,1,0.9)),pch=c(NA,15,15,15,15,15,15,15,NA),bty="n",y.intersp=1,cex=1.5)

dev.off()

write.csv(DrugWIC50,"Mean Drug IC50.csv")
write.csv(DrugWIC50B,"Mean Drug IC50 Batch.csv")
write.csv(DrugNames,"DrugNames.csv")
write.csv(names,"levenes_test.csv")


names_gt  = read.csv("levenes_test_greaterthan05.csv", 
                  header = T, 
                  #row.names = 2, 
                  stringsAsFactors = T, 
                  fileEncoding="UTF-8-BOM")



###############################################################################
install.packages("inferr")
library(inferr)
library(car)
data(PlantGrowth)
drug <- read.csv("Mean Drug IC50.csv",header=TRUE
                 ,row.names = 2,
                 stringsAsFactors = T, 
                 fileEncoding="UTF-8-BOM")


drugB  = read.csv("Mean Drug IC50 Batch.csv", 
                    header = TRUE, 
                    #row.names = 2, 
                    stringsAsFactors = T, 
                    fileEncoding="UTF-8-BOM")

names  = read.csv("drug_names.csv", 
                  header = T, 
                  #row.names = 2, 
                  stringsAsFactors = F, 
                  fileEncoding="UTF-8-BOM")

for (i in 1:nrow(names)){
  dname = names[i,1]
  print(paste("res = leveneTest(",dname," ~ Batch, data = drugB)", sep=""))
  print(paste("names[i,2] = res[1,3]"))
  print(paste("i=i+1"))
}

p_value = res2$p_lev
p_value2 = res[1,3]
names[i,2] = p_value

res2 = infer_levene_test(drug, NCGC00015251_13, group_var = Batch)

infer_levene_test(drugB, drugB$dname, group_var = Batch)

vecx = colnames(drugB[3:length(drugB)])

dt.out = drugB[,lapply(.SD, leveneTest, group = Batch, .SDcols=vecx)]
as.factor(dname)

i=1
res = leveneTest(NCGC00015251_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00015256_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00015323_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00015334_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00015580_22 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00015819_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016256_24 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016373_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016421_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016612_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016789_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016913_14 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00016961_16 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00017363_31 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00017376_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00018301_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00022019_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00022678_37 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00024246_30 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00024415_46 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00024995_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00025060_36 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00025125_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00025155_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00025341_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00052276_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00090168_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00090208_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00091112_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00091231_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00091944_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00093976_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00094926_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00095055_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00095196_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00098104_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00098192_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00100513_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00100738_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00101991_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00108995_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00114972_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00115685_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00118982_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00159483_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00159519_14 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00159544_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00159574_09 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00160163_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00160391_09 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00160546_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00161621_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00161634_20 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00161679_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00161831_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00162341_09 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00162398_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00162400_09 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163451_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163469_16 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163548_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163699_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163700_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163818_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00163847_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00164493_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00164559_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00164591_10 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00164631_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00165865_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00166121_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167328_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167481_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167490_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167507_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167518_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00167558_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00168085_26 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00173359_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00179034_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00179250_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00179500_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00180461_10 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00181109_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00181319_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00181375_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00181776_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00181913_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00182996_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183024_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183121_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183285_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183656_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183844_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00183878_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00185000_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00185323_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00189073_17 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00189220_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00229704_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00238453_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00238454_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00242051_10 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00242478_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00242506_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00244252_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00244256_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00246187_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00247878_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00247954_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00249372_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00249613_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00249896_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00249897_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00250381_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00250388_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00250400_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00250405_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00250408_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00253607_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00253624_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00253626_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00253770_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00256767_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00262526_14 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00262608_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263087_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263117_06 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263152_15 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263156_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263177_14 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263180_18 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263183_14 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00263785_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00271789_18 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00319019_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00344515_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00344584_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00344630_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00345539_13 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346444_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346459_10 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346475_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346484_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346486_07 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346487_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346493_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346505_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346542_10 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346590_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346626_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346654_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346659_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346673_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346681_12 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346692_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346716_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346717_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346718_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346808_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346877_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00346940_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00347048_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00347284_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00351481_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00351607_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00356154_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00371897_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00371898_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00378460_09 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00378599_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00378614_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00378895_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00378978_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00379241_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00380396_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386263_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386281_08 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386310_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386313_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386367_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386398_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386425_11 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00386616_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00387037_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00387811_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00387877_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00388587_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389104_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389265_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389323_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389337_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389594_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389695_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00389773_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00390442_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00390643_03 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00390699_04 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00412395_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00412536_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00412552_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00413871_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00417201_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00421178_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00441011_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00454902_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00465444_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00479176_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00479600_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00481606_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00482483_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00483024_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00483924_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00484987_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00485899_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00487207_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00509859_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00589068_02 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00589069_05 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00655580_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00655581_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1
res = leveneTest(NCGC00686000_01 ~ Batch, data = drugB)
names[i,2] = res[1,3]
i=i+1






