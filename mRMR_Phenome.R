setwd("C:/Users/megan/Downloads")
library(mRMRe)
library(tidyverse)
correctedt4 <- read.table("corrected_norm_noreps_T4_PlateCRCSV_.txt", header = FALSE, stringsAsFactors = FALSE)
corrected30 <- read.table("corrected_norm_noreps_T30_PlateCRCSV_.txt", header = TRUE, stringsAsFactors = FALSE)
corrected44 <- read.table("corrected_norm_noreps_T44_PlateCRCSV_.txt", header = TRUE, stringsAsFactors = FALSE)
x <- read.table("corrected_norm_noreps_T4_PlateCRC_.txt", header = FALSE, stringsAsFactors = FALSE)
correct30 <- read.table("corrected_norm_noreps_T30_PlateCRC_.txt", header = FALSE, stringsAsFactors = FALSE)
correct44 <- read.table("corrected_norm_noreps_T44_PlateCRC_.txt", header = FALSE, stringsAsFactors = FALSE)

mRMR <- function(x){
rsa <- read.csv("NF54xNHP4026.csv", stringsAsFactors = FALSE)
rsa <- subset(rsa, select = c(X.1, RSA_120_PHENOME))
rsaWithParents <- rsa[grep("N", rsa$X.1),]
rsa <- rsa[grep("AC", rsa$X.1),]
rsa <- rbind(rsaWithParents, rsa)
header <- c("start", "end", "id", "gid", "strd", "GF_PL01_AC025_4hpi_708", "GF_PL02_AC004_4hpi_711", "GF_PL02_AC006_4hpi_711", "GF_PL02_AC008_4hpi_708", "GF_PL02_AC027_4hpi_517", 
            "GF_PL02_AC028_4hpi_517", "GF_PL02_AC033_4hpi_708", "GF_PL02_AC034_4hpi_517", "GF_PL02_AC049_4hpi_517", "GF_PL02_AC050_4hpi_517", "GF_PL02_AC088_4hpi_1016", 
            "GF_PL02_AC094_4hpi_711", "GF_PL02_AC096_4hpi_711", "GF_PL02_AC098_4hpi_1016", "GF_PL02_AC100_4hpi_1016", "GF_PL02_AC110_4hpi_1016", "GF_PL02_AC120_4hpi_517", "GF_PL02_AC124_4hpi_1016", 
            "GF_PL03_AC005_4hpi_708", "GF_PL03_AC023_4hpi_228", "GF_PL03_AC024_4hpi_228", "GF_PL03_AC041_4hpi_604", "GF_PL03_AC042_4hpi_604", "GF_PL03_AC043_4hpi_702", "GF_PL03_AC044_4hpi_604", 
            "GF_PL03_AC045_4hpi_604", "GF_PL03_AC046_4hpi_228", "GF_PL03_AC047_4hpi_702", "GF_PL03_AC048_4hpi_610", "GF_PL03_AC056_4hpi_228"	, "GF_PL03_AC057_4hpi_702", "GF_PL03_AC070_4hpi_610" ,
            "GF_PL03_AC077_4hpi_228" , "GF_PL03_AC079_4hpi_228" , 	"GF_PL03_AC093_4hpi_604" ,	"GF_PL03_AC097_4hpi_604" ,	"GF_PL03_AC103_4hpi_702"	, "GF_PL03_AC109_4hpi_604" , 
            "GF_PL03_AC118_4hpi_604" , 	"GF_PL03_AC127_4hpi_702" ,	"GF_PL03_A_NHP4026_4hpi_708"	, "GF_PL03_C_NF54gfp_4hpi_708" , 	"GF_PL05_AC125_4hpi_517" ,	"GF_PL05a_AC075_4hpi_1216" ,	"GF_PL05a_AC082_4hpi_1216" ,
            "GF_PL05a_AC130_4hpi_1216" ,	"GF_PL05b_AC030_4hpi_1216" ,	"GF_PL05b_AC032_4hpi_1216" , "GF_PL05b_AC038_4hpi_1216"	, "GF_PL05b_AC074_4hpi_1216")
x <- data.frame(x)
fixed <- data.frame(t(x[-1]))
fixed <- cbind(fixed, row.names(fixed))
colnames(fixed) <- fixed[4,]
data <- fixed[-c(1:5),]
for(i in 1:nrow(rsa)){
  j <- grep(rsa[i,1], data[,5156])
  data[j,5157]<- rsa[i,2]
}
datas <- subset(data, select = -c(gid))
dataset <- data.frame()
for (i in 1:nrow(datas)){
  for (j in 1:ncol(datas)){
    dataset[i,j] <- as.numeric(datas[i,j])
  }
}
colnames(dataset) <- colnames(datas)
rownames(dataset) <- rownames(datas)
mpackage <- mRMR.data(data = dataset)
mdpackage <- mRMR.ensemble(data = mpackage, target_indices = c(1), solution_count = 1, feature_count = 100)

compare <- data.frame(mdpackage@filters[["1"]], mdpackage@scores[["1"]])
for (i in 1:nrow(compare)){
  compare[i,3] <- x[compare[i,1]+1,5]
  
}
write.csv(compare, "T44FeatureCountsNonCV.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
mRMR(correct44)
write.csv(compare, "T4FeatureCount.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
#Add in data for RSA Phenome
#Grep and match RSA to the samples
#make geneid as colnames 
#make mRMR ensemble 
t4 <- read.csv("T4FeatureCount.csv", header = TRUE, stringsAsFactors = FALSE)
t4PlateCRC <- read.csv("T4FeatureCountPlateCRCFixed.csv", header = TRUE, stringsAsFactors = FALSE)
t44 <- read.csv("T44FeatureCountsCorrected.csv", header = TRUE, stringsAsFactors = FALSE)
t44PlateCRC <- read.csv("T44FeatureCountsNonCV.csv", header = TRUE, stringsAsFactors = FALSE)
common <- c()
for (i in 1:nrow(t4)){
  for (j in 1:nrow(t4PlateCRC)){
  if (t4[i,1] == t4PlateCRC[j,1]){
    common <- c(common, t4[i,3])
  }
  }  
}
sameLocation <- c()
for (i in 1:nrow(t4)){
  if (t4[i,1] == t4PlateCRC[i,1]){
    sameLocation <- c(sameLocation, t4[i,3])
  }
}
t4[31,3]
write.csv(common, "commonGenest4.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.csv(sameLocation, "sameLoactiont4.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
'''
for(i in 1:ncol(correctedt4)){
  for (j in 1:length(header)){
  correctedt4 <- rename(correctedt4, j = V)
  }
}
'''
t4 <- data.frame(t4)
t4fix <- data.frame(t(t4[-1]))
colnames(t4fix) <- t4[,1]

t4dd <- mRMR.data(data = t4fix)
t <- mRMR.classic(data = t4dd, target_indices = c(1), feature_count = 100)

network <- new("mRMRe.Network", data = demo, target_indices = c(1, 2),
               levels = c(2, 1), layers = 1)
set.thread.count(2)
data(cgps)
data.annot <- data.frame(cgps.annot)
data <- data.frame(cgps.ge)
data.cgps <- data.frame(cgps.ic50, cgps.ge)

demo <- mRMR.data(data = data.cgps)
demod <- mRMR.ensemble(data = demo, target_indices = c(1), solution_count = 1, feature_count = 50)
?solutions
scores(demo, solutions)
#solutions(object, mi_threshold, causality_threshold, with_fixed_features)
solution <- solutions(t)
target(t)
visualize(network)
 e <- solutions(network)
lc <- c(1:1001)
lt <- data.frame(t(sapply(l, lc)))
library(dplyr)

