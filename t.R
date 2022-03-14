library(ggplot2)
library(ggvenn)
library(ggVennDiagram)
setwd("C:/Users/megan/Ferdig-Lab/Fixed")
T4C <- read.csv("T4FeatureCount_150.csv", header = TRUE, stringsAsFactors = FALSE)
T4N <- read.csv("T4FeatureCount_ELO.csv", header = TRUE, stringsAsFactors = FALSE)
common <- intersect(T4C[,3], T4N[,3])
t <- as.data.frame(common)
write.csv(t, "commonT4_ELO_RSA.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
only4C <- c()
t <- T4C[,3]%in%common
for (i in 1:nrow(T4C)){
  if(t[i] == FALSE){
    only4C <- c(only4C, T4C[i,3])
  }
}
only4C
only4N <- c()
tN <- T4N[,3]%in%common
for (i in 1:nrow(T4N)){
  if(tN[i] == FALSE){
    only4N <- c(only4N, T4N[i,3])
  }
}
A <- list(T4C[,3],T4N[,3])
ggVennDiagram(A, category.names = c("T4 RSA", "T4 ELO"),color = 1, lwd = 0.7)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
T30C <- read.csv("T30FeatureCount_150.csv", header = TRUE, stringsAsFactors = FALSE)
T30N <- read.csv("T30FeatureCount_ELO.csv", header = TRUE, stringsAsFactors = FALSE)
common30 <- intersect(T30C[,3], T30N[,3])
t <- as.data.frame(common30)
write.csv(t, "commonT30_RSA_ELO.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
only30C <- c()
t <- T30C[,3]%in%common30
for (i in 1:nrow(T30C)){
  if(t[i] == FALSE){
    only30C <- c(only30C, T30C[i,3])
  }
}
only30C
only30N <- c()
tN <- T30N[,3]%in%common30
for (i in 1:nrow(T30N)){
  if(tN[i] == FALSE){
    only30N <- c(only30N, T30N[i,3])
  }
}
B <- list(T30C[,3],T30N[,3])
ggVennDiagram(B, category.names = c("T30 RSA", "T30 ELO"),color = 1, lwd = 0.7)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")


T44C <- read.csv("T44FeatureCount_150.csv", header = TRUE, stringsAsFactors = FALSE)
T44N <- read.csv("T44FeatureCount_ELO.csv", header = TRUE, stringsAsFactors = FALSE)
common44 <- intersect(T44C[,3], T44N[,3])
t <- as.data.frame(common44)
write.csv(t, "commonT44_RSA_ELO.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)
only44C <- c()
t <- T44C[,3]%in%common44
for (i in 1:nrow(T44C)){
  if(t[i] == FALSE){
    only44C <- c(only44C, T44C[i,3])
  }
}
only44C
only44N <- c()
tN <- T44N[,3]%in%common44
for (i in 1:nrow(T44N)){
  if(tN[i] == FALSE){
    only44N <- c(only44N, T44N[i,3])
  }
}
only44N

C <- list(T44C[,3],T44N[,3])
ggVennDiagram(C, category.names = c("T44 RSA", "T44 ELO"),color = 1, lwd = 0.7)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

D <- list(T4C[,3],T30C[,3],T44C[,3])
ggVennDiagram(D, category.names = c("T4 150", "T30", "T44"),color = 1, lwd = 0.7)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
common4.30 <- intersect(T4C[,3],T30C[,3])
common4.44 <- intersect(T4C[,3], T44C[,3])
common30.44 <- intersect(T30C[,3], T44C[,3])
E <- list(T4N[,3],T30N[,3],T44N[,3])
ggVennDiagram(E, category.names = c("T4 ELO", "T30", "T44"),color = 1, lwd = 0.7)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
common4.30N <- intersect(T4N[,3],T30N[,3])
common4.44N <- intersect(T4N[,3], T44N[,3])
common30.44N <- intersect(T30N[,3], T44N[,3])
as.data.frame(common30.44N)
