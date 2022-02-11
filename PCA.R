library(tidyverse)
summary(USArrests)
head("USArrests", 10)
setwd("/Users/macke/Documents/PCA")

# Import data
ma.pollen.raw <- read.csv("NHP1337_MKK2835_IC50_ONLY.csv", header=TRUE, row.names=1, sep=",", check.names=FALSE)
head(ma.pollen.raw)
# Tidy/subset imported data
# Remove the first four columns
ma.pollen <- ma.pollen.raw[-1:-2]
# Remove total AP/NAP and pollen concentration columns
#ma.pollen <- ma.pollen[-49:-51]
# Remove columns where total abundanace is less than 10%
ma.sum <- colSums(ma.pollen) # Calculate column sums
#ma.pollen1 <- ma.pollen[, which(ma.sum > 10)] # Subset data

#################

# PCA using base function - prcomp()
ma.pollen.ot = data.frame(t(na.omit(t(ma.pollen))))
ma.pollen.ot.log = log(ma.pollen.ot)
p <- prcomp(ma.pollen.ot, scale=TRUE)
# Summary
s <- summary(p)

layout(matrix(1:2, ncol=2))
screeplot(p)
screeplot(p, type="lines")

layout(matrix(1:1, ncol=2))
biplot(p)

# Create groups
pch.group <- c(rep(21, times=2), rep(22, times=5), rep(24, times=11), rep(22, times=5), rep(24, times=11), rep(22, times=4), rep(24, times=10), rep(22, times=4), rep(24, times=10))
col.group <- c(rep("black", times=2), rep("skyblue2", times=16), rep("darkorchid1", times=16), rep("gold", times=14), rep("green2", times=14))

# Plot individuals
plot(p$x[,1], 
     p$x[,2], 
     xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), 
     ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), 
     pch=pch.group, 
     col="black", 
     bg=col.group, cex=2, las=1, asp=1)

abline(v=0, lty=2, col="grey50")
abline(h=0, lty=2, col="grey50")
# Add labels
text(p$x[,1], p$x[,2], 
     labels=row.names(p$x), 
     pos=c(1,3,4,2), font=2)

# Get co-ordinates of variables (loadings), and multiply by 10
l.x <- p$rotation[,1]*50
l.y <- p$rotation[,2]*50
# Draw arrows
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)

# Label position
l.pos <- l.y # Create a vector of y axis coordinates
lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
hi <- which(l.y > 0) # Get variables on the top half
# Replace values in the vector
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")

# Variable labels
text(l.x, l.y, labels=row.names(p$rotation), col="red", pos=l.pos)
# Add legend
legend("topleft", legend=c("0", "1", "2","3","4"), col="black", pt.bg=c("red"," skyblue2", "black", "gold", "green2"), pch=25, pt.cex=1.5)


# Get individuals (observations) as a matrix
tab <- matrix(c(p$x[,1], p$x[,2]), ncol=2)
# Calculate correlations
c0 <- cor(tab[1:2,])
c1 <- cor(tab[3:18,])
c2 <- cor(tab[19:34,])
c3 <- cor(tab[35:48,])
c4 <- cor(tab[49:62,])


library(ellipse)
# Plot ellipse
polygon(ellipse(c0*(max(abs(p$rotation))*10), centre=colMeans(tab[1:2,]), level=0.95), col=adjustcolor("red", alpha.f=0.25), border="red")
polygon(ellipse(c1*(max(abs(p$rotation))*10), centre=colMeans(tab[3:18,]), level=0.95), col=adjustcolor("skyblue2", alpha.f=0.25), border="skyblue")
polygon(ellipse(c2*(max(abs(p$rotation))*10), centre=colMeans(tab[19:34,]), level=0.95), col=adjustcolor("black", alpha.f=0.25), border="black")
polygon(ellipse(c3*(max(abs(p$rotation))*10), centre=colMeans(tab[35:48,]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold")
polygon(ellipse(c4*(max(abs(p$rotation))*10), centre=colMeans(tab[49:62,]), level=0.95), col=adjustcolor("green2", alpha.f=0.25), border="green2")

# Plot biplot
par(mar=c(4.5, 4.5, 1, 1))
plot(p$x[,1], p$x[,2], xlim=c(-13, 18), ylim=c(-14, 11), xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), 
     pch=pch.group, col="black", bg=col.group, cex=2.5, cex.axis=1.5, cex.lab=1.5, las=1,
     panel.first= {
       polygon(ellipse(c0*(max(abs(p$rotation))*10), centre=colMeans(tab[1:2,]), level=0.95), col=adjustcolor("black", alpha.f=0.25), border="red")
       polygon(ellipse(c1*(max(abs(p$rotation))*10), centre=colMeans(tab[3:18,]), level=0.95), col=adjustcolor("skyblue2", alpha.f=0.25), border="skyblue")
       polygon(ellipse(c2*(max(abs(p$rotation))*10), centre=colMeans(tab[19:34,]), level=0.95), col=adjustcolor("darkorchid1", alpha.f=0.25), border="darkorchid1")
       polygon(ellipse(c3*(max(abs(p$rotation))*10), centre=colMeans(tab[35:48,]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold")
       polygon(ellipse(c4*(max(abs(p$rotation))*10), centre=colMeans(tab[49:62,]), level=0.95), col=adjustcolor("green2", alpha.f=0.25), border="green2")       
       abline(v=0, lty=2, col="grey50") 
       abline(h=0, lty=2, col="grey50")
     },
     panel.last=arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red3", length=0.1, lwd=1.5))
# Labels
#text(p$x[,1], p$x[,2], labels=row.names(p$x), pos=c(1,3,4,2), font=2)
#text(l.x, l.y, labels=row.names(p$rotation), col="red", pos=l.pos)
# Add legend
#legend("topleft", legend=c("Tislit", "Sidi Ali", "Michliffen"), title="Sample Area", inset=c(0.01, 0.01), col="black", 
#       pt.bg=c("skyblue2", "gold", "green2"), pch=c(21, 22, 24), pt.cex=2, cex=1.5, bty="n")
legend("topleft", legend=c("0", "1", "2","3","4"), col="black", pt.bg=c("black"," skyblue2", "darkorchid1", "gold", "green2"), pch=25, pt.cex=1.5)


# Plot biplot
par(mar=c(4.5, 4.5, 1, 1))
plot(p$x[,1], p$x[,3], xlim=c(-13, 18), ylim=c(-14, 11), xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 3 (", round(s$importance[8]*100, 1), "%)", sep = ""), 
     pch=pch.group, col="black", bg=col.group, cex=2.5, cex.axis=1.5, cex.lab=1.5, las=1,
     panel.first= {
       polygon(ellipse(c0*(max(abs(p$rotation))*10), centre=colMeans(tab[1:2,]), level=0.95), col=adjustcolor("black", alpha.f=0.25), border="red")
       polygon(ellipse(c1*(max(abs(p$rotation))*10), centre=colMeans(tab[3:18,]), level=0.95), col=adjustcolor("skyblue2", alpha.f=0.25), border="skyblue")
       polygon(ellipse(c2*(max(abs(p$rotation))*10), centre=colMeans(tab[19:34,]), level=0.95), col=adjustcolor("darkorchid1", alpha.f=0.25), border="darkorchid1")
       polygon(ellipse(c3*(max(abs(p$rotation))*10), centre=colMeans(tab[35:48,]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold")
       polygon(ellipse(c4*(max(abs(p$rotation))*10), centre=colMeans(tab[49:62,]), level=0.95), col=adjustcolor("green2", alpha.f=0.25), border="green2")       
       abline(v=0, lty=2, col="grey50") 
       abline(h=0, lty=2, col="grey50")
     },
     panel.last=arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red3", length=0.1, lwd=1.5))
# Labels
#text(p$x[,1], p$x[,2], labels=row.names(p$x), pos=c(1,3,4,2), font=2)
#text(l.x, l.y, labels=row.names(p$rotation), col="red", pos=l.pos)
# Add legend
#legend("topleft", legend=c("Tislit", "Sidi Ali", "Michliffen"), title="Sample Area", inset=c(0.01, 0.01), col="black", 
#       pt.bg=c("skyblue2", "gold", "green2"), pch=c(21, 22, 24), pt.cex=2, cex=1.5, bty="n")
legend("topleft", legend=c("0", "1", "2","3","4"), col="black", pt.bg=c("black"," skyblue2", "darkorchid1", "gold", "green2"), pch=25, pt.cex=1.5)


# Plot biplot
par(mar=c(4.5, 4.5, 1, 1))
plot(p$x[,2], p$x[,3], xlim=c(-13, 18), ylim=c(-14, 11), xlab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), ylab=paste("PCA 3 (", round(s$importance[8]*100, 1), "%)", sep = ""), 
     pch=pch.group, col="black", bg=col.group, cex=2.5, cex.axis=1.5, cex.lab=1.5, las=1,
     panel.first= {
       polygon(ellipse(c0*(max(abs(p$rotation))*10), centre=colMeans(tab[1:2,]), level=0.95), col=adjustcolor("black", alpha.f=0.25), border="red")
       polygon(ellipse(c1*(max(abs(p$rotation))*10), centre=colMeans(tab[3:18,]), level=0.95), col=adjustcolor("skyblue2", alpha.f=0.25), border="skyblue")
       polygon(ellipse(c2*(max(abs(p$rotation))*10), centre=colMeans(tab[19:34,]), level=0.95), col=adjustcolor("darkorchid1", alpha.f=0.25), border="darkorchid1")
       polygon(ellipse(c3*(max(abs(p$rotation))*10), centre=colMeans(tab[35:48,]), level=0.95), col=adjustcolor("gold", alpha.f=0.25), border="gold")
       polygon(ellipse(c4*(max(abs(p$rotation))*10), centre=colMeans(tab[49:62,]), level=0.95), col=adjustcolor("green2", alpha.f=0.25), border="green2")       
       abline(v=0, lty=2, col="grey50") 
       abline(h=0, lty=2, col="grey50")
     },
     panel.last=arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red3", length=0.1, lwd=1.5))
# Labels
#text(p$x[,1], p$x[,2], labels=row.names(p$x), pos=c(1,3,4,2), font=2)
#text(l.x, l.y, labels=row.names(p$rotation), col="red", pos=l.pos)
# Add legend
#legend("topleft", legend=c("Tislit", "Sidi Ali", "Michliffen"), title="Sample Area", inset=c(0.01, 0.01), col="black", 
#       pt.bg=c("skyblue2", "gold", "green2"), pch=c(21, 22, 24), pt.cex=2, cex=1.5, bty="n")
legend("topleft", legend=c("0", "1", "2","3","4"), col="black", pt.bg=c("black"," skyblue2", "darkorchid1", "gold", "green2"), pch=25, pt.cex=1.5)





library("factoextra")
# Create groups
group <- c(rep("0", times=2), rep("1", times=16), rep("2", times=16), rep("3", times=14), rep("4", times=14))
# Plot
thing = fviz_pca_biplot(p, axes=c(1,2), label="none", 
                repel=T, 
                pointsize=6, 
                pointshape=21,
                col.var="grey60", 
                fill.ind=group,
                arrowsize=0.6, 
                labelsize=0.5, 
                col.ind=group, 
                palette=c("black", "red", "orange", "green", "blue"), 
                addEllipses=TRUE, ellipse.type="confidence")


png(file="2835x1337 pca 1v2 grey.png", height=8, width=11, units="in", res=220)
print(thing)
dev.off()




