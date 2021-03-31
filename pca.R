
pcadata <- read.csv("pca.csv")

names <- c(rep("CNG30", 4), rep("TX30", 4), rep("WT30", 3))

pcadata$genotype <- names

mycolors <- c('royalblue1', 'darkcyan', 'oldlace')

library(rgl)
library(plot3D)
library(car)
library("RColorBrewer")
colors <- brewer.pal(n=3, name="Dark2")


scatter3d(x = pcadata$PCA.1, y = pcadata$PCA.2, z = pcadata$PCA.3,
          groups = as.factor(pcadata$genotype),
          grid = FALSE, surface = FALSE,
          axis.col = c("black", "black", "black"),
          surface.col = colors,
          xlab = "PC1 - 30%", ylab = "PC2 - 17%",
          zlab = "PC3 - 14%",
          labels = pcadata$genotype,
          axis.scales = FALSE)
text3d(x=1.1, y=c(.9,1,1.1), z=1.1, text = unique(pcadata$genotype) ,col="black")
points3d(x=1.2,y=c(.9,1,1.1),z=1.1, col = colors, size=5)
