library(dplyr)
library(phyloregion)
library(ape)
library(rgdal)
library(openxlsx)
rm(list = ls())
comm <- read.xlsx("comm_matrix.xlsx")
comm <- comm[,colSums(comm) > 0]
comm2 <- comm[2:ncol(comm)]
rownames(comm2) <- comm$layer
sparse <- phyloregion::dense2sparse(comm2)
tree <-read.tree("moss_chrono_noboots.tre")
prunedTree <- prune.sample(comm2,tree)
phylob <- phylobeta(sparse, prunedTree)

df <- dist(phylob[[1]])
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

cells<-readOGR("eq_occgri.shp")
cellsdat <- cells@data
rownames(cellsdat) <- cellsdat$layer
colnames(cellsdat) <- "grids"
cells@data <- cellsdat
y <- phyloregion(phylob[[1]], shp=cells, k=6)
plot(y, palette="NMDS")

y <- select_linkage(phylob[[1]])
barplot(y, horiz = TRUE, las = 1)
y <- phyloregion(phylob[[1]], shp=cells, k=6, method = "ward.D2" )
plot(y, palette="NMDS")
plot_NMDS(y, cex=3)
text_NMDS(y)

shape <- y[["shp"]]
writeOGR(shape, "phyloregion", layer="phyloregions", driver="ESRI Shapefile")
info <- y[["region.df"]]
write.csv(info, "~/Documents/bigBryogeography/round2/phyloregion/info.csv", row.names = F)
att <-  y[["shp"]]@data
write.csv(att, "~/Documents/bigBryogeography/round2/phyloregion/attributes.csv", row.names = F)
