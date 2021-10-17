#first make equal area grids....##########
mossDF <- read.xlsx("mossDF.xlsx")
dat <- data.frame(species = mossDF$speciesFinal, longitude = mossDF$lon, latitude = mossDF$lat)
# Define projections
wgs1984 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
# reproject the occurrence records
repro <- dat[, c("longitude", "latitude")] %>% SpatialPoints(proj4string = wgs1984) %>% 
  spTransform(behr)
# create a aqual area template in behrman projection
extent(repro)

be <- raster(ncol = 180, nrow = 142/2, ext = extent(repro), crs = behr)
extent(be)
extent(repro) == extent(be)

#be <- crop(be, extent(repro))

pts <- data.frame(dat$species, coordinates(repro))

# Equal area Occurrence number
eq_occgri <- RichnessGrid(x = pts, ras = be, type = "spnum")

plot(eq_occgri)
polys1 = rasterToPolygons(eq_occgri, na.rm=T)
polys3 <- polys1
str(polys3@data[["layer"]])
polys3@data[["layer"]] <- 1:length(polys3@data[["layer"]])
str(polys3@data[["layer"]])
raster::shapefile(polys3, "eq_occgri.shp", overwrite=F)

#check
mossDF <- cbind(mossDF[,-c(17,18)], repro@coords)
#write.xlsx(mossDF, "mossDF.xlsx", asTable=F)
check <- mossDF[sample(nrow(mossDF), 200), c(1,4,5,17,18)]

ggplot() + geom_map(data = world.inp, map = world.inp, aes(x = long, y = lat, map_id = region), fill = "grey80") + 
  geom_point(data = check, aes(x = lon, y = lat), colour = "darkblue", size = 1) + coord_fixed() + theme_bw() + 
  theme(axis.title = element_blank())

world.inp
r <- world.inp[, c("long", "lat")] %>% SpatialPoints(proj4string = wgs1984) %>% 
  spTransform(behr)
world.inp2 <- cbind(world.inp, r)
world.inp2 <- world.inp2[, -c(1:2)]


ggplot() + geom_map(data = world.inp2, map = world.inp2, aes(x = long, y = lat, map_id = region), fill = "grey80") + 
  geom_point(data = check, aes(x = longitude, y = latitude), colour = "darkblue", size = 1) + coord_fixed() + theme_bw() + 
  theme(axis.title = element_blank())

#get the data#######
cells<-readOGR("eq_occgri.shp")
cellsdat <- cells@data
rownames(cellsdat) <- cellsdat$layer
cells@data <- cellsdat
mossDF <- read.xlsx("mossDF.xlsx")



coords <- mossDF[,c(17,18)] %>% as.matrix()
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
point_shp <- SpatialPoints(coords,proj4string=behr)
pointsSPDF <- SpatialPointsDataFrame(coords=point_shp,data=mossDF)
#writeOGR(pointsSPDF, "points", layer="points2", driver="ESRI Shapefile")

#check crs: should be "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
cells@proj4string
pointsSPDF@proj4string

#do the overlay
over <- over(pointsSPDF, cells)
over$species <- pointsSPDF@data$treeNames
over2 <- over %>% distinct()
write.csv(over2, "over_raw.csv", row.names=F)

#transform the data to a community matrix
data <- unique(over) %>% filter(!is.na(layer))
names(data)
data$count=1
head(data)
anyNA(data$layer)
data %>% filter(is.na(layer))

com1=dcast(data=data,layer~species,value.var="count")
com2 <- com1
com2[is.na(com2)] <- 0
anyNA(com2)
write.xlsx(com2, "comm_matrix_raw.xlsx", asTable = F)




comFull <- over %>% filter(!is.na(layer))
#freq <- frequency(comFull)
comFull=dcast(data=comFull,layer~species)
anyNA(comFull)
write.xlsx(comFull, "comm_matrix_full_raw.xlsx", asTable = F)


#PD Test##########
rm(list = ls())
library(ape)
library(dplyr)
library(openxlsx)
library(picante)
comm <- read.xlsx("comm_matrix_raw.xlsx") %>% dplyr::select(-"Plagiothecium_novaeseelandiae")
max(comm[,2:ncol(comm)])
rownames(comm) <- comm$layer
commLayers <- comm$layer
comm <- comm[,-1]
comm[1:5, 1:5]
treeOG <- read.tree("moss_chrono_noboots.tre")
mossDF <- read.xlsx("mossDF.xlsx")

treenames <- treeOG[["tip.label"]]
commnames <- unique(colnames(comm))
mossDFnames <- unique(mossDF$treeNames)
mossDFnames[! mossDFnames%in%commnames] 
commnames[! commnames%in%mossDFnames] 
commnames[! commnames%in%treenames] 
drop <- treenames[! treenames%in%commnames] 
tree <- drop.tip(treeOG, drop)
tree

prunedTree <- prune.sample(comm,tree)
pd <- pd(comm, prunedTree)
pd$layer <- as.numeric(rownames(pd))
write.csv(pd, "pd_sr.csv", row.names=F)

phydist <- cophenetic(tree)
mpd <- mpd(comm, phydist)
View(mpd)

metrics <- pd
metrics$MPD <- mpd
metrics$MPD_0 <- mpd
metrics$MPD_0[is.na(metrics$MPD)] <- 0
str(metrics)
#ses.pd.result <- ses.pd(comm, tree, null.model = "taxa.labels", runs = 1000)
#ses.pd.result <- read.csv("pdRedo_ses.pd.csv")
#ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels",  abundance.weighted = FALSE, runs = 1000)
#####Get rid of Cells with too few spp ##########
comm <- read.xlsx("comm_matrix_raw.xlsx") %>% dplyr::select(-"Plagiothecium_novaeseelandiae")
pd <- read.csv("pd_sr.csv")
hist(pd$SR, breaks=50)
comFull <- read.xlsx("comm_matrix_full_raw.xlsx") %>% dplyr::select(-"Plagiothecium_novaeseelandiae")
sums <- comFull %>% filter(rowSums(.[2:ncol(comFull)]) < 50)
sums <- rowSums(sums[2:ncol(sums)])
hist(sums, breaks=50)

pd <- pd[pd$SR > 5,] 
write.csv(pd, "pd_sr_5.csv", row.names=F)
getridofme <- pd %>% dplyr::select(layer)
comm_done <- left_join(getridofme, comm)
comFull_done <- left_join(getridofme, comFull)
sums <- comFull_done %>% filter(rowSums(.[2:ncol(comFull_done)]) < 50)
sums <- rowSums(sums[2:ncol(sums)])
hist(sums, breaks=50)
write.xlsx(comFull_done, "comm_matrix_full.xlsx", asTable = F)
write.xlsx(comm_done, "comm_matrix.xlsx", asTable = F)


#
comm <- read.xlsx("comm_matrix.xlsx") 
comm <- comm[,colSums(comm) > 0]
comm2 <- comm[2:ncol(comm)]
tree <-read.tree("moss_chrono_noboots.tre")
prunedTree <- prune.sample(comm2,tree)

pd <- ses.pd(comm2, prunedTree, prunedTree, null.model = "taxa.labels", runs = 1000)
pd$layer <- comm$layer
write.csv(pd, "sesPD5.csv", row.names=F)

phydist <- cophenetic(prunedTree)
mpd<- ses.mpd(comm2, phydist, null.model = "taxa.labels", runs = 1000)
mpd$layer <- comm$layer
####interpolation, rarefaction for SR############
library(iNEXT)
comFull <- read.xlsx("comm_matrix_full.xlsx")
sums <- comFull %>% filter(rowSums(.[2:ncol(comFull)]) < 200)
sums <- rowSums(sums[2:ncol(sums)])
hist(sums, breaks=1000)


t <- t(comFull[2:ncol(comFull)])
m <- c(1, 2,10,20,50, 100, 200, 500, 5000, max(rowSums(comFull)))
inex <- iNEXT(t, datatype = "abundance", size = m)
length(inex[["iNextEst"]])

inex_df <- inex[["iNextEst"]][[1]]
inex_df$layer <- paste(comFull$layer[1])
for(i in 2:length(inex[["iNextEst"]])){
  df.tmp <- inex[["iNextEst"]][[i]]
  df.tmp$layer <- paste(comFull$layer[i])
  inex_df <- bind_rows(inex_df, df.tmp)
  print(i)
}
inex_df <- inex_df %>% distinct()
inex_df <- inex_df %>% mutate(layer = as.numeric(layer))
write.csv(inex_df, "inexOut.csv", row.names=F)

inex_df %>% filter(m==500)

i=1
inex_df2 <- inex_df[inex_df$m==m[i],]
colnames(inex_df2) <- c("m", paste("method", m[i], sep="_"),  paste("order", m[i], sep="_"),  paste("qD", m[i], sep="_"),  paste("qD.LCL", m[i], sep="_"),  paste("qD.UCL", m[i], sep="_"),  paste("SC", m[i], sep="_"),  paste("SC.LCL", m[i], sep="_"),  paste("SC.UCL", m[i], sep="_"), "layer")
inex_df2 <- inex_df2 %>% dplyr::select(-m)
for(i in 2:length(m)){
  df.tmp <- inex_df[inex_df$m==m[i],]
  colnames(df.tmp) <- c("m", paste("method", m[i], sep="_"),  paste("order", m[i], sep="_"),  paste("qD", m[i], sep="_"),  paste("qD.LCL", m[i], sep="_"),  paste("qD.UCL", m[i], sep="_"),  paste("SC", m[i], sep="_"),  paste("SC.LCL", m[i], sep="_"),  paste("SC.UCL", m[i], sep="_"), "layer")
  nrow(df.tmp) == 3087
  df.tmp <- df.tmp %>% dplyr::select(-m)
  inex_df2 <- full_join(inex_df2, df.tmp)
  print(i)
}

anyNA(inex_df2)
write.csv(inex_df2, "inexOut_cols.csv", row.names=F)
method <- inex_df2 %>% select(method_200, layer)
method$method <- NA
method$method[method$method_200 == "interpolated"] <- 1
method$method[method$method_200 == "extrapolated"] <- 2
method$method[method$method_200 == "observed"] <- 3
write.csv(method, "inex_method.csv", row.names=F)



###interpolation, rarefaction for PD######
library(iNEXTPD2)
library(ape)
library(openxlsx)
library(picante)
treeOG <- read.tree("moss_chrono_noboots.tre")
comFull <- read.xlsx("comm_matrix_full.xlsx")
tree <- prune.sample(comFull,treeOG)
data <- as.data.frame(t(comFull[2:ncol(comFull)]))

bad <- NULL
for(i in 1:length(data)){
  print(i)
  data3 <- data[,c(i, i+1)]
  inex <- iNEXTPD(data3, tree=tree, datatype = "abundance", size = 200, nboot = 0)
  ans <- is.null(inex)
  if (ans==TRUE){
    print(ans)
    bad <- c(bad,i)
  }
}

data3 <- data[,-c(bad)]

# returns NULL if nboot != 0
inex <- iNEXTPD(data= data3, tree=tree, datatype = "abundance", size = 200, nboot = 0)
inex_nullSize <- iNEXTPD(data= data3, tree=tree, datatype = "abundance", nboot = 0)
#inex_boot <- iNEXTPD(data= data3, tree=tree, datatype = "abundance", size = 200, nboot = 5)


# okay let/s see if we can identify which assemblages are causing problems
badtoo <- NULL
for(i in 1:length(data)){
  print(i)
  data3 <- data[,c(i, i+1)]
  inex <- iNEXTPD(data3, tree=tree, datatype = "abundance", size = 200)
  ans <- is.null(inex)
  if (ans==TRUE){
    print(ans)
    badtoo <- c(badtoo,i)
  }
}
badtoo <- unique(badtoo)
data4 <- data[,-c(badtoo)]
inex_boot <- iNEXTPD(data=data4, tree=tree, datatype = "abundance", size = 200)

#output files:
data_fix <- data[1,]
data_fix[2,] <- comFull$layer
data_fix2 <- data_fix[,-c(badtoo)]
colz <- as.numeric(data_fix2[2,])
result_200$layer <- colz
write.csv(result_200, "inextPD_200.csv", row.names = F)

df <- full_join(data.frame(SR_rare = inextSR$qD,
                           layer = inextSR$layer),
                data.frame(PD_rare = result_200$qPD,
                           layer = result_200$layer))

pd <- read.csv("pd_sr_5.csv")
df2 <- full_join(df, pd)
write.csv(df2, "pdSR_rares.csv", row.names = F)


#WE and CWE############
#comm <- read.xlsx("comm_matrix_raw.xlsx") %>% dplyr::select(-"Plagiothecium_novaeseelandiae")
comm <- read.xlsx("comm_matrix.xlsx") 
comm <- comm[,colSums(comm) > 0]

weighted.endemism <- function(x){
  df <-reshape2::melt(x, id.vars="layer") %>% filter(value == 1) %>% select(-"value") %>% distinct() %>% arrange(layer)
  colnames(df) <- c("layer", "species")
  range <- data.frame(table(df$species))
  colnames(range) <- c("species", "range")
  df <- full_join(df, range)
  WE <- df %>% 
    group_by(layer) %>% 
    summarise(WE = sum(1/range))
  WE
}

SpecRich <- function(x){
  mm1 <- data.frame(table(x$layer))
  names(mm1) <- c("layer", "SR")
  mm1
}

corrected.weighted.endemism <- function(x){
  df <-reshape2::melt(x, id.vars="layer") %>% filter(value == 1) %>% select(-"value") %>% distinct() %>% arrange(layer)
  colnames(df) <- c("layer", "species")
  tmp <- SpecRich(df) %>% arrange(layer)
  we <- weighted.endemism(x) %>% arrange(layer)
  we$CWE <- we$WE/tmp$SR
  we
}

WE <- weighted.endemism(comm)
CWE <- corrected.weighted.endemism(comm)
pd <- read.csv("pd_sr_5.csv")
metrics <- full_join(pd, CWE)
anyNA(metrics)

##WE and CWE rand######
ses.we <- function(samp, null.model = c( "richness", "frequency", "independentswap", "trialswap"),
                   runs = 999, iterations = 1000){
  we.obs <- as.vector(weighted.endemism(samp)$WE)
  null.model <- match.arg(null.model)
  we.rand <- switch(null.model,
                    richness = t(replicate(runs, as.vector(weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="richness")))$WE))),
                    frequency = t(replicate(runs, as.vector(weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="frequency")))$WE))),
                    independentswap  = t(replicate(runs, as.vector(weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="independentswap", iterations)))$WE))),
                    trialswap = t(replicate(runs, as.vector(weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="trialswap", iterations)))$WE)))
  )
  
  we.rand.mean <- apply(X = we.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  we.rand.sd <- apply(X = we.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  we.obs.z <- (we.obs - we.rand.mean)/we.rand.sd
  we.obs.rank <- apply(X = rbind(we.obs, we.rand), MARGIN = 2,
                       FUN = rank)[1, ]
  we.obs.rank <- ifelse(is.na(we.rand.mean),NA,we.obs.rank)
  return(test <- data.frame(ntaxa=specnumber(samp[2:nrow(samp)]),we.obs, we.rand.mean, we.rand.sd, we.obs.rank,
                            we.obs.z, we.obs.p=we.obs.rank/(runs+1),runs=runs, layer = samp$layer))
  
  
  
}
sesWE <- ses.we(comm, null.model="frequency", runs = 1000)

ses.cwe <- function(samp, null.model = c( "richness", "frequency", "independentswap", "trialswap"),
                    runs = 999, iterations = 1000){
  cwe.obs <- as.vector(corrected.weighted.endemism(samp)$CWE)
  null.model <- match.arg(null.model)
  cwe.rand <- switch(null.model,
                     richness = t(replicate(runs, as.vector(corrected.weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="richness")))$CWE))),
                     frequency = t(replicate(runs, as.vector(corrected.weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="frequency")))$CWE))),
                     independentswap  = t(replicate(runs, as.vector(corrected.weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="independentswap", iterations)))$CWE))),
                     trialswap = t(replicate(runs, as.vector(corrected.weighted.endemism( data.frame(layer = as.numeric(rownames(samp)), randomizeMatrix(samp[2:ncol(samp)],null.model="trialswap", iterations)))$CWE)))
  )
  
  cwe.rand.mean <- apply(X = cwe.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  cwe.rand.sd <- apply(X = cwe.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  cwe.obs.z <- (cwe.obs - cwe.rand.mean)/cwe.rand.sd
  cwe.obs.rank <- apply(X = rbind(cwe.obs, cwe.rand), MARGIN = 2,
                        FUN = rank)[1, ]
  cwe.obs.rank <- ifelse(is.na(cwe.rand.mean),NA,cwe.obs.rank)
  return(test <- data.frame(ntaxa=specnumber(samp[2:nrow(samp)]),cwe.obs, cwe.rand.mean, cwe.rand.sd, cwe.obs.rank,
                            cwe.obs.z, cwe.obs.p=cwe.obs.rank/(runs+1),runs=runs, layer = samp$layer))
  
  
  
}

sesCWE <- ses.cwe(comm, null.model="frequency", runs = 1000)

write.csv(sesWE, "WE_5.csv", row.names = F)
write.csv(sesCWE, "CWE_5.csv", row.names = F)
####PE and CPWE####################
rm(list = ls())
#comm <- read.xlsx("comm_matrix_raw.xlsx") %>% dplyr::select(-"Plagiothecium_novaeseelandiae")
comm <- read.xlsx("comm_matrix.xlsx")
comm <- comm[,colSums(comm) > 0]

library(canaper)
library(picante)
comm2 <- comm[2:ncol(comm)]
tree <-read.tree("moss_chrono_noboots.tre")
prunedTree <- prune.sample(comm2,tree)

library(phyloregion)
sparse <- phyloregion::dense2sparse(comm2)
pe_obs <- phylo_endemism(sparse, prunedTree, 
                         weighted = TRUE)

pe_df <- data.frame(layer = comm$layer, pe = pe_obs)
write.csv(pe_df, "PE.csv", row.names=F)


pcwe_obs <- function(sparse, y, tree){
  pe <- phylo_endemism(sparse, tree, weighted = TRUE)
  df <- pd(y, tree, include.root=T)
  df$pe <- pe
  CWE <- df$pe/df$PD
  CWE
}

pddf <- read.csv("pd_sr_5.csv") %>% select(layer, PD) %>% arrange(layer)
pcwe_obs <- function(sparse, y, pd, tree){
  pe <- phylo_endemism(sparse, tree, weighted = TRUE)
  CWE <- pe/pd
  CWE
}

pcwe.obs <- pcwe_obs(sparse, comm2, prunedTree)
pcwe_df <- data.frame(layer = comm$layer, pcwe = pcwe.obs)
write.csv(pcwe_df, "PCWE.csv", row.names=F)


###PE and CPWE rand ########
ses.pe <- function(samp, tree, null.model = c("taxa.labels", "richness", "frequency",
                                              "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), runs = 999, iterations = 1000){
  pe.obs <- as.vector(phyloregion::phylo_endemism(samp, tree, weighted=TRUE))
  null.model <- match.arg(null.model)
  pe.rand <- switch(null.model,
                    taxa.labels = t(replicate(runs, as.vector(phyloregion::phylo_endemism(samp, tipShuffle(tree), weighted=TRUE)))),
                    richness = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="richness"), tree, weighted=TRUE)))),
                    frequency = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="frequency"), tree, weighted=TRUE)))),
                    sample.pool = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="richness"), tree, weighted=TRUE)))),
                    phylogeny.pool = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="richness"),
                                                                                             tipShuffle(tree), weighted=TRUE)))),
                    independentswap = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="independentswap", iterations), tree, weighted=TRUE)))),
                    trialswap = t(replicate(runs, as.vector(phyloregion::phylo_endemism(randomizeMatrix(samp,null.model="trialswap", iterations), tree, weighted=TRUE))))
  )
  pe.rand.mean <- apply(X = pe.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  pe.rand.sd <- apply(X = pe.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  pe.obs.z <- (pe.obs - pe.rand.mean)/pe.rand.sd
  pe.obs.rank <- apply(X = rbind(pe.obs, pe.rand), MARGIN = 2,
                       FUN = rank)[1, ]
  pe.obs.rank <- ifelse(is.na(pe.rand.mean),NA,pe.obs.rank)
  return(data.frame(ntaxa=specnumber(samp),pe.obs, pe.rand.mean, pe.rand.sd, pe.obs.rank,
                    pe.obs.z, pe.obs.p=pe.obs.rank/(runs+1),runs=runs, row.names = row.names(samp)))
  
  
}

#sesPE <- ses.pe(sparse, prunedTree, null.model = "taxa.labels", runs = 1000)
sesPE <- ses.pe(sparse, prunedTree, null.model = "richness", runs = 1000)
sesPE$layer <- comm$layer
write.csv(sesPE, "PE_1000_richness.csv", row.names = F)

ses.cwpe <- function(sparse, y, tree, null.model = c("taxa.labels", "richness", "frequency",
                                                     "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
                     runs = 999, iterations = 1000){
  cwpe.obs <- as.vector(pcwe_obs(sparse, y, tree))
  null.model <- match.arg(null.model)
  cwpe.rand <- switch(null.model,
                      taxa.labels = t(replicate(runs, as.vector(pcwe_obs(sparse, y, tipShuffle(tree))))),
                      richness = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="richness"), tree)))),
                      frequency = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="frequency"), tree)))),
                      sample.pool = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="richness"), tree)))),
                      phylogeny.pool = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="richness"),
                                                                            tipShuffle(tree))))),
                      independentswap = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="independentswap", iterations), tree)))),
                      trialswap = t(replicate(runs, as.vector(pcwe_obs(randomizeMatrix(samp,null.model="trialswap", iterations), tree))))
  )
  cwpe.rand.mean <- apply(X = cwpe.rand, MARGIN = 2, FUN = mean, na.rm=TRUE)
  cwpe.rand.sd <- apply(X = cwpe.rand, MARGIN = 2, FUN = sd, na.rm=TRUE)
  cwpe.obs.z <- (cwpe.obs - cwpe.rand.mean)/cwpe.rand.sd
  cwpe.obs.rank <- apply(X = rbind(cwpe.obs, cwpe.rand), MARGIN = 2,
                         FUN = rank)[1, ]
  cwpe.obs.rank <- ifelse(is.na(cwpe.rand.mean),NA,cwpe.obs.rank)
  return(data.frame(ntaxa=specnumber(samp),cwpe.obs, cwpe.rand.mean, cwpe.rand.sd, cwpe.obs.rank,
                    cwpe.obs.z, cwpe.obs.p=cwpe.obs.rank/(runs+1),runs=runs, row.names = row.names(samp)))
  
  
}
sesCWPE <- ses.cwpe(sparse, comm2, prunedTree, null.model = "taxa.labels", runs = 1000)
sesCWPE$layer <- comm$layer
write.csv(sesCWPE, "sesCWPE.csv", row.names = F)

