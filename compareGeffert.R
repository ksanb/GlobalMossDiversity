#geffert shit
rm(list = ls())
library(picante)
library(openxlsx)
library(ggplot2)
library(raster)
library(rgdal)
library(reshape2)
library(dplyr)
#

#combining datasets - matching my names to his#######
d4 <-read.table("~/Documents/bigBryogeography/round2/GeffertData/Dist_Rank4.txt", header=T) 
df <- cbind(d4, data.frame(TropID_Final=NA,
                           speciesFinal=NA,
                           sciNameFinal=NA,
                           ID_syn=NA,
                           species_syn=NA,
                           sciName_syn=NA,
                           treeNames=NA))
InSyns <-data.frame(species_syn = word(d4$FullName, 1,2))
InSyns <- left_join(InSyns, syns) %>% distinct()
for(i in 1:nrow(df)){
  spec <- word(df$FullName[i], 1,2)
  if(!is.na(InSyns$homonym[InSyns$species_syn==spec])){
    if(InSyns$homonym[InSyns$species_syn==spec] == FALSE){
      tmp <- InSyns %>% filter(species_syn==spec)
      df$TropID_Final[i] <- paste(tmp$TropID_Final)
      df$speciesFinal[i] <- paste(tmp$speciesFinal)
      df$sciNameFinal[i] <- paste(tmp$sciNameFinal)
      df$ID_syn[i] <- paste(tmp$ID_syn)
      df$species_syn[i] <- paste(tmp$species_syn)
      df$sciName_syn[i] <- paste(tmp$sciName_syn)
      df$treeNames[i] <- paste(tmp$treeNames)
    }
    if(InSyns$homonym[InSyns$species_syn==spec] == TRUE){
      sci <- df$FullName[i]
      tmp <- InSyns %>% filter(sciName_syn==sci)
      if(nrow(tmp) > 0){
        df$TropID_Final[i] <- paste(tmp$TropID_Final)
        df$speciesFinal[i] <- paste(tmp$speciesFinal)
        df$sciNameFinal[i] <- paste(tmp$sciNameFinal)
        df$ID_syn[i] <- paste(tmp$ID_syn)
        df$species_syn[i] <- paste(tmp$species_syn)
        df$sciName_syn[i] <- paste(tmp$sciName_syn)
        df$treeNames[i] <- paste(tmp$treeNames)
      }
    }
  }
  
  print(i)
}

data_save <- df %>% filter(!is.na(TropID_Final))
write.csv(data, "GeffertDat_D4.csv", row.names = F)



#
#redoing geffert:#############
d4 <-read.table("~/Documents/bigBryogeography/round2/GeffertData/Dist_Rank4.txt", header=T)  %>% 
  dplyr::select(FullName, AreaCode) %>% 
  filter(!is.na(AreaCode)) %>% distinct()
d4$AreaCode[70853] <- "ZW"
d4 <- d4 %>% distinct()
OGUshp <- readOGR('~/Documents/Shapefiles/OGU.shp')
data <- d4
data$count = 1
anyNA(data)
com1=dcast(data=data,AreaCode~FullName,value.var="count") 
com2 <- com1
com2[is.na(com2)] <- 0
anyNA(com2)
max(rowSums(com2[,2:ncol(com2)]))
write.xlsx(com2, "comm_matrix_geffert_OG_d4.xlsx", asTable = F)
rownames(com2) <- com2$AreaCode
comm <- com2[,-1]
sr <- data.frame(rowSums(comm))
colnames(sr) <- "SR"
anyNA(sr)
write.csv(sr, "sr_geffert_OG_d4.csv", row.names=T)

#using the combined datasets#######
rm(list = ls())
OGUshp <- readOGR('~/Documents/Shapefiles/OGU.shp')
data <- read.csv('GeffertDat_d4.csv') %>% dplyr::select(treeNames, AreaCode) %>% 
  filter(!is.na(AreaCode))%>% distinct()
#data$AreaCode[data$AreaCode == "ZW "] <- "ZW"
data <- data %>% distinct()
data$count = 1
com1=dcast(data=data,AreaCode~treeNames,value.var="count") 
com2 <- com1
com2[is.na(com2)] <- 0
anyNA(com2)

write.xlsx(com2, "comm_matrix_geffert_treespp_d4.xlsx", asTable = F)
rownames(com2) <- com2$AreaCode
comm <- com2[,-1]
treeOG <- read.tree("moss_chrono_noboots.tre")
prunedTree <- prune.sample(comm,treeOG)
pd <- pd(comm, prunedTree)
anyNA(pd)
write.csv(pd, "pd_sr_geffert_treespp_d4.csv", row.names=T)

phydist <- cophenetic(prunedTree)
mpd <- mpd(comm, phydist)

length(mpd)
length(rownames(comm))
mpddf <- data.frame(MPD = mpd, AreaCode = rownames(comm)) %>% filter(!is.na(MPD))
write.csv(mpddf, "mpd_geffert_treespp_d4.csv", row.names=F)


pd$mpd <- mpd
anyNA(pd$mpd)
write.csv(pd, "mpd_pd_sr_geffert_treespp_d4.csv", row.names=T)


#
#regressions: tree spp and OG dataset####
rm(list = ls())
pd <- read.csv("pd_sr_geffert_treespp_d4.csv")
colnames(pd) <- c("AreaCode", "PD", "SR")
sr <- read.csv("sr_geffert_OG_d4.csv")
colnames(sr) <- c("AreaCode", "OG")
dat <- full_join(pd, sr)
anyNA(dat$OG)

length(unique(data$treeNames)) #2622
length(unique(dLegit$FullName)) #9210
# test for correlation between original species richness and tree species richness,
#first just plot:
library(ggplot2)
ggplot(dat) + geom_point(aes(x=SR, y = OG)) + xlab("Richness (Tree)") + ylab("Richness (Geffert Data)") 
ggplot(dat) + geom_point(aes(x=SR, y = PD)) + xlab("Richness (Tree)") + ylab("PD") 
ggplot(dat) + geom_point(aes(x=PD, y = OG)) + xlab("PD") + ylab("Richness (Geffert Data)") 

r <- lm(OG ~ SR, dat)
summary(r)
x <- dat$SR
y <- dat$OG
cor.test(x,y) #0.9581355 

x <- dat$SR
y <- dat$PD
cor.test(x,y) #0.9781879

x <- dat$OG
y <- dat$PD
cor.test(x,y) #0.9374088 

m <- lm(OG ~ PD, dat)
summary(m)
res <- resid(m)
plot(fitted(m), res)
abline(0,0)
#The x-axis displays the fitted values and the y-axis displays the residuals.
m <- lm(OG ~ PD, dat)
summary(m)
res <- resid(m)
plot(fitted(m), res)
abline(0,0)

OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()
ggplot(dat2) + geom_point(aes(x=SR, y = OG, color =Continent)) + xlab("Richness (Tree)") + ylab("Richness (Geffert Data)")

m <- lm(OG ~ PD, dat2)
summary(m)
res <- resid(m)
plotdat <- cbind(dat2, res)
plotdat$fit <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit, y = res, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - PD")

m <- lm(OG ~ SR, dat2)
summary(m)
res_SR <- resid(m)
plotdat <- cbind(dat2, res_SR)
plotdat$fit_SR <-  fitted(m)
plotdat <- plotdat %>% filter(!is.na(Continent))
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~SR (tree spp.)") +
  annotate("text", x=700, y=350, label= "Corr = 0.96")

labelz <- c("Antarctica", 
            "Africa", 
            "Asia", 
            "SE Asia", 
            "Australia", 
            "Central America", 
            "Europe", 
            "North America",
            "South America")
#valuez <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#6A869F", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 9, name = 'Paired')
set.seed(001) 
valuez <- sample(brewer.pal(n = 9, name = 'Paired'))


ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) +
  # xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~SR") +
  annotate("text", x=700, y=350, label= "Corr = 0.96") + labs(title = "Residuals vs. Fitted: OG ~SR\n", x = "Fitted Values", y = "Residuals", color = "Continent\n") +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent), size = 2) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()

#geff vs treespp
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()

pdf(file = "~/Documents/bigBryogeography/round2/geffFigz/regressions/geff_vs_treespp.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
dev.off()



#ggplot(plotdat, aes(x = SR, y = OG)) + geom_point(aes(color =Continent)) + geom_abline(slope = coef(m)[[2]], intercept = coef(m)[[1]])
ggplot(plotdat) + geom_point(aes(x=SR, y = PD, color =Continent)) + xlab("Richness (Tree)") + ylab("PD") 
ggplot(plotdat) + geom_point(aes(x=SR, y = OG, color =Continent)) + xlab("Richness (Tree)") + ylab("Richness (Geffert Data)") 
ggplot(plotdat) + geom_point(aes(x=PD, y = OG, color =Continent)) + xlab("PD") + ylab("Richness (Geffert Data)") 

resDat <- plotdat
resDat$resMinusFit_OGSR <- resDat$res_SR - resDat$fit_SR
resDat$fitMinusRes_OGSR <- resDat$fit_SR - resDat$res_SR

m <- lm(OG ~ PD, dat2)
summary(m)
res_PD <- resid(m)
resDat <- cbind(resDat, res_PD)
resDat$fit_PD <-  fitted(m)
resDat$resMinusFit_OGPD <- resDat$res_PD - resDat$fit_PD
resDat$fitMinusRes_OGPD <- resDat$fit_PD - resDat$res_PD
write.csv(resDat, "treespp_OG_residuals_PD_SR_d4.csv", row.names=F)

#regressions gbif vs OG dataset####
rm(list = ls())
pd <- read.csv("pd_sr_geffert_gbif.csv")
colnames(pd) <- c("AreaCode", "PD", "SR")
sr <- read.csv("sr_geffert_OG_d4.csv")
colnames(sr) <- c("AreaCode", "OG")
dat <- full_join(pd, sr)
anyNA(dat$OG)

length(unique(data$treeNames)) #2622
length(unique(dLegit$FullName)) #9210
# test for correlation between original species richness and tree species richness,
#first just plot:
library(ggplot2)
ggplot(dat) + geom_point(aes(x=SR, y = OG)) + xlab("Richness (GBIF)") + ylab("Richness (Geffert Data)") 
ggplot(dat) + geom_point(aes(x=SR, y = PD)) + xlab("Richness (GBIF)") + ylab("PD") 
ggplot(dat) + geom_point(aes(x=PD, y = OG)) + xlab("PD") + ylab("Richness (Geffert Data)") 

r <- lm(OG ~ SR, dat)
summary(r)
x <- dat$SR
y <- dat$OG
cor.test(x,y) #0.6651194

x <- dat$SR
y <- dat$PD
cor.test(x,y) #0.9796159 

x <- dat$OG
y <- dat$PD
cor.test(x,y) #0.6652471 

m <- lm(OG ~ PD, dat)
summary(m)
res <- resid(m)
plot(fitted(m), res)
abline(0,0)
#The x-axis displays the fitted values and the y-axis displays the residuals.
m <- lm(OG ~ PD, dat)
summary(m)
res <- resid(m)
plot(fitted(m), res)
abline(0,0)

OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()
ggplot(dat2) + geom_point(aes(x=SR, y = OG, color =Continent)) + xlab("Richness (GBIF)") + ylab("Richness (Geffert Data)")

m <- lm(OG ~ PD, dat2)
summary(m)
res <- resid(m)
plotdat <- dat2 %>% filter(!is.na(OG)) %>% filter(!is.na(PD))
plotdat <- cbind(plotdat, res)
plotdat$fit <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit, y = res, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - PD")
ggplot(plotdat) + geom_point(aes(x= fit, y = res, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~ PD (gbif)") +
  annotate("text", x=600, y=350, label= "Corr = 0.62")


m <- lm(OG ~ SR, dat2)
summary(m)
res_SR <- resid(m)
plotdat <- dat2 %>% filter(!is.na(OG)) %>% filter(!is.na(SR))
plotdat <- cbind(plotdat, res_SR)
plotdat$fit_SR <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - SR")
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~ SR (gbif)") +
  annotate("text", x=600, y=350, label= "Corr = 0.62")

#ggplot(plotdat, aes(x = SR, y = OG)) + geom_point(aes(color =Continent)) + geom_abline(slope = coef(m)[[2]], intercept = coef(m)[[1]])
ggplot(plotdat) + geom_point(aes(x=SR, y = PD, color =Continent)) + xlab("Richness (Tree)") + ylab("PD") 
ggplot(plotdat) + geom_point(aes(x=SR, y = OG, color =Continent)) + xlab("Richness (Tree)") + ylab("Richness (Geffert Data)") 
ggplot(plotdat) + geom_point(aes(x=PD, y = OG, color =Continent)) + xlab("PD") + ylab("Richness (Geffert Data)") 

resDat <- plotdat
resDat$resMinusFit_OGSR <- resDat$res_SR - resDat$fit_SR
resDat$fitMinusRes_OGSR <- resDat$fit_SR - resDat$res_SR

m <- lm(OG ~ PD, dat2)
summary(m)
res_PD <- resid(m)
resDat <- cbind(resDat, res_PD)
resDat$fit_PD <-  fitted(m)
resDat$resMinusFit_OGPD <- resDat$res_PD - resDat$fit_PD
resDat$fitMinusRes_OGPD <- resDat$fit_PD - resDat$res_PD
write.csv(resDat, "gbif_OG_residuals_PD_SR.csv", row.names=F)

#regressions pt 2: gbif vs. treespp#########
rm(list = ls())
gbifspp <- read.csv("pd_sr_geffert_gbif.csv")
colnames(gbifspp) <- c("AreaCode", "PD_gbif", "SR_gbif")
gbifmpd <- read.csv("mpd_geffert_gbif.csv")
colnames(gbifmpd) <- c("MPD_gbif","AreaCode")
dat <- full_join(gbifmpd, gbifspp)

treespp <- read.csv("mpd_pd_sr_geffert_treespp.csv")
colnames(treespp) <- c("AreaCode", "PD_tree", "SR_tree", "MPD_tree")
dat <- full_join(dat, treespp)

r <- lm(SR_tree ~ SR_gbif, dat)
summary(r)
x <- dat$PD_gbif
y <- dat$PD_tree
cor.test(x,y) #0.6981533

r <- lm(PD_tree ~ PD_gbif, dat)
summary(r)
x <- dat$PD_gbif
y <- dat$PD_tree
cor.test(x,y) #0.6981533

r <- lm(MPD_tree ~ MPD_gbif, dat)
summary(r)
x <- dat$MPD_gbif
y <- dat$MPD_tree
cor.test(x,y) #0.7596206 


OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()

m <- lm(PD_tree ~ PD_gbif, dat)
summary(m)
res_PD <- resid(m)
plotdat <- dat2 %>% filter(!is.na(PD_gbif)) %>% filter(!is.na(PD_tree))
plotdat <- cbind(plotdat, res_PD)
plotdat$fit_PD <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_PD, y = res_PD, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - PD")
ggplot(plotdat) + geom_point(aes(x= fit_PD, y = res_PD, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: PD (tree) ~ PD (gbif)") +
  annotate("text", x=17000, y=10000, label= "Corr = 0.7")


resDat <- plotdat
resDat$resMinusFit_PD <- resDat$res_PD - resDat$fit_PD
resDat$fitMinusRes_PD <- resDat$fit_PD - resDat$res_PD
write.csv(resDat, "treespp_gbif_residuals_PD.csv", row.names=F)

m <- lm(MPD_tree ~ MPD_gbif, dat)
summary(m)
res_MPD <- resid(m)
plotdat <- dat2 %>% filter(!is.na(MPD_gbif)) %>% filter(!is.na(MPD_tree))
plotdat <- cbind(plotdat, res_MPD)
plotdat$fit_MPD <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - MPD")
ggplot(plotdat) + geom_point(aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: MPD (tree) ~ MPD (gbif)") +
  annotate("text", x=450, y=80, label= "Corr = 0.76")

resDat <- plotdat
resDat$resMinusFit_MPD <- resDat$res_MPD - resDat$fit_MPD
resDat$fitMinusRes_MPD <- resDat$fit_MPD - resDat$res_MPD
write.csv(resDat, "treespp_gbif_residuals_MPD.csv", row.names=F)

#
##Plots!!!##############
#plots: geff vs. treesppp##########
rm(list = ls())
library(ggplot2)
pd <- read.csv("pd_sr_geffert_treespp_d4.csv")
colnames(pd) <- c("AreaCode", "PD", "SR")
sr <- read.csv("sr_geffert_OG_d4.csv")
colnames(sr) <- c("AreaCode", "OG")
dat <- full_join(pd, sr)
anyNA(dat$OG)

OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()

m <- lm(OG ~ SR, dat2)
summary(m)
res_SR <- resid(m)
plotdat <- cbind(dat2, res_SR)
plotdat$fit_SR <-  fitted(m)
plotdat <- plotdat %>% filter(!is.na(Continent))
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~SR (tree spp.)") +
  annotate("text", x=700, y=350, label= "Corr = 0.96")

labelz <- c("Antarctica", 
            "Africa", 
            "Asia", 
            "SE Asia", 
            "Australia", 
            "Central America", 
            "Europe", 
            "North America",
            "South America")

valuez <- c("#FDBF6F", "#E31A1C", "#A6CEE3", "#33A02C", "#CAB2D6", "#B2DF8A", "#1F78B4", "#FF7F00", "#FB9A99")

#SR:
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()

pdf(file = "~/Documents/bigBryogeography/round2/geffFigz/regressions/geff_vs_treespp_SR.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
dev.off()

#Plots: geff vs. gbif:#####
rm(list = ls())
pd <- read.csv("pd_sr_geffert_gbif.csv")
colnames(pd) <- c("AreaCode", "PD", "SR")
sr <- read.csv("sr_geffert_OG_d4.csv")
colnames(sr) <- c("AreaCode", "OG")
dat <- full_join(pd, sr)
anyNA(dat$OG)

OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()


m <- lm(OG ~ SR, dat2)
summary(m)
res_SR <- resid(m)
plotdat <- dat2 %>% filter(!is.na(OG)) %>% filter(!is.na(SR))
plotdat <- cbind(plotdat, res_SR)
plotdat$fit_SR <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - SR")
ggplot(plotdat) + geom_point(aes(x= fit_SR, y = res_SR, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: OG ~ SR (gbif)") +
  annotate("text", x=600, y=350, label= "Corr = 0.62")

labelz <- c("Antarctica", 
            "Africa", 
            "Asia", 
            "SE Asia", 
            "Australia", 
            "Central America", 
            "Europe", 
            "North America",
            "South America")

valuez <- c("#FDBF6F", "#E31A1C", "#A6CEE3", "#33A02C", "#CAB2D6", "#B2DF8A", "#1F78B4", "#FF7F00", "#FB9A99")


#SR:
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()

pdf(file = "~/Documents/bigBryogeography/round2/geffFigz/regressions/geff_vs_gbif_SR.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(plotdat, aes(x= fit_SR, y = res_SR, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
dev.off()




#Plots: treespp vs. gbif: #####
rm(list = ls())
gbifspp <- read.csv("pd_sr_geffert_gbif.csv")
colnames(gbifspp) <- c("AreaCode", "PD_gbif", "SR_gbif")
gbifmpd <- read.csv("mpd_geffert_gbif.csv")
colnames(gbifmpd) <- c("MPD_gbif","AreaCode")
dat <- full_join(gbifmpd, gbifspp)

treespp <- read.csv("mpd_pd_sr_geffert_treespp.csv")
colnames(treespp) <- c("AreaCode", "PD_tree", "SR_tree", "MPD_tree")
dat <- full_join(dat, treespp)

OGU <- read.table("~/Documents/bigBryogeography/round2/GeffertData/OGU.txt", sep="\t", header=T)
OGU <- OGU[,1:4]
OGU$Continent[is.na(OGU$Continent)] <- "Nam"
OGU <- OGU %>% filter(!is.na(Continent))
dat2 <- left_join(dat, OGU) %>% distinct()

#PD
m <- lm(PD_tree ~ PD_gbif, dat)
res_PD <- resid(m)
plotdat <- dat2 %>% filter(!is.na(PD_gbif)) %>% filter(!is.na(PD_tree))
plotdat <- cbind(plotdat, res_PD)
plotdat$fit_PD <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_PD, y = res_PD, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - PD")
ggplot(plotdat) + geom_point(aes(x= fit_PD, y = res_PD, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: PD (tree) ~ PD (gbif)") +
  annotate("text", x=17000, y=10000, label= "Corr = 0.7")

labelz <- c("Antarctica", 
            "Africa", 
            "Asia", 
            "SE Asia", 
            "Australia", 
            "Central America", 
            "Europe", 
            "North America",
            "South America")

valuez <- c("#FDBF6F", "#E31A1C", "#A6CEE3", "#33A02C", "#CAB2D6", "#B2DF8A", "#1F78B4", "#FF7F00", "#FB9A99")

ggplot(plotdat, aes(x= fit_PD, y = res_PD, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()


pdf(file = "~/Documents/bigBryogeography/round2/geffFigz/regressions/treespp_vs_gbif_PD.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(plotdat, aes(x= fit_PD, y = res_PD, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
dev.off()


#MPD
m <- lm(MPD_tree ~ MPD_gbif, dat)
summary(m)
res_MPD <- resid(m)
plotdat <- dat2 %>% filter(!is.na(MPD_gbif)) %>% filter(!is.na(MPD_tree))
plotdat <- cbind(plotdat, res_MPD)
plotdat$fit_MPD <-  fitted(m)
ggplot(plotdat) + geom_point(aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_hline(aes(yintercept=0)) + xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted - MPD")
ggplot(plotdat) + geom_point(aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("Residuals vs. Fitted: MPD (tree) ~ MPD (gbif)") +
  annotate("text", x=450, y=80, label= "Corr = 0.76")

ggplot(plotdat, aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()

pdf(file = "~/Documents/bigBryogeography/round2/geffFigz/regressions/treespp_vs_gbif_MPD.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches
ggplot(plotdat, aes(x= fit_MPD, y = res_MPD, color =Continent)) + geom_point(alpha=.7, size = 2.5) + geom_hline(aes(yintercept=0)) +
  xlab("Fitted Values") + ylab("Residuals") + ggtitle("")  +
  scale_color_manual(labels = labelz, values = valuez) +theme_classic()
dev.off()

#





