#combine the two datasets into one and clean it up#########
rm(list = ls())
library(dplyr)
library(countrycode)
library(openxlsx)
library(stringr)
library(CoordinateCleaner)
gbif <- read.xlsx("gbif_done.xlsx")
idigbio <- read.xlsx("idig_done.xlsx")

#these datasets have different column names. here, i standardize them for easy combining.
gbifDF <- data.frame(speciesFinal = gbif$speciesFinal,
                     sciNameFinal = gbif$sciNameFinal,
                     tropFamFinal = gbif$tropFamFinal,
                     lat = gbif$decimalLatitude,
                     lon = gbif$decimalLongitude,
                     coordinateUncertainty = gbif$coordinateUncertaintyInMeters,
                     issues = gbif$issues,
                     databaseID = gbif$key,
                     country = gbif$country,
                     database = paste("gbif"),
                     species_syn = gbif$species_syn,
                     sciName_syn = gbif$sciName_syn)
gbifDF <- gbifDF %>% distinct()


idigDF <- data.frame(speciesFinal = idigbio$speciesFinal,
                     sciNameFinal = idigbio$sciNameFinal,
                     tropFamFinal = idigbio$tropFamFinal,
                     lat = idigbio$geopoint.lat,
                     lon = idigbio$geopoint.lon,
                     coordinateUncertainty = idigbio$coordinateuncertainty,
                     issues = idigbio$flags,
                     databaseID = idigbio$uuid,
                     country = str_to_title(idigbio$country),
                     database = paste("idigbio"),
                     species_syn = idigbio$species_syn,
                     sciName_syn = idigbio$sciName_syn)
idigDF <- idigDF %>% distinct()

distDFbad <- rbind(gbifDF,idigDF)

#JANUARY FIX###########
syns <- read.xlsx("synonyms_done.xlsx")
fix <- syns %>% filter(ID_syn ==  35173398 | ID_syn == 35186284 | ID_syn == 35124560) %>% select(speciesFinal, sciNameFinal, species_syn, sciName_syn)
fixDF <- distDFbad %>% filter(species_syn =="Dicranum spadiceum var. subscabrifolium" | species_syn== "Orthotrichum fastigiatum" | species_syn=="Plagiothecium novae-seelandiae") %>% select(-speciesFinal, -sciNameFinal)
goodDF <- anti_join(distDFbad, fixDF)
nrow(goodDF) + nrow(fixDF) == nrow(distDFbad)

fixDF <- full_join(fixDF, fix) %>% filter(!is.na(database))
distDF <- rbind(goodDF, fixDF)

check <- distDF %>% filter(species_syn =="Dicranum spadiceum var. subscabrifolium" | species_syn== "Orthotrichum fastigiatum" | species_syn=="Plagiothecium novae-seelandiae") 

###TRENAME#######
#adding a column for matching this to tips of phylogenetic tree (which have to have "_" rather than spaces.)
distDF$treeNames <- str_replace_all(distDF$speciesFinal, " ", "_")

#manually fix.
#distDF$treeNames[distDF$treeNames=="Cryptodicranum_armitii_"] <- "Cryptodicranum_armitii"
#distDF$treeNames[distDF$treeNames=="Plagiomnium_novae-zealandiae_"] <- "Plagiomnium_novae-zealandiae"
distDF$treeNames[distDF$treeNames=="Plagiothecium_novae-seelandiae"] <- "Plagiothecium_novaeseelandiae"


#check 
library(ape)
tree <- read.tree("~/Documents/bigBryogeography/newtrees/bestTree_jan.tre")
namez <- data.frame(treeNames = tree[["tip.label"]], inTree = TRUE)
namezDF <- distDF %>% select(treeNames) %>% distinct()
check <- left_join(namez, namezDF)
what <- left_join(namezDF, namez) %>% filter(is.na(inTree))


write.xlsx(distDF, "Distributions.xlsx", as.table = F)
write.csv(distDF, "Distributions.csv", row.names = F)
##Clean#
rm(list = ls())
distDF <- read.xlsx("Distributions.xlsx")
library(CoordinateCleaner)
distDF_save <- distDF

specs <- distDF %>% select(speciesFinal) %>% distinct()

distDF <- left_join(specs, distDF_save)

#cc_dupl Identify Duplicated Records
dup1DF <- cc_dupl(distDF, lon = "lon", lat = "lat", species = "speciesFinal")

#cc_inst Identify Records in the Vicinity of Biodiversity Institutions
inst_DF <- cc_inst(dup1DF, lon = "lon", lat = "lat", species = "speciesFinal")

#cc_cap Identify Coordinates in Vicinity of Country Capitals.
capDF <- cc_cap(inst_DF, lon = "lon", lat = "lat", species = "speciesFinal", verify = T)

#cc_cen Identify Coordinates in Vicinity of Country and Province Centroids
cenDF <- cc_cen(capDF, lon = "lon", lat = "lat", species = "speciesFinal", verify = T)

#cc_equ Identify Records with Identical lat/lon
eqDF <- cc_equ(cenDF, lon = "lon", lat = "lat")

#cc_val Identify Invalid lat/lon Coordinates
valDF <- cc_val(eqDF, lon = "lon", lat = "lat")

#cc_zero Identify Zero Coordinates
#zeroDF <- cc_zero(valDF, lon = "lon", lat = "lat")

#cc_coun Identify Coordinates Outside their Reported Country
naDF <-  valDF %>% filter(is.na(countrycode))
gDF <- valDF %>% filter(!is.na(countrycode))
check <- rbind(naDF, gDF)
rm(check)
counDF <- cc_coun(gDF, lon = "lon", lat = "lat", iso3 = "countrycode")
counDF <- rbind(counDF, naDF)

#cc_gbif Identify Records Assigned to GBIF Headquarters
gbifDF <- cc_gbif(valDF, lon = "lon", lat = "lat")
gbifDF <- cc_gbif(counDF, lon = "lon", lat = "lat")


#cc_sea Identify Non-terrestrial Coordinates
library(rgdal)
ref <- readOGR('~/Documents/Shapefiles/buffers/buff.2.shp')
#seaDF_1 <- cc_sea(gbifDF, lon = "lon", lat = "lat", ref = ref)
seaDF <- cc_sea(gbifDF, lon = "lon", lat = "lat",  ref = ref)

sea <- cc_sea(gbifDF, lon = "lon", lat = "lat", value = "flagged",  ref = ref)
seaDF <- gbifDF
seaDF$flag_sea <- sea


str(seaDF)
counts <- data.frame(table(seaDF$speciesFinal)) %>% mutate(speciesFinal = Var1) %>% select(-Var1)
seaDF_fin <- full_join(seaDF, counts)
write.xlsx(seaDF_fin, "seaDF.xlsx", asTable = F)
write.csv(seaDF_fin, "seaDF.csv", row.names = F)

#make one only for mosses:
genera <- read.xlsx("/Users/ranunculus/Documents/bigBryogeography/genera.xlsx")
genz <- genera %>% select(genusFinal, Division.taxon.name, Subdivision.taxon.name)
seaDF$genusFinal <- word(seaDF$speciesFinal, 1, 1)
dist <- left_join(seaDF, genz)
check <- dist %>% filter(is.na(Division.taxon.name))
unique(dist$Division.taxon.name)
dist_moss <- dist %>% filter(Division.taxon.name == "Bryophyta")
write.xlsx(dist_moss, "mossDF.xlsx", asTable = F)
write.csv(dist_moss, "mossDF.csv", row.names=F)

