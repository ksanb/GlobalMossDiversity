#get data from GBIF and iDigBio
#get GBIF Data########
rm(list = ls())
library(rgbif)
spec <- read.xlsx("syns_done.xlsx") %>% select(!V1) %>% distinct()
spec_smol <- spec %>% select(TropID_Final, speciesFinal, sciNameFinal, tropFamFinal, ID_syn, species_syn, sciName_syn , fam_syn, taxonRank, syn_author, homonym, genusFinal, Division.taxon.name, Subdivision.taxon.name, Order.taxon.name, Family.taxon.name)
accs <- read.xlsx("accessionsFin.xlsx") %>% filter(fin==FALSE) %>% select(speciesFinal) %>% distinct()
sp <- left_join(accs, spec) #this is done to exclude certain species whose taxonomy is unresolved. I only had this information in the accessions file.

#this function (occ_search) downloads data from gbif based on a name.
gbif <- occ_search(scientificName=as.character(sp$species_syn[1]), hasCoordinate=T,hasGeospatialIssue=F, basisOfRecord="PRESERVED_SPECIMEN", limit=100000)
col <- c("speciesFinal", "species_syn", "key", "issues",  "scientificName", "decimalLatitude", "decimalLongitude", "issues", "datasetKey", "publishingOrgKey", "publishingCountry", "basisOfRecord", "taxonKey", "kingdomKey", "phylumKey", "classKey", "orderKey", "familyKey", "genusKey", "speciesKey", "acceptedTaxonKey", "acceptedScientificName", "kingdom", "phylum", "order", "family", "genus", "species", "genericName", "specificEpithet", "taxonRank", "taxonomicStatus", "elevation", "year", "month", "day", "eventDate", "class", "countryCode", "country", "identifier", "recordedBy", "catalogNumber", "institutionCode", "fieldNumber", "gbifID", "collectionCode", "occurrenceID", "identifiedBy", "name", "stateProvince", "locality", "continent", "nomenclaturalStatus",  "nomenclaturalCode", "datasetName", "taxonID", "ownerInstitutionCode", "collectionID", "higherGeography", "higherClassification", "coordinateUncertaintyInMeters", "county", "elevationAccuracy", "habitat", "verbatimEventDate", "endDayOfYear", "startDayOfYear", "verbatimElevation", "sp", "municipality", "identificationID", "coordinatePrecision", "occurrenceRemarks", "verbatimCoordinateSystem", "typeStatus", "islandGroup")
gbif <- gbif[["data"]] %>% as.data.frame() %>% select(any_of(col))
gbif$species_syn <- as.character(sp$species_syn[1])
gbif$speciesFinal <- paste(sp$speciesFinal[1])

#the following is a for loop that does the exact same thing as above for each species.
for(i in 2:nrow(sp)){
  query <- as.character(sp$species_syn[i])
  raw<-occ_search(scientificName=query, hasCoordinate=T, hasGeospatialIssue=F, basisOfRecord="PRESERVED_SPECIMEN", limit=100000)
  dat <- data.frame(raw[["data"]])
  print("gbif data done")
  if(nrow(dat)==0){}
  else{
    dat$species_syn <- paste(query)
    dat$speciesFinal <- paste(sp$speciesFinal[i])
    data <- dat %>% select(any_of(col))
    gbif<-bind_rows(gbif,data)
  }
  print(i)
}
gbif_save <- gbif #i'm just saving this in case I need to redo something.

gbif <- gbif_save
gbif$gbif_taxonRank <- gbif$taxonRank #I need to include taxon ranks for later steps.
#gbif$taxonRank[gbif$taxonRank == "SPECIES"] <- "sp."
#gbif$taxonRank[gbif$taxonRank == "FORM"] <- "f."
#gbif$taxonRank[gbif$taxonRank == "SUBSPECIES"] <- "subsp."
#gbif$taxonRank[gbif$taxonRank == "VARIETY"] <- "var."
gbif <- gbif %>% select(-species_syn, -speciesFinal, -taxonRank)
gbif$sciName_syn <- gbif$scientificName
good <- left_join(gbif, spec_smol) %>% filter(!is.na(speciesFinal)) #this thing combines the species list with the gbif data by SCIENTIFIC NAME.
bad <- left_join(gbif, spec_smol) %>% filter(is.na(speciesFinal)) %>% select(-TropID_Final, -speciesFinal, -sciNameFinal, -tropFamFinal, -ID_syn, - homonym, -species_syn, -sciName_syn, -fam_syn, -taxonRank, -syn_author, -genusFinal, -Division.taxon.name, -Subdivision.taxon.name, -Order.taxon.name, -Family.taxon.name)
#'bad' are where scientific names don't match between the gbif dataset and the species list. This could be because of spelling or because of homonyms. 

bad_sciNames <- bad %>% select(scientificName) %>% distinct() #pulling out only the scientific names that don't match
bad_sciNames$sciName_syn <- NA
#the following loop uses a string distance function to partially match names, using a threshold limit for matching.
for(i in 1:nrow(sp)){
  pattern <- sp$sciName_syn[i]
  if(is.na(pattern)){}
  else{
    for(j in 1:nrow(bad_sciNames)){
      string <- bad_sciNames$scientificName[j]
      dist <- stringdist(pattern, string)
      if(dist < 16){
        bad_sciNames$sciName_syn[j] <- pattern
      }
    }
    print(i)
  }
}


###THESE SCIENTIFIC NAMES DON'T MATCH. I'VE DONE THAT HERE MANUALLY###
bad_sciNames <- data.frame(scientificName = unique(bad$scientificName), sciName_syn = c("Camptothecium lutescens var. fallax (H. Philib. ex Schimp.) Breidl.",
                                                                                        "Homalothecium lutescens (Hedw.) H. Rob.",
                                                                                        "Homalothecium fulgescens (Mitt. ex Müll. Hal.) A. Jaeger",
                                                                                        "Homalothecium lutescens subsp. fulgescens (Mitt. ex Müll. Hal.) Heike Hofm.",
                                                                                        "Camptothecium fulgescens (Mitt. ex Müll. Hal.) Paris",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Homalothecium lutescens var. fallax (H. Philib. ex Schimp.) Düll",
                                                                                        "delete",
                                                                                        "Camptothecium fallax (H. Philib. ex Schimp.) Schiffn.",
                                                                                        "Campylopus capillatus Hook. f. & Wilson",
                                                                                        "Campylopus sparksii R. Br. bis",
                                                                                        "Dicranum pulvinatum (Hedw.) Sw. ex Lag., D. García & Clemente",
                                                                                        "Campylopus pyriformis (Schultz) Brid.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Campylopus woollsii (Müll. Hal.) Paris",
                                                                                        "Campylopus kunkelii E.B. Bartram",
                                                                                        "Campylopus novae-zealandiae E.B. Bartram & Dixon",
                                                                                        "Campylopus tijucae Broth.",
                                                                                        "Campylopus tijucae Broth.",
                                                                                        "Campylopus torfaceus var. fallaciosus Thér.",
                                                                                        "Diphyscium pilmaiquen (Crosby) Magombo",
                                                                                        "Muscoflorschuetzia pilmaiquen (Crosby) Crosby",
                                                                                        "Florschuetzia pilmaiquen Crosby",
                                                                                        "Bryobartramia novae-valesiae (Broth. ex G. Roth) I.G. Stone & G.A.M. Scott",
                                                                                        "Zygodon spathulifolius Besch.",
                                                                                        "delete",
                                                                                        "Bryomaltaea obtusifolia (Hook.) Goffinet",
                                                                                        "Orthotrichum crassifolium Hook. f. & Wilson",
                                                                                        "delete",
                                                                                        "Muelleriella crassifolia (Hook. f. & Wilson) Dusén",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Orthomitrium tuberculatum Lewinsky & Crosby",
                                                                                        "Gemmabryum apiculatum (Schwägr.) J.R. Spence & H.P. Ramsay",
                                                                                        "Brachymenium wattsii Broth.",
                                                                                        "Bryum apiculatum Schwägr.",
                                                                                        "Bryum nitens Hook. ex Harv.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Bryum nitens Hook. ex Harv.",
                                                                                        "Gemmabryum barnesii (J.B. Wood ex Schimp.) J.R. Spence",
                                                                                        "Bryum barnesii J.B. Wood ex Schimp.",
                                                                                        "Gemmabryum caespiticium (Hedw.) J.R. Spence",
                                                                                        "Bryum caespiticioides Müll. Hal.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Ptychostomum imbricatulum (Müll. Hal.) Holyoak & N. Pedersen",
                                                                                        "Gemmabryum coronatum (Schwägr.) J.R. Spence & H.P. Ramsay",
                                                                                        "Bryum coronatum Schwägr.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Bryum pachytheca Müll. Hal.",
                                                                                        "Bryum otahapaense R. Br. bis",
                                                                                        "delete",
                                                                                        "Gemmabryum radiculosum (Brid.) J.R. Spence & H.P. Ramsay",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Bryum balticum Nyholm & Hedenäs",
                                                                                        "Bryum bellii Müll. Hal.",
                                                                                        "delete",
                                                                                        "Bryum ovatocarpum R. Br. bis",
                                                                                        "Bryum ovatothecium R. Br. bis",
                                                                                        "Bryum petriei R. Br. bis",
                                                                                        "Bryum urbanskyi Broth.",
                                                                                        "Bryum waikariense R. Br. bis",
                                                                                        "Bryum microglobum Müll. Hal. & Kindb.",
                                                                                        "Bryum dichotomum Hedw.",
                                                                                        "Bryum barnesii J.B. Wood ex Schimp.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Bryum webbianum R. Br. bis",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Bryum gemmiferum R. Wilczek & Demaret",
                                                                                        "Gemmabryum gemmilucens (R. Wilczek & Demaret) J.R. Spence",
                                                                                        "Bryum gemmilucens R. Wilczek & Demaret",
                                                                                        "Gemmabryum klinggraeffii (Schimp.) J.R. Spence & H.P. Ramsay",
                                                                                        "Bryum klinggraeffii Schimp.",
                                                                                        "Gemmabryum radiculosum (Brid.) J.R. Spence & H.P. Ramsay",
                                                                                        "Bryum radiculosum Brid.",
                                                                                        "Bryum murale Wilson ex Hunt",
                                                                                        "Bryum subdecursivum Müll. Hal.",
                                                                                        "Bryum fendleri Müll. Hal.",
                                                                                        "delete",
                                                                                        "Gemmabryum ruderale (Crundw. & Nyholm) J.R. Spence",
                                                                                        "Bryum ruderale Crundw. & Nyholm",
                                                                                        "Gemmabryum subapiculatum (Hampe) J.R. Spence & H.P. Ramsay",
                                                                                        "Bryum microerythrocarpum Müll. Hal. & Kindb.",
                                                                                        "Bryum subapiculatum Hampe",
                                                                                        "Bryum microerythrocarpum Müll. Hal. & Kindb.",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "delete",
                                                                                        "Gemmabryum violaceum (Crundw. & Nyholm) J.R. Spence",
                                                                                        "Bryum violaceum Crundw. & Nyholm",
                                                                                        "Philonotis bartramioides (Griff.) D.G. Griffin & W.R. Buck",
                                                                                        "Bartramidula bartramioides (Griff.) Wijk & Margad.",
                                                                                        "Glossadelphus glossoides (Bosch & Sande Lac.) M. Fleisch.",
                                                                                        "Phyllodon lingulatus (Cardot) W.R. Buck",
                                                                                        "delete",
                                                                                        "Myurella acuminata Lindb. & Arnell",
                                                                                        "Pseudotrachypus wallichii (Brid.) W.R. Buck",
                                                                                        "Pseudobarbella attenuata (Thwaites & Mitt.) Nog.",
                                                                                        "Syrrhopodon banksii Müll. Hal.",
                                                                                        "Syrrhopodon bornensis (Hampe) A. Jaeger",
                                                                                        "Syrrhopodon involutus Schwägr.",
                                                                                        "delete",
                                                                                        "Tayloria cochabambae Müll. Hal.",
                                                                                        "Brachymitrion cochabambae (Müll. Hal.) A.K. Kop.",
                                                                                        "Brachymitrion immersum Goffinet",
                                                                                        "Tayloria jamesonii (Taylor) Müll. Hal.",
                                                                                        "Brachymitrion jamesonii Taylor",
                                                                                        "delete",
                                                                                        "Tayloria moritziana Müll. Hal.",
                                                                                        "Brachymitrion moritzianum (Müll. Hal.) A.K. Kop.",
                                                                                        "Tayloria moritziana var. carbonaria Müll. Hal.",
                                                                                        "Bryonoguchia molkenboeri (Sande Lac.) Z. Iwats. & Inoue",
                                                                                        "Tortula lanceola R.H. Zander",
                                                                                        "Pottia lanceolata (Hedw.) Müll. Hal.",
                                                                                        "Pottia lanceolata (Hedw.) Müll. Hal.",
                                                                                        "delete",
                                                                                        "Tortula modica R.H. Zander",
                                                                                        "Pottia truncata var. major (F. Weber & D. Mohr) Bruch & Schimp.",
                                                                                        "Tortula pallida (Lindb.) R.H. Zander",
                                                                                        "Tortula wilsonii (Hook.) R.H. Zander",
                                                                                        "Pottia wilsonii (Hook.) Bruch & Schimp.",
                                                                                        "delete",
                                                                                        "Pottia crinita Wilson ex Bruch & Schimp.",
                                                                                        "Ulota perichaetialis (Sainsbury) Goffinet",
                                                                                        "Leptodontiopsis fragilifolia Broth.",
                                                                                        "delete"))
#
bad_fixed <- full_join(bad, bad_sciNames) %>% filter(!sciName_syn == "delete") #now the old dataframe has the new sciNames
check <- bad_fixed %>% filter(is.na(sciName_syn))
bad_done <- left_join(bad_fixed, spec_smol) #again, combining the fixed gbif data with the species list.
check <- bad_done %>% filter(is.na(sciName_syn))

bif <- bind_rows(good, bad_done)
check <- bif %>% filter(is.na(speciesFinal))
bif <- bif %>% select(TropID_Final, speciesFinal, sciNameFinal, tropFamFinal, 
                      ID_syn, species_syn, sciName_syn, taxonRank, syn_author, 
                      homonym, key, issues, scientificName, decimalLatitude, decimalLongitude, 
                      datasetKey, publishingCountry, taxonKey, kingdomKey, phylumKey, 
                      classKey, orderKey, familyKey, genusKey, speciesKey, acceptedTaxonKey, 
                      acceptedScientificName, kingdom, phylum, order, family, genus, 
                      species, genericName, specificEpithet, elevation, year, month, day, 
                      countryCode, country, catalogNumber, institutionCode, gbifID, 
                      occurrenceID, nomenclaturalStatus, nomenclaturalCode, taxonID, 
                      coordinateUncertaintyInMeters, elevationAccuracy, verbatimEventDate, 
                      verbatimElevation, 
                      verbatimCoordinateSystem, gbif_taxonRank, fam_syn)
bif <- bif %>% mutate(key = as.numeric(key), gbifID = as.numeric(gbifID)) %>% distinct()

check <- bif %>% filter(is.na(speciesFinal)) #making sure nothing is wrong with the dataframe
check <- bif %>% filter(duplicated(key))


#THESE NAMES ARE IN MY GROUP AND I KNOW THE PROPER TAXONOMY
bif_save <- bif
specz <- read.xlsx("syns_done.xlsx")
#specz$bind <- paste(specz$ID_syn, specz$sciName_syn, sep = "")

#bif$bind <- paste(bif$ID_syn, bif$sciName_syn, sep = "")
bif$ID_syn[bif$species_syn == "Lewinskya speciosa var. speciosa"] <- 1000
bif$ID_syn[bif$species_syn == "Macrocoma perrottetii"] <- 2000
bif$ID_syn[bif$species_syn == "Mastopoma haidensis"] <- 3000
bif$ID_syn[bif$species_syn == "Plagiobryum algovicum"] <- 4000
bif$ID_syn[bif$species_syn == "Plagiothecium silvaticum var. rhynchostegioides"] <- 5000
bif$ID_syn[bif$species_syn == "Plagiothecium silvaticum var. latifolium"] <- 6000
bif$ID_syn[bif$species_syn == "Platygyrium leptohymenioides"] <- 7000
bif$ID_syn[bif$species_syn == "Sphagnum palenae"] <- 8000
bif$ID_syn[bif$species_syn == "Sphagnum sjorsii"] <- 9000
bif$ID_syn[bif$species_syn == "Trichosteleum cuspidatum"] <- 10000
bif$ID_syn[bif$species_syn == "Ulota membranacea "] <- 11000
bif$ID_syn[bif$species_syn == "Limbella pachylomata"] <- 12000
check <- bif %>% filter(is.na(ID_syn))

gbif <- bif %>% select(ID_syn,
                       key,
                       issues,
                       scientificName,
                       decimalLatitude,
                       decimalLongitude,
                       datasetKey,
                       publishingCountry,
                       taxonKey,
                       kingdomKey,
                       phylumKey,
                       classKey,
                       orderKey,
                       familyKey,
                       genusKey,
                       speciesKey,
                       acceptedTaxonKey,
                       acceptedScientificName,
                       kingdom,
                       phylum,
                       order,
                       family,
                       genus,
                       species,
                       genericName,
                       specificEpithet,
                       elevation,
                       year,
                       month,
                       day,
                       countryCode,
                       country,
                       catalogNumber,
                       institutionCode,
                       gbifID,
                       occurrenceID,
                       nomenclaturalStatus,
                       nomenclaturalCode,
                       taxonID,
                       coordinateUncertaintyInMeters,
                       elevationAccuracy,
                       verbatimEventDate,
                       verbatimElevation,
                       identificationID,
                       coordinatePrecision,
                       occurrenceRemarks,
                       verbatimCoordinateSystem,
                       islandGroup,
                       gbif_taxonRank)
done <- left_join(gbif, specz) %>% distinct() 
check <- done %>% filter(is.na(speciesFinal))
done <- done %>% filter(!is.na(speciesFinal))
done_save <- done


check <- bif %>% select(key, speciesFinal) %>% distinct()
check$dup <- duplicated(check$key)

#doing some more manual checks and cleaning...
wtf <- done %>% select(key, speciesFinal) %>% distinct()
wtf$dup <- duplicated(wtf$key)
df <- wtf %>% filter(dup==TRUE) %>% select(key) %>% distinct()
bad_ones <- left_join(df, done)
good_ones <- anti_join(done, bad_ones)

bad <- bad_ones %>% select(key,
                           issues,
                           scientificName,
                           decimalLatitude,
                           decimalLongitude,
                           datasetKey,
                           publishingCountry,
                           taxonKey,
                           kingdomKey,
                           phylumKey,
                           classKey,
                           orderKey,
                           familyKey,
                           genusKey,
                           speciesKey,
                           acceptedTaxonKey,
                           acceptedScientificName,
                           kingdom,
                           phylum,
                           order,
                           family,
                           genus,
                           species,
                           genericName,
                           specificEpithet,
                           elevation,
                           year,
                           month,
                           day,
                           countryCode,
                           country,
                           catalogNumber,
                           institutionCode,
                           gbifID,
                           occurrenceID,
                           nomenclaturalStatus,
                           nomenclaturalCode,
                           taxonID,
                           coordinateUncertaintyInMeters,
                           elevationAccuracy,
                           verbatimEventDate,
                           verbatimElevation,
                           identificationID,
                           coordinatePrecision,
                           occurrenceRemarks,
                           verbatimCoordinateSystem,
                           islandGroup,
                           gbif_taxonRank)
namez <- bad %>% select(scientificName) %>% distinct()
namez$sciName_syn <- c("Camptothecium lutescens var. fallax (H. Philib. ex Schimp.) Breidl.", 
                       "Campylopus pallidus Hook. f. & Wilson", 
                       "Bryum caespiticium Hedw.", 
                       "Gemmabryum dichotomum (Hedw.) J.R. Spence & H.P. Ramsay", 
                       "Syrrhopodon revolutus Dozy & Molk.",
                       "Syrrhopodon involutus Schwägr.")
bad2 <- full_join(bad, namez)
bad_done <- left_join(bad2, specz)
bad_done <- bad_done[-c(6711, 6713), ]

donedone <- bind_rows(good_ones, bad_done)
check <- donedone %>% select(speciesFinal, key) %>% distinct()
check$dup <- duplicated(check$key)
check <- donedone %>% filter(is.na(speciesFinal))

write.xlsx(donedone, "gbif_done.xlsx", asTable = F)


####GET DATA FROM IDIGBIO #############
#Sys.setenv(TROPICOS_KEY="") #NEED to input your tropicose API Key here.
rm(list = ls())
library(openxlsx)
library(dplyr)
library(stringdist)
library(tm)
library(stringr)
library(taxize)
library(ridigbio)
spec <- read.xlsx("syns_done.xlsx") %>% select(!V1) %>% distinct()
accs <- read.xlsx("accessionsFin.xlsx") %>% filter(fin==FALSE) %>% select(speciesFinal) %>% distinct()
#getting rid of species with bad taxonomy, this info is in the accessions file.
sp <- left_join(accs, spec)

species <- unique(sp$species_syn) 

#this loop grabs data from idigbio
idigbio <- NULL
for(i in 1:length(species)){
  print(i)
  query <- as.character(species[i])
  raw<-idig_search_records(rq=list("scientificname"=query), fields = "all")
  if(nrow(raw)=="0"){} else {
    raw$sp<- paste(query)
    data <- raw %>% filter(!is.na(geopoint.lon)) %>% filter(!is.na(geopoint.lat))
    idigbio<-rbind(idigbio,data)
  }
}
idigbio_save <- idigbio
cols <- c("uuid", "associatedsequences", "basisofrecord", "catalognumber", "class", "collectionid", "collectionname", "collector", "continent", "coordinateuncertainty", "country", "countrycode", "county", "datasetid", "datecollected",  "eventdate", "family", "fieldnumber", "flags", "genus", "geopoint.lon", "geopoint.lat", "group", "highertaxon", "highestbiostratigraphiczone", "individualcount", "infraspecificepithet", "institutioncode", "institutionid", "institutionname", "kingdom", "locality", "maxelevation", "mediarecords", "member", "minelevation", "municipality", "occurrenceid", "order", "phylum", "recordids", "recordnumber", "recordset", "scientificname", "specificepithet", "stateprovince", "taxonid", "taxonrank", "typestatus", "verbatimeventdate", "verbatimlocality", "waterbody", "sp")
idigbio <- idigbio_save %>% select(any_of(cols)) #grabbing only the columns i am interested in
idigbio$species_syn <- idigbio$sp
idigbio <- idigbio %>% distinct() #removing duplicates

#THESE ARE downloaded from idigbio: simply search "moss", "hornwort", "liverwort","bryophyta", "hepatophyta", and "musci"
#this is done because the idig_search function doesn't return scientific names.
idig_d_moss <- read.csv("occurrence_raw.csv")
idig_d_horn <- read.csv("occurrence_raw_horn.csv") #%>% mutate(dwc.coordinateUncertaintyInMeters = as.character(dwc.coordinateUncertaintyInMeters)) %>% mutate(dwc.day = as.character(dwc.day)) %>% mutate(dwc.eventID = as.character(dwc.eventID)) %>% mutate(dwc.month = as.character(dwc.month))
idig_d_liv <- read.csv("occurrence_raw_liver.csv") #%>% mutate(dwc.coordinateUncertaintyInMeters = as.character(dwc.coordinateUncertaintyInMeters)) 
idig_try <- read.csv("occurrence_raw_try.csv")
idig_liver2 <- read.csv("occurrence_raw_liver2.csv") 
idig_musci <- read.csv("occurrence_raw_musci.csv") 

#and combining because they are very large...
idig_direct <- rbind(idig_d_moss, idig_d_liv)
idig_direct <- rbind(idig_direct, idig_d_horn)
idig_direct <- rbind(idig_direct, idig_try)
idig_direct <- rbind(idig_direct, idig_liver2)
idig_direct <- rbind(idig_direct, idig_musci)


#grabbing the scientific names from the downloaded datasets for uuid's that match in my dataset
idigbio$scientificNameAuthorship <- idig_direct$dwc.scientificNameAuthorship[match(idigbio$uuid, idig_direct$coreid)]
idigbio$sciName <- paste(idigbio$sp, idigbio$scientificNameAuthorship, sep = " ")
head(idigbio)
probs <- idigbio %>% filter(scientificNameAuthorship == "")
dig <- anti_join(idigbio, probs)
check <- rbind(dig, probs)

##Fixing scientific names
namez <- dig %>% select(sciName, species_syn) %>% distinct()
#use a string distance function to partial match scientific names.
namez$sciName_syn <- NA
for(i in 1:nrow(sp)){
  pattern <- sp$sciName_syn[i]
  if(is.na(pattern)){}
  else{
    for(j in 1:nrow(namez)){
      string <- namez$sciName[j]
      dist <- stringdist(pattern, string)
      if(dist < 4){
        namez$sciName_syn[j] <- pattern
      }
    }
    print(i)
  }
}
namezsave <- namez

#Doing this mannually
namez$sciName_syn[grep(0, namez$sciName)] <- paste("DELETE")  #some of them have zeros... get rid of these.
namez$sciName_syn[grep("musgo", namez$sciName)] <- paste("DELETE")  #some of them have "musgo" instead of an author... get rid of these.
namez$sciName_syn[namez$sciName == "Camptothecium lutescens (Hedw.) H. Rob."] <-paste("Homalothecium lutescens (Hedw.) H. Rob.")
namez$sciName_syn[namez$sciName == "Campylopus pyriformis (F.W.Schultz) Brid."] <- paste("Campylopus pyriformis (Schultz) Brid.")
namez$sciName_syn[namez$sciName == "Zygodon araucariae C. M."] <- paste("Zygodon araucariae Müll. Hal.")
namez$sciName_syn[namez$sciName == "Gemmabryum apiculatum (Schwägr.) Spence & Ransay."] <- paste("Gemmabryum apiculatum (Schwägr.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Gemmabryum apiculatum (Schwägr.) Spence & H.P.Ramsay"] <- paste("Gemmabryum apiculatum (Schwägr.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Bryum barnesii J.B. Wood"] <- paste("Bryum barnesii J.B. Wood ex Schimp.")
namez$sciName_syn[namez$sciName == "Gemmabryum caespiticium J.R. Spence"] <- paste("Gemmabryum caespiticium (Hedw.) J.R. Spence")
namez$sciName_syn[namez$sciName == "Gemmabryum caespiticium (Hedw.) Spence"] <- paste("Gemmabryum caespiticium (Hedw.) J.R. Spence")
namez$sciName_syn[namez$sciName == "Gemmabryum dichotomum J.R. Spence & H.P. Ramsay"] <- paste("Gemmabryum dichotomum (Hedw.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Gemmabryum klinggraeffii J.R. Spence & H.P. Ramsay"] <- paste("Gemmabryum klinggraeffii (Schimp.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Gemmabryum radiculosum J.R. Spence & H.P. Ramsay"] <- paste("Gemmabryum radiculosum (Brid.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Gemmabryum radiculosum (Brid.) spence & ramsay"] <- paste("Gemmabryum radiculosum (Brid.) J.R. Spence & H.P. Ramsay")
namez$sciName_syn[namez$sciName == "Bryum radiculosum Brid., 1817"] <- paste("Bryum radiculosum Brid.")
namez$sciName_syn[namez$sciName == "Bryum microerythrocarpum C. Müll. & Kindb. ex Mac. & Kindb."] <- paste("Bryum microerythrocarpum Müll. Hal. & Kindb.")
namez$sciName_syn[namez$sciName == "Bryum microerythrocarpum C. Müll. & Kindb. in Mac. & Kindb."] <- paste("Bryum microerythrocarpum Müll. Hal. & Kindb.")
namez$sciName_syn[namez$sciName == "Pottia intermedia R.H. Zander"] <- paste("Pottia intermedia (Turner) Fürnr.")
namez$sciName_syn[namez$sciName == "Homalothecium lutescens var. fallax (H. Philib.) Hedenäs & L. Söderstr."] <- paste("Homalothecium lutescens var. fallax (H. Philib. ex Schimp.) Düll")
namez$sciName_syn[namez$sciName == "Homalothecium lutescens var. fallax (H.Philib.) Hedenäs & L.Söderstr."] <- paste("Homalothecium lutescens var. fallax (H. Philib. ex Schimp.) Düll")




namez <- namez %>% select(-species_syn)
dig2 <- full_join(dig, namez) %>% filter(!sciName_syn == "DELETE") %>% select(-species_syn)

digFin <- left_join(dig2, sp) %>% distinct()
check <- digFin %>% filter(duplicated(uuid))

check <- digFin %>% filter(is.na(TropID_Final))


digFin <- digFin %>% select(-homo_sci1, -homo1_ID, -homo1_auth, -homo_sci2, -homo2_ID, -homo2_auth, -homo_sci3, -homo3_ID, -homo3_auth, -homo_sci4, -homo4_ID, -homo4_auth, -homo_sci5, -homo5_ID, -homo5_auth, -homo_sci6, -homo6_ID, -homo6_auth)
done2 <- left_join(digFin, spec)
check <- done2 %>% filter(is.na(TropID_Final))

donedone <- rbind(done2, digFin) %>% distinct()
check <- done2 %>% filter(is.na(speciesFinal))

check <- donedone %>% filter(duplicated(uuid))

write.xlsx(donedone, "idig_done.xlsx", asTable = F)


#now combine with original species list.
rm(list = ls())
idig <- read.xlsx("idig_done.xlsx") 
specz <- read.xlsx("syns_done.xlsx")
idigbio <- idig %>% select(sciName_syn, 
                           uuid, 
                           associatedsequences, 
                           basisofrecord, 
                           catalognumber, 
                           class, 
                           collectionid, 
                           collectionname, 
                           collector, 
                           continent, 
                           coordinateuncertainty, 
                           country, 
                           countrycode, 
                           county, 
                           datasetid, 
                           datecollected, 
                           eventdate, 
                           family, 
                           fieldnumber, 
                           flags, 
                           genus, 
                           geopoint.lon, 
                           geopoint.lat, 
                           group, 
                           highertaxon, 
                           highestbiostratigraphiczone, 
                           individualcount, 
                           infraspecificepithet, 
                           institutioncode, 
                           institutionid, 
                           institutionname, 
                           kingdom, 
                           locality, 
                           maxelevation, 
                           mediarecords, 
                           member, 
                           minelevation, 
                           municipality, 
                           occurrenceid, 
                           order, 
                           phylum, 
                           recordids, 
                           recordnumber, 
                           recordset, 
                           scientificname, 
                           specificepithet, 
                           stateprovince, 
                           taxonid, 
                           taxonrank, 
                           typestatus, 
                           verbatimeventdate, 
                           verbatimlocality, 
                           waterbody, 
                           sp, 
                           scientificNameAuthorship, 
                           sciName)
donedone <- left_join(idigbio, specz) %>% distinct() 
check <- donedone %>% filter(is.na(speciesFinal))
donedone <- donedone %>% filter(!is.na(speciesFinal))

write.xlsx(donedone, "idig_done.xlsx", asTable = F)



























