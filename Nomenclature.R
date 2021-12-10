#species
rm(list = ls())
setwd("/Users/ranunculus/Documents/bigBryogeography")
species <- read.csv("species.csv")
library(taxize)
library(dplyr)
species <- species %>% mutate(sp = as.character(sp))
query <- species$sp
g <-get_gbifid_(query)
ids <- lapply(g, function(z) z[1,])
df <- bind_rows(ids, .id = "species.sp")
df$V1 <- 1:nrow(df)
write.csv(df, "species.csv", row.names=F)

#NAMES
rm(list = ls())
library(dplyr)
library(taxize)
library(rentrez)
#Sys.setenv(TROPICOS_KEY=)
#Sys.setenv(ENTREZ_KEY=)
species <- read.csv("species.csv")
#species <-  species %>% slice(1:15)

df <- species %>% mutate(species = as.character(species))
df$TropID= NA
df$tropSciName= NA 
df$sciNameWithAuthors= NA
df$tropFam= NA
df$rank= NA
df$nomenclaturestatusid= NA
df$nomenclaturestatusname= NA
df$symbol= NA
df$author=NA
df <- df %>% mutate(author = as.character(author))
#df$NCBIfam=NA
#df$NCBIorder=NA
df$gbifID=NA
sapply(df, class)

#df <- data.frame(species=species$species, TropID= numeric(), tropSciName= character() , sciNameWithAuthors= character(), tropFam= character(), rank= character(), nomenclaturestatusid= numeric(), nomenclaturestatusname= character(), symbol= character() , author=character(), NCBIfam=character(), NCBIorder=character, gbifID=numeric())


for(i in 1:nrow(df)){
  query <- df$species[i]
  tropID <- tp_search(sci=query, type="exact")
  if(tropID[1,1]=="No names were found"){
    df[i,3]<-paste(tropID[1,1])
  } else {
    df[i,2]<-paste(tropID[1,1])
    df[i,3]<-paste(tropID[1,2])
    df[i,4]<-paste(tropID[1,3])
    df[i,5]<-paste(tropID[1,4])
    df[i,6]<-paste(tropID[1,5])
    if(tropID[1,6]=="No opinion"){
      df[i,7]<- paste(tropID[1,6])
      df[i,10] <- paste(tropID[1,7])
    } else {
      df[i,7]<- paste(tropID[1,6])
      df[i,8]<- paste(tropID[1,7])
      df[i,9]<- paste(tropID[1,8])
      df[i,10]<- paste(tropID[1,9])
    }}
  print(paste("values for species",i,"done",sep=" "))
}


#species <-  species %>% slice(1:13)
#list2 <- list()
list3 <- list()
list <- split( species , cut(1:nrow(species), 3936) )
#list <- split( species , cut(1:nrow(species), 13) )
data <- NULL
for(l in  1:length(list)){
  df.tmp <- list[[l]]   %>% mutate(species = as.character(species))
  for(i in 1:nrow(df.tmp)){
    query <- df.tmp$species[i]
    NCBI <- tax_name(query = query, get = c("family", "order"), db = "ncbi") 
    list3[[l]] <- NCBI
    print(paste("df.tmp",i,"done",sep=" "))
  }
  print(paste("list component",l,"done",sep=" "))
}

big_ncbi = do.call(rbind, list3)
big_ncbi <- big_ncbi %>% mutate(species = as.character(query)) %>% mutate(NCBIfam = as.character(family)) %>% mutate(NCBIorder = as.character(order)) %>% select(species, NCBIfam, NCBIorder) 

big_ncbi[2351,]
df[2351,]

df <- df %>% select(!NCBIfam) %>% select(!NCBIorder)
big_ncbi <- big_ncbi %>% select(!species)
df <- cbind(df, big_ncbi)


query <- df$species
gbif <-get_gbifid(query, ask=F)
gbif_ <- as.data.frame(gbif)  %>% mutate(gbifID = as.numeric(ids)) %>% select(gbifID)

df <- df %>% select(!gbifID)
df <- cbind(df, gbif_) 

write.csv(df, "species_ids.csv", row.names = F)


df <- read.csv("species_ids.csv") %>% mutate(species=as.character(species))
noNames <- df %>% filter(tropSciName== "No names were found")
df_no <- anti_join(df, noNames)
###For misspellings and such###
subsp <-noNames[ grep("subsp.", noNames$species),]
var <-noNames[ grep("var.", noNames$species),]
noNames <-noNames[- grep("subsp.", noNames$species),]
noNames <-noNames[- grep("var.", noNames$species),]
noNames$Genus <- word(noNames$species, 1)
noNames$specificEpithet <- word(noNames$species, 2)
for(i in 1:nrow(noNames)){
  tmp <- noNames$specificEpithet[i]
  if(tmp =="cf."){
    noNames$specificEpithet[i] <- paste(word(noNames$species[i], 3))
  }
}
noNames$sp <-  paste(noNames$Genus, substr(noNames$specificEpithet, start=1, stop=4))

for(i in 1:nrow(noNames)){
  query <- noNames$sp[i]
  tropID <- tp_search(sci=query)
  if(tropID[1,1]=="No names were found"){
    noNames[i,3]<-paste(tropID[1,1])
  } else {
    noNames[i,2]<-paste(tropID[1,1])
    noNames[i,3]<-paste(tropID[1,2])
    noNames[i,4]<-paste(tropID[1,3])
    noNames[i,5]<-paste(tropID[1,4])
    noNames[i,6]<-paste(tropID[1,5])
    if(tropID[1,6]=="No opinion"){
      noNames[i,7]<- paste(tropID[1,6])
      noNames[i,10] <- paste(tropID[1,7])
    } else {
      noNames[i,7]<- paste(tropID[1,6])
      noNames[i,8]<- paste(tropID[1,7])
      noNames[i,9]<- paste(tropID[1,8])
      noNames[i,10]<- paste(tropID[1,9])
    }}
  print(paste("values for species",i,"done",sep=" "))
}
###
manually_check <- noNames %>% filter(tropSciName== "No names were found") %>% select(!Genus)  %>% select(!specificEpithet)  %>% select(!sp) 

no_check <- anti_join(noNames, manually_check) %>% select(!Genus) %>% select(!specificEpithet)  %>% select(!sp)

manually_checkFin <- rbind(manually_check, var, subsp) %>% distinct()

done <- rbind(df_no, no_check) %>% mutate(TropID=as.character(TropID))

write.csv(done, "species_ids_done.csv", row.names = F)
write.csv(manually_checkFin, "species_ids_man.csv", row.names = F)

man <- read.csv("species_ids_man.csv")
man_tropID <- man %>% filter(!is.na(TropID)) %>% mutate(tropSciName=as.character(tropSciName)) %>% mutate(sciNameWithAuthors=as.character(sciNameWithAuthors))

for(i in 1:nrow(man_tropID)){
  query <- man_tropID$TropID[i]
  tropID <- tp_search(nameid=query)
  man_tropID[i,3]<-paste(tropID[1,2])
  man_tropID[i,4]<-paste(tropID[1,3])
  man_tropID[i,5]<-paste(tropID[1,4])
  man_tropID[i,6]<-paste(tropID[1,5])
  if(tropID[1,6]=="No opinion"){
    man_tropID[i,7]<- paste(tropID[1,6])
    man_tropID[i,10] <- paste(tropID[1,7])
  } else {
    man_tropID[i,7]<- paste(tropID[1,6])
    man_tropID[i,8]<- paste(tropID[1,7])
    man_tropID[i,9]<- paste(tropID[1,8])
    man_tropID[i,10]<- paste(tropID[1,9])
  }
  print(paste("values for species",i,"done",sep=" "))
}



man_name <- man %>% filter(is.na(TropID))  %>% mutate(tropSciName=as.character(tropSciName)) %>% mutate(sciNameWithAuthors=as.character(sciNameWithAuthors))
for(i in 1:nrow(man_name)){
  query <- man_name$tropSciName[i]
  tropID <- tp_search(sci=query)
  if(tropID[1,1]=="No names were found"){
    man_name[i,3]<-paste(tropID[1,1])
  } else {
    man_name[i,2]<-paste(tropID[1,1])
    man_name[i,5]<-paste(tropID[1,4])
    man_name[i,6]<-paste(tropID[1,5])
    if(tropID[1,6]=="No opinion"){
      man_name[i,7]<- paste(tropID[1,6])
      man_name[i,10] <- paste(tropID[1,7])
    } else {
      man_name[i,7]<- paste(tropID[1,6])
      man_name[i,8]<- paste(tropID[1,7])
      man_name[i,9]<- paste(tropID[1,8])
      man_name[i,10]<- paste(tropID[1,9])
    }}
  print(paste("values for species",i,"done",sep=" "))
}

man_done <- rbind(man_name, man_tropID)
done_done <- rbind(done, man_done)
write.csv(done_done, "species_ids_done.csv", row.names = F)
done_done$Genus <- word(done_done$species, 1)
nas <- done_done %>% filter(is.na(NCBIorder))
query <- nas %>% select(Genus) %>% as.matrix()
NCBI <- tax_name(query = query, get = c("family", "order"), db = "ncbi") 
n <- NCBI %>% mutate(NCBIfam=as.character(family)) %>% mutate(NCBIorder=as.character(order)) %>% select(NCBIfam, NCBIorder)
nas <- nas %>% select(!NCBIorder) %>% select(!NCBIfam)
nas <- cbind(nas, n)
not_nas <- done_done %>% filter(!is.na(NCBIorder))
done_done2 <- rbind(not_nas, nas)

query <- done_done2 %>% select(NCBIorder) %>%distinct() %>% as.matrix()
NCBI <- tax_name(query = query, get = c("class"), db = "ncbi") 
n <- NCBI %>% mutate(NCBIorder=as.character(query)) %>% mutate(NCBIclass=as.character(class)) %>% select(NCBIorder, NCBIclass)
done_done3 <- full_join(done_done2, n)

classification <- read.csv("classification.csv")
done_done4 <- full_join(done_done3, classification)



###checking some things
done_done4 <- done_done4  %>% mutate(tropSciName=as.character(tropSciName)) %>% mutate(sciNameWithAuthors=as.character(sciNameWithAuthors))
for(i in 1:nrow(done_done4)){
  if(is.na(done_done4$tropSciName[i])){
    query <- done_done4$TropID[i]
    tropID <- tp_search(nameid=query)
    done_done4[i,3]<-paste(tropID[1,2])
    done_done4[i,4]<-paste(tropID[1,3])
    done_done4[i,5]<-paste(tropID[1,4])
    done_done4[i,6]<-paste(tropID[1,5])
    if(tropID[1,6]=="No opinion"){
      done_done4[i,7]<- paste(tropID[1,6])
      done_done4[i,10] <- paste(tropID[1,7])
    } else {
      done_done4[i,7]<- paste(tropID[1,6])
      done_done4[i,8]<- paste(tropID[1,7])
      done_done4[i,9]<- paste(tropID[1,8])
      done_done4[i,10]<- paste(tropID[1,9])
    }
    print(paste("values for species",i,"done",sep=" "))
  }
}



write.csv(done_done4, "species_ids_cleaned.csv", row.names = F)


####a few more things.... #####
species <- read.csv("species_ids_cleaned.csv")

nom <- species %>% select(nomenclaturestatusid) %>% distinct()
one <- species %>% filter(nomenclaturestatusid==1)
two <- species %>% filter(nomenclaturestatusid==2)
three <- species %>% filter(nomenclaturestatusid==3)

species$tropAccName <- NA
species$tropAccID <- NA
species$sciAccName <- NA
species$tropAccFam <- NA

for(i in 1:nrow(species)){
  if(!is.na(species$TropID[i])){
    query <- species$TropID[i]
    trop <- tp_accnames(id=query)
    t <- as.data.frame(trop[[1]])
    if(t[1,1]=="No accepted names were found"){
      species[i,17] <- paste(t[1,1])
    } else {
      t2 <- as.data.frame(trop[[2]])
      species[i,17]<-paste(t2[1,2])
      species[i,18]<-paste(t2[1,1])
      species[i,19]<-paste(t2[1,3])
      species[i,20]<-paste(t2[1,4])
    }
    print(paste("values for species",i,"done",sep=" "))
  }
}

twothree <- species %>% filter(nomenclaturestatusid==3 | nomenclaturestatusid==2) %>% filter(is.na(tropAccID))
write.csv(twothree, "species_nomenclatureMAN.csv", row.names=F)
not <- anti_join(species, twothree)
twothree <- read.csv("species_nomenclatureMAN.csv")
for(i in 1:nrow(twothree)){
  query <- twothree$TropID[i]
  tropID <- tp_search(nameid=query)
  twothree[i,3]<-paste(tropID[1,2])
  twothree[i,4]<-paste(tropID[1,3])
  twothree[i,5]<-paste(tropID[1,4])
  twothree[i,6]<-paste(tropID[1,5])
  if(tropID[1,6]=="No opinion"){
    twothree[i,7]<- paste(tropID[1,6])
    twothree[i,10] <- paste(tropID[1,7])
  } else {
    twothree[i,7]<- paste(tropID[1,6])
    twothree[i,8]<- paste(tropID[1,7])
    twothree[i,9]<- paste(tropID[1,8])
    twothree[i,10]<- paste(tropID[1,9])
  }
  print(paste("values for species",i,"done",sep=" "))
}

twothree <- twothree %>% mutate(sciAccName=as.character(sciAccName)) %>% mutate(tropAccName=as.character(tropAccName)) %>% mutate(tropAccFam=as.character(tropAccFam))
for(i in 1:nrow(twothree)){
  if(is.na(twothree$tropAccID[i])){
    query <- twothree$TropID[i]
    trop <- tp_accnames(id=query)
    t <- as.data.frame(trop[[1]])
    if(t[1,1]=="No accepted names were found"){
      twothree[i,17] <- paste(t[1,1])
    } else {
      t2 <- as.data.frame(trop[[2]])
      twothree[i,17]<-paste(t2[1,2])
      twothree[i,18]<-paste(t2[1,1])
      twothree[i,19]<-paste(t2[1,3])
      twothree[i,20]<-paste(t2[1,4])
    }
    print(paste("values for species",i,"done",sep=" "))
  }
}

hmm <- rbind(twothree, not)
write.csv(hmm, "species_ids_cleaned_bad.csv", row.names = F)


#Fixing mistakes that were made...#####
rm(list = ls())
library(dplyr)
library(taxize)
library(rentrez)
#Sys.setenv(TROPICOS_KEY="")
#Sys.setenv(ENTREZ_KEY="")
species <- read.csv("species_ids_cleaned_bad.csv")
species$V1 <- 1:nrow(species)

withAcc <- species  %>% filter(!is.na(tropAccID)) %>% select(species, TropID, tropAccName, tropAccID)
withAcc$nomenclaturestatus <- NA
withAcc$nomenclaturestatusACC <- NA

for(i in 1:nrow(withAcc)){
  query <- withAcc$TropID[i]
  tropID <- tp_search(nameid =query, type="exact")
  if(tropID[1,1]=="No names were found"){
    withAcc[i,5]<- paste(tropID[1,1])
  } else {
    if(tropID[1,6]=="No opinion"){
      withAcc[i,5]<- paste(tropID[1,6])
    } else {
      withAcc[i,5]<- paste(tropID[1,6])
    }}
  query2 <- withAcc$tropAccID[i]
  tropID2 <- tp_search(nameid =query2, type="exact")
  if(tropID[1,1]=="No names were found"){
    withAcc[i,6]<- paste(tropID2[1,1])
  } else {
    if(tropID2[1,6]=="No opinion"){
      withAcc[i,6]<- paste(tropID2[1,6])
    } else {
      withAcc[i,6]<- paste(tropID2[1,6])
    }}
  print(paste("values for species",i,"done",sep=" "))
}


#okay...
withAcc$determination <- NA
#test <- df_acc_checks %>% slice(1:20)
for(i in 1:nrow(withAcc)){
  if(withAcc[i,5]=="No opinion"	& withAcc[i,6]==1){
    withAcc[i,7]<- "acc"
  }
  if(withAcc[i,5]>1	& withAcc[i,6]==1){
    withAcc[i,7]<- "acc"
  }
  if(withAcc[i,5]>1	& withAcc[i,6]=="No opinion"){
    withAcc[i,7]<- "acc"
  }
  if(withAcc[i,5]==1	& withAcc[i,6]=="No opinion"){
    withAcc[i,7]<- "OG"
  }
  if(withAcc[i,5]==1	& withAcc[i,6]>1){
    withAcc[i,7]<- "OG"
  }
  if(withAcc[i,5]=="No opinion"	& withAcc[i,6]>1){
    withAcc[i,7]<- "OG"
  }
  print(paste("species",i,"done",sep=" "))
  if(withAcc[i,5]==withAcc[i,6]){
    withAcc[i,7]<- "equal"
  }
}

species_noacc <- species  %>% filter(is.na(tropAccID)) %>% mutate(nomenclaturestatus = as.integer(nomenclaturestatusid)) %>% select(!nomenclaturestatusid)
species_noacc$determination <- "OG"
species_noacc$nomenclaturestatusACC <- "NA"


species_acc <- species  %>% filter(!is.na(tropAccID)) %>% select(!nomenclaturestatusid)
acc <- withAcc %>% filter(determination == "acc")
acc_sp <- left_join(acc, species_acc)
#acc_sp[, 5]== acc_sp[, 12]
OG <- withAcc %>% filter(determination == "OG")
OG_sp <- left_join(OG, species_acc)
equal <- withAcc %>% filter(determination == "equal")
eq_sp <- left_join(equal, species_acc)
#check <- rbind(OG, acc, equal)
#check <- rbind(OG_sp, acc_sp, eq_sp)
df_sp <- rbind(OG_sp, acc_sp, eq_sp, species_noacc)
#write.csv(df_sp, "species_ids_cleaned.csv", row.names = F)
df_sp2 <- df_sp
df_sp2$TropID_Final <- NA #24
df_sp2$speciesFinal <- NA #25
df_sp2$sciNameFinal <- NA #26
df_sp2$tropFamFinal <- NA #27

for(i in 1:nrow(df_sp2)){
  if(df_sp2$determination[i]=="OG"){
    df_sp2[i,24]<- paste(df_sp2[i, 2]) #ID
    df_sp2[i,25]<- paste(df_sp2[i, 8]) #name
    df_sp2[i,26]<- paste(df_sp2[i, 9]) #nameauthors
    df_sp2[i,27]<- paste(df_sp2[i, 10]) #fam
  }
  if(df_sp2$determination[i]=="acc"){
    df_sp2[i,24]<- paste(df_sp2[i, 4]) #ID
    df_sp2[i,25]<- paste(df_sp2[i, 3]) #name
    df_sp2[i,26]<- paste(df_sp2[i, 21]) #nameauthors
    df_sp2[i,27]<- paste(df_sp2[i, 22]) #fam
  }
  if(df_sp2$determination[i]=="equal"){
    df_sp2[i,24]<- paste(df_sp2[i, 4]) #ID
    df_sp2[i,25]<- paste(df_sp2[i, 3]) #name
    df_sp2[i,26]<- paste(df_sp2[i, 21]) #nameauthors
    df_sp2[i,27]<- paste(df_sp2[i, 22]) #fam
  }
  print(paste("species",i,"done",sep=" "))
}


write.csv(df_sp2, "species_ids_cleaned.csv", row.names = F)


############################
species <- read.csv("species_ids_cleaned.csv")

accs <- species %>% select(V1, tropAccID) %>% filter(!is.na(tropAccID))
#accs$Syn <- NA
#accs$SynID <- NA
#accs$SynSci <- NA
#accs$SynFam <- NA
#accs <-accs %>% mutate(Syn=as.character(Syn)) %>% mutate(SynSci=as.character(SynSci)) %>% mutate(SynFam=as.character(SynFam)) %>% mutate(SynID=as.integer(SynID))


synsofAcc <- NULL
for(i in 1:nrow(accs)){
  query <- accs[i,2]
  syns <- synonyms(query, db="tropicos")
  synsdf <- syns[[1]]
  synsdf$ID <- query
  synsofAcc <- rbind(synsdf, synsofAcc)
  print(paste("values for species",i,"done",sep=" "))
}
sapply(accs, class)
a <- data.frame(tropAccID=synsofAcc$ID, 
                SynID=synsofAcc$nameid, 
                Syn=synsofAcc$scientificname, 
                SynSci = synsofAcc$scientificnamewithauthors, 
                SynFam=synsofAcc$family)
a <-a %>% mutate(Syn=as.character(Syn)) %>% mutate(SynSci=as.character(SynSci)) %>% mutate(SynFam=as.character(SynFam)) %>% mutate(SynID=as.integer(SynID))

sapply(a, class)

acc_join <- full_join(accs, a) %>% distinct()



ogs <- species %>% select(V1, TropID) %>% filter(!is.na(TropID))
synsofOG <- NULL
for(i in 1:nrow(ogs)){
  query <- ogs[i,2]
  syns <- synonyms(query, db="tropicos")
  synsdf <- syns[[1]]
  synsdf$ID <- query
  synsofOG <- rbind(synsdf, synsofOG)
  print(paste("values for species",i,"done",sep=" "))
}
sapply(ogs, class)
o<- data.frame(TropID=synsofOG$ID, 
               SynID=synsofOG$nameid, 
               Syn=synsofOG$scientificname, 
               SynSci = synsofOG$scientificnamewithauthors, 
               SynFam=synsofOG$family)
o <-o %>% mutate(Syn=as.character(Syn)) %>% mutate(SynSci=as.character(SynSci)) %>% mutate(SynFam=as.character(SynFam)) %>% mutate(SynID=as.integer(SynID))

sapply(o, class)

og_join <- full_join(ogs, o) %>% distinct()

df1 <- data.frame(V1 = species$V1, 
                  TropID_Final = species$TropID_Final, 
                  speciesFinal = species$speciesFinal, 
                  sciNameFinal = species$sciNameFinal, 
                  tropFamFinal = species$tropFamFinal,
                  ID_OG = species$TropID, 
                  species_OG = species$tropSciName, 
                  sciName_OG = species$sciNameWithAuthors)
#ID_syn = species$tropSynIDofAcc, 
#species_syn = species$tropSynofAcc, 
#sciName_syn = species$synSciNameofAcc)
sapply(df1, class)

try1 <- og_join  %>% mutate(TropID_Final = as.integer(TropID)) %>% select(V1, SynID, Syn, SynSci, TropID_Final) 
df_OG <- left_join(df1, try1)
try2 <- acc_join %>% mutate(TropID_Final = as.integer(tropAccID)) %>% select(V1, SynID, Syn, SynSci, TropID_Final) 
df_acc <- left_join(df1, try2)

syns_df_all<- full_join(og_join, acc_join)
syns_df_hmm <-  rbind(df_acc, df_OG) %>% distinct()

for(i in 1:nrow(syns_df)){
  if(is.na(syns_df$SynID[i])){
    syns_df[i,9] <- paste(syns_df[i,2]) #ID
    syns_df[i,10] <- paste(syns_df[i,3]) #name
    syns_df[i,11] <- paste(syns_df[i,4]) #sci
  }
  print(paste("species",i,"done",sep=" "))
}

syns_df_fin <- syns_df %>% distinct()
write.csv(syns_df, "synonyms.csv", row.names = F)

####make a complete species list:
l1 <- data.frame(V1=species$V1, 
                 TropID=species$TropID_Final, 
                 sciName=species$sciNameFinal, 
                 species=species$speciesFinal)
l2 <- data.frame(V1=syns_df$V1, 
                 TropID=syns_df$ID_OG, 
                 sciName=syns_df$sciName_OG, 
                 species=syns_df$species_OG)
l3 <- data.frame(V1=syns_df$V1, 
                 TropID=syns_df$SynID, 
                 sciName=syns_df$SynSci, 
                 species=syns_df$Syn)

spList <- rbind(l1,l2,l3) %>% distinct()

write.csv(spList, "species_list_bad.csv", row.names=F)


### FIX NO NAMES
FIXSP <- read.csv("species_list_bad.csv")
fix1 <- FIXSP %>% filter(species=="No names were found")
fix2 <- anti_join(FIXSP, fix1)
FIXIDS <- read.csv("species_ids_cleaned.csv") %>% select(species, V1)
fix1$sciName <- fix1$species
fix1 <- fix1 %>% select(!species)
fixed <- left_join(fix1, FIXIDS, by="V1")

spFixed <- rbind(fix2, fixed)
spFixed <- spFixed %>% distinct()
write.csv(spFixed, "species_list.csv", row.names=F)



#SYNS ######
#start with different spellings
species <- read.csv("species_ids_cleaned.csv") %>% select(!speciesFinal) %>% select(!sciNameFinal) %>% select(!tropFamFinal) %>% select(!TropID_Final)

syns <- data.frame(V1 = species$V1, 
                   TropID_Final = species$tropID_NCBI, 
                   speciesFinal = species$NCBIspecies, 
                   sciNameFinal = species$sciName_NCBI, 
                   tropFamFinal = species$NCBIfam,
                   ID_OG = species$TropID, 
                   species_OG = species$tropSciName, 
                   sciName_OG = species$sciNameWithAuthors)
syns$ID_syn <- NA #9
syns$species_syn <- NA #10
syns$sciName_syn<- NA #11
syns$fam_syn <- NA#12 
syns$synType <- NA  #13

for(i in 1:nrow(syns)){
  query <- syns$TropID_Final[i]
  if(!is.na(query)){
    tropID <- tp_search(nameid =query, type="exact")
    tropname <- tropID[1,2]
    if(!tropname == syns$speciesFinal[1]){
      syns[i,9] <- paste(tropID[1,1])
      syns[i, 10] <- paste(tropID[1,2])
      syns[i, 11] <- paste(tropID[1,3])
      syns[i, 12] <- paste(tropID[1,4])
      syns[i, 13] <- "alternate spelling"
    }
  }
  print(paste("species",i,"done",sep=" "))
}
syns_save <- syns

altsyns <- syns %>% filter(!is.na(synType)) %>% select(V1, ID_syn, species_syn, sciName_syn, fam_syn, synType, TropID_Final)



#synsof NCBI (i.e. synsyns)
synsofNCBI <- NULL
for(i in 1:nrow(syns)){
  query <- syns$TropID_Final[i]
  if(!is.na(query)){
    syn <- synonyms(query, db="tropicos")
    synsdf <- syn[[1]]
    synsdf$ID <- query
    synsdf$V1 <- paste(syns$V1[i])
    synsofNCBI <- rbind(synsdf, synsofNCBI)
  }
  print(paste("values for species",i,"done",sep=" "))
}
synsofNCBI$synType <- "NCBI_syn"
write.csv(synsofNCBI, "synsofNCBI.csv", row.names=F)
synsyns <- data.frame(V1 = as.integer(synsofNCBI$V1), 
                      ID_syn = synsofNCBI$nameid, 
                      species_syn = synsofNCBI$scientificname, 
                      sciName_syn = synsofNCBI$scientificnamewithauthors, 
                      fam_syn = synsofNCBI$family,
                      synType= synsofNCBI$synType,
                      TropID_Final=synsofNCBI$ID)


accsyns <- species %>% filter(!tropAccName=="No accepted names were found")
accsyns <- data.frame(V1 = accsyns$V1, 
                      TropID_Final = accsyns$tropID_NCBI,
                      ID_syn = accsyns$tropAccID,
                      species_syn = accsyns$tropAccName,
                      sciName_syn = accsyns$sciAccName,
                      fam_syn = accsyns$tropAccFam)
accsyns$synType <- "acc"

syns_of_accs <-  NULL
for(i in 1:nrow(accsyns)){
  query <- accsyns$ID_syn[i]
  if(!is.na(query)){
    syn <- synonyms(query, db="tropicos")
    synsdf <- syn[[1]]
    synsdf$V1 <- paste(accsyns$V1[i])
    syns_of_accs <- rbind(synsdf, syns_of_accs)
  }
  print(paste("values for species",i,"done",sep=" "))
}
syns_of_accs$synType <- "syn of acc"
synaccs <- data.frame(V1 = as.integer(syns_of_accs$V1), 
                      ID_syn = syns_of_accs$nameid, 
                      species_syn = syns_of_accs$scientificname, 
                      sciName_syn = syns_of_accs$scientificnamewithauthors, 
                      fam_syn = syns_of_accs$family,
                      synType= syns_of_accs$synType)
V1 <- species %>% mutate(TropID_Final = tropID_NCBI) %>% select(V1, TropID_Final)
synaccs2 <- left_join(synaccs, V1)

#synsyns, accsyns, synaccs2, and altsyns
head(accsyns)
head(synsyns)
x <-rbind(synsyns, accsyns)
y <- rbind(synaccs2, x)
z <- rbind(altsyns, y)
joinme <- syns %>% select(V1, TropID_Final, speciesFinal, sciNameFinal, tropFamFinal, ID_OG, species_OG, sciName_OG)
donesyns <- full_join(joinme, z) %>% filter(!is.na(synType)) %>% arrange(V1)
write.csv(donesyns, "synonyms_almost.csv", row.names=F)



###Syns CHECK
rm(list = ls())
synonyms <- read.csv("synonyms_almost.csv")
synonyms <- synonyms[!grepl("cf.", synonyms$speciesFinal),]
#FIX FOR VARIETIES AND SUBSP.
var <- synonyms[grepl(" var. ", synonyms$species_OG),]
subsp <- synonyms[grepl(" subsp. ", synonyms$species_OG),]
varsnsubs <- rbind(var, subsp) 
not <- anti_join(synonyms, varsnsubs)
check <- rbind(not, varsnsubs)
varsnsubs <- varsnsubs %>% select(!speciesFinal) %>% select(!TropID_Final) %>% select(!sciNameFinal)
varsnsubs$speciesFinal <- varsnsubs$species_OG
varsnsubs$TropID_Final <- varsnsubs$ID_OG
varsnsubs$sciNameFinal <- varsnsubs$sciName_OG
synonyms <- rbind(not, varsnsubs)



both2 <- data.frame(both =intersect(unique(syn$ID_syn), unique(syn$TropID_Final)))



syn2 <- syn %>% distinct(V1, TropID_Final, speciesFinal, sciNameFinal, tropFamFinal, ID_OG, species_OG, sciName_OG, ID_syn, species_syn, sciName_syn, fam_syn, check_sp, check_sci, .keep_all = T)
both2 <- data.frame(both =intersect(unique(syn2$ID_syn), unique(syn2$TropID_Final)))
bad2 <- data.frame(ID_syn = both2$both)
bad2 <- left_join(bad2, syn2) %>% select(ID_syn, species_syn, sciName_syn) %>% distinct()
write.csv(syn2, "fixmesyns2.csv", row.names=F)
write.csv(bad2, "bad2.csv", row.names=F)

bad <- data.frame(ID_syn = both2$both)
bad1 <- data.frame(TropID_Final = both2$both)
bad2 <- left_join(bad, syn2) %>% distinct()
bad3 <- left_join(bad1, syn2) %>% distinct()
baddest <- rbind(bad2, bad3)  %>% distinct() %>% arrange(speciesFinal)
write.csv(baddest, "baddest_old.csv", row.names=F)

notbad <- anti_join(syn2, baddest)
check <- rbind(notbad, baddest)

#CHECK BADDEST
baddest_done <- read.csv("baddest.csv")
#write.csv(baddest_done, "baddest_justincase.csv", row.names=F)
#write.csv(notbad, "notbad.csv", row.names = F)
notbad <- read.csv("notbad.csv")
#check <- data.frame(both =intersect(unique(baddest_done$ID_syn), unique(baddest_done$TropID_Final)))
#check1 <- data.frame(ID_syn = check$both)
#check2 <- data.frame(TropID_Final = check$both)
#check3 <- left_join(check1, baddest_done) %>% distinct()
#check4 <- left_join(check2, baddest_done) %>% distinct()
#check5 <- rbind(check3, check4)

#badcheck <- check5 %>% filter(!species_syn == speciesFinal)
ref <- baddest_done %>% select(speciesFinal) %>% distinct()
badtropID <- baddest_done %>% select(TropID_Final, speciesFinal) %>% distinct()
badtropID$dup <- duplicated(badtropID$speciesFinal)
badtropID <- badtropID %>% filter(dup == TRUE) %>% select(speciesFinal) %>% distinct()

#now the same thing for the OG names
ref <- baddest_done %>% select(species_OG) %>% distinct()
badtropID <- baddest_done %>% select(ID_OG, species_OG) %>% distinct()
badtropID$dup <- duplicated(badtropID$species_OG)
badtropID <- badtropID %>% filter(dup == TRUE) %>% select(species_OG) %>% distinct()

#combine them lolz
#syns <- rbind(notbad, baddest_done) %>% arrange(speciesFinal)
write.csv(syns, "almostdonesyns.csv", row.names=F)
syns <- read.csv("almostdonesyns.csv")




#nice


#alrighty....check dups
dup_syn <- syns %>% select(speciesFinal, species_syn, ID_syn) %>% distinct()
dup_syn$dup <- duplicated(dup_syn$ID_syn)
bad_syn <- dup_syn %>% filter(dup == TRUE) %>% select(ID_syn) %>% distinct()
fixme <- left_join(bad_syn, syns) %>% arrange(speciesFinal)
fine <- anti_join(syns, fixme)
check <- rbind(fixme, fine)


write.csv(fixme, "morebad.csv", row.names = F)
write.csv(fine, "fine.csv", row.names = F)


morebad <- read.csv("morebad.csv")
done <- read.csv("fine.csv")

synz <- rbind(morebad, done) %>% select(!synType)  %>% select(!check_sci) %>% select(!check_sp) %>% distinct()
#double check
dup <- synz %>% select(TropID_Final, ID_syn) %>% distinct()
dup$dup <- duplicated(dup$ID_syn)
bad_syn <- dup %>% filter(dup == TRUE) %>% select(ID_syn) %>% distinct()

synz$species_OG_hmm <- synz$species_OG
synz <- synz %>% select(!species_OG)
sp <- read.csv("species_ids_cleaned.csv")
spp <- data.frame(species_OG = sp$species)
spp$V1 = sp$V1
amidoneyet <- left_join(synz, spp)
write.csv(amidoneyet, "synonyms_done.csv")

#alrighty. NOW.
rename_seqs <- amidoneyet %>% select(V1, TropID_Final, speciesFinal, ID_OG, species_OG,  species_OG_hmm) %>% distinct()
#okay wtf happened

okay <- data.frame(V1 = sp$V1, 
                   ID_OG = sp$tropID_NCBI, 
                   species_OG_hmm = sp$NCBIspecies, 
                   sciName_OG = sp$sciName_NCBI, 
                   species_OG = sp$species,
                   speciesFinal = sp$speciesFinal,
                   sciNameFinal = sp$sciNameFinal,
                   TropID_Final = sp$TropID_Final,
                   tropFamFinal = sp$tropFamFinal)


hmm <- anti_join(spp, rename_seqs)
hmm <- hmm[!grepl("cf.", hmm$species_OG),]
hmm <- hmm[!grepl("aff.", hmm$species_OG),]
hmm <- hmm[!grepl("sp.", hmm$species_OG),]

damn <- left_join(hmm, okay)
damn[2,3] <- 100477246
damn[12,3] <- 35190163


syns_new <- NULL
for(i in 1:nrow(damn)){
  query <- damn[i,3]
  V1 <- damn[i, 2]
  if(!is.na(query)){
    syns <- synonyms(query, db="tropicos")
    synsdf <- syns[[1]]
    synsdf$ID <- query
    synsdf$V1 <- V1
    syns_new <- rbind(synsdf, syns_new)
    print(paste("values for species",i,"done",sep=" ")) 
  }
}

syns_new_save <- syns_new
syns_new$syntype <- paste("syn")



acc_new <- NULL
for(i in 1:nrow(damn)){
  query <- damn[i,3]
  V1 <- damn[i, 2]
  if(!is.na(query)){
    acc <- tp_accnames(query)
    synsofacc <- acc[[1]]
    if(!synsofacc == "No accepted names were found"){
      synsofacc$ID <- query
      synsofacc$syntype <- paste("synofacc")
      synsofacc$V1 <- V1
      accsdf <- acc[[2]]
      accsdf$ID <- query
      accsdf$syntype <- paste("acc")
      accsdf$V1 <- V1
      accs <- bind_rows(accsdf, synsofacc)
      acc_new <- bind_rows(accs, acc_new)
    }
    print(paste("values for species",i,"done",sep=" ")) 
  }
}

acc_new_save <- acc_new


both <- bind_rows(syns_new, acc_new)
both2 <- data.frame(V1 = both$V1, species_syn = both$scientificname, sciName_syn=both$scientificnamewithauthors, fam_syn = both$family, ID_syn = both$nameid, syntype = both$syntype)
wowsoclose <- full_join(damn, both2) %>% arrange(species_syn)
wowsoclose_save <- wowsoclose
write.csv(wowsoclose, "SODAMNCLOSE.csv", row.names=F)
#check <- data.frame(both =intersect(unique(wowsoclose$ID_syn), unique(wowsoclose$TropID_Final)))

close <- read.csv("SODAMNCLOSE.csv") %>% select(!syntype) %>% distinct()

#didn't catch all syns of accs
moresyns <- NULL
for(i in 1:nrow(close)){
  query <- close[i,2]
  V1 <- close[i, 1]
  if(!is.na(query)){
    syns <- synonyms(query, db="tropicos")
    synsdf <- syns[[1]]
    synsdf$ID <- query
    synsdf$V1 <- V1
    moresyns <- rbind(synsdf, moresyns)
    print(paste("values for species",i,"done",sep=" ")) 
  }
}

more <- data.frame(V1 = moresyns$V1, species_syn = moresyns$scientificname, sciName_syn=moresyns$scientificnamewithauthors, fam_syn = moresyns$family, ID_syn = moresyns$nameid)
ugh <- left_join(more, damn)

again <- full_join(ugh, close) %>% distinct()
#intersect <- data.frame(both =intersect(unique(again$ID_syn), unique(again$TropID_Final)))




#combine and check it:
godwhy <- read.csv("synonyms_done.csv")
wow <- bind_rows(godwhy, again) %>% select(!X)
rename_seqs <- wow %>% select(V1, TropID_Final, speciesFinal, ID_OG, species_OG,  species_OG_hmm) %>% distinct()

hmm <- anti_join(spp, rename_seqs)

###NOOOO
hmm <- hmm[grepl("sp.", hmm$species_OG),]
hmm <- hmm[!grepl("cf.", hmm$species_OG),]

damn2 <- left_join(hmm, okay)
damn2 <- damn2[-c(2,3,4,5),]


syns_new2 <- NULL
for(i in 1:nrow(damn2)){
  query <- damn2[i,3]
  V1 <- damn2[i, 2]
  if(!is.na(query)){
    syns <- synonyms(query, db="tropicos")
    synsdf <- syns[[1]]
    synsdf$ID <- query
    synsdf$V1 <- V1
    syns_new2 <- rbind(synsdf, syns_new2)
    print(paste("values for species",i,"done",sep=" ")) 
  }
}
syns_new$syntype <- paste("syn")
acc_new2 <- NULL
for(i in 1:nrow(damn2)){
  query <- damn2[i,3]
  V1 <- damn2[i, 2]
  if(!is.na(query)){
    acc <- tp_accnames(query)
    synsofacc <- acc[[1]]
    if(!synsofacc == "No accepted names were found"){
      synsofacc$ID <- query
      synsofacc$syntype <- paste("synofacc")
      synsofacc$V1 <- V1
      accsdf <- acc[[2]]
      accsdf$ID <- query
      accsdf$syntype <- paste("acc")
      accsdf$V1 <- V1
      accs <- bind_rows(accsdf, synsofacc)
      acc_new2 <- bind_rows(accs, acc_new2)
    }
    print(paste("values for species",i,"done",sep=" ")) 
  }
}
moresyns2 <- NULL
for(i in 1:nrow(damn2)){
  query <- damn2[i,8]
  V1 <- damn2[i, 2]
  if(!is.na(query)){
    syns <- synonyms(query, db="tropicos")
    synsdf <- syns[[1]]
    synsdf$ID <- query
    synsdf$V1 <- V1
    moresyns2 <- rbind(synsdf, moresyns2)
    print(paste("values for species",i,"done",sep=" ")) 
  }
}
both_ <- bind_rows(syns_new2, moresyns2)
both2_ <- data.frame(V1 = both_$V1, species_syn = both_$scientificname, sciName_syn=both_$scientificnamewithauthors, fam_syn = both_$family, ID_syn = both_$nameid)
no <- full_join(damn2, both2_) %>% arrange(species_syn)
#combine and check it... again.
done <- bind_rows(wow, no) %>% arrange(V1)
write.csv(done, "synonyms_done.csv", row.names = F)
done <- read.csv("synonyms_done.csv") %>% arrange(speciesFinal)
write.csv(done, "synonyms_done.csv", row.names = F)

#check:
dup <- done %>% select(TropID_Final, ID_syn) %>% distinct()
dup$dup <- duplicated(dup$ID_syn)
bad_syn <- dup %>% filter(dup == TRUE) %>% select(ID_syn) %>% distinct()







done <- read.csv("synonyms_done.csv") %>% select(!x)
check <- done %>% select(V1, species_OG) %>% distinct()
check$dup <- duplicated(check$species_OG)
#syns_only <- done %>% select(V1, TropID_Final, speciesFinal, species_syn, ID_syn) %>% distinct()
#rename_seqs <- done %>% select(V1, TropID_Final, speciesFinal, ID_OG, species_OG, species_OG_hmm) %>% distinct()
df <- data.frame(V1 = done$V1, TropID_Final = done$TropID_Final, speciesFinal = done$speciesFinal)


#FINAL CHECKz...
V1s <- doneFin %>% select(V1, speciesFinal) %>% distinct()
V1s$dup <- duplicated(V1s$V1)


finalz <- done %>% select(V1, TropID_Final, speciesFinal, sciNameFinal, species_OG, species_OG_hmm, ID_OG, sciName_OG, tropFamFinal) %>% distinct()
fin <- data.frame(species_OG = finalz$species_OG,
                  V1=finalz$V1, 
                  TropID_Final= finalz$TropID_Final, 
                  speciesFinal = finalz$speciesFinal, 
                  sciNameFinal= finalz$sciNameFinal, 
                  tropFamFinal= finalz$tropFamFinal, 
                  ID_OG = finalz$ID_OG, 
                  species_OG_hmm =finalz$species_OG_hmm, 
                  sciName_OG= finalz$sciName_OG,
                  ID_syn = finalz$TropID_Final, 
                  species_syn = finalz$speciesFinal, 
                  sciName_syn= finalz$sciNameFinal, 
                  fam_syn= finalz$tropFamFinal)
doneFin <- bind_rows(fin, done) %>% distinct()
write.csv(doneFin, "synonyms_done.csv", row.names = F) 
doneFin <- read.csv("synonyms_done.csv") #%>% arrange(speciesFinal)

#IGNORED PROBLEM CHILDREN e.g. Phaeoceros microsporus
doneFinFIN <- doneFin[-c(568, 1731, 2921, 3156, 5572, 7918,4727, 4718, 8770, 8824, 8825, 8826, 9538, 9840, 11155, 11164, 11172, 12495, 13595, 14673, 14705, 18900, 19127, 19167, 19138, 19303, 19304, 19830, 19953), ]
write.csv(doneFinFIN, "synonyms_done.csv", row.names = F) 

library(openxlsx)
library(dplyr)
setwd("/Users/ranunculus/Documents/bigBryogeography/")
doneFin <- read.xlsx("synonyms_done.xlsx")
syns_only <- doneFin %>% select(V1, TropID_Final, speciesFinal,sciNameFinal, species_syn, sciName_syn, ID_syn) %>% distinct()
write.xlsx(syns_only, "syns_done.xlsx", row.names=F)
rename_seqs <- doneFin %>% select(V1, TropID_Final, speciesFinal, ID_OG, species_OG, species_OG_hmm) %>% distinct()
write.xlsx(rename_seqs, "rename_seqs.xlsx", row.names=F)

