acc <- read.csv("accessions.csv") #see bryogegraphy.sh
acc$species <- NA
acc$product <- NA
acc$gene <- NA
acc$note <- NA


for (i in 1:nrow(acc)){
  print(i)
  query<-paste(acc$accession[i])
  see <- entrez_fetch(db="nuccore", id=query, rettype="gb")
  species="/organism=(.*?)\n"
  resultSp <- regmatches(see,regexec(species,see))
  dfsp <- data.frame(authors=resultSp)
  acc[i,2]<-paste(dfsp[2,1])
  product="/product=(.*?)\n"
  resultproduct <- regmatches(see,regexec(product,see))
  dfproduct <- data.frame(authors=resultproduct)
  acc[i,3]<-paste(dfproduct[2,1])
  gene="/gene=(.*?)\n"
  resultgene <- regmatches(see,regexec(gene,see))
  dfgene <- data.frame(authors=resultgene)
  acc[i,4]<-paste(dfgene[2,1])
  note="/note=(.*?)\n"
  resultnote <- regmatches(see,regexec(note,see))
  dfnote <- data.frame(authors=resultnote)
  acc[i,5]<-paste(dfnote[2,1])
  print(paste("values for accession",i,"calculated",sep=" "))
}


acc2 <- as.data.frame(sapply(acc, function(x) gsub("\"", "", x)))
acc2[] <- lapply(acc2, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})
acc2$final <- NA
for (i in 1:nrow(acc2)){
  query<-paste(acc2$gene[i])
  query2<-paste(acc2$product[i])
  query3<-paste(acc2$note[i])
  if(!query=="NA"){
    acc2[i, 6]<-query
  }
  if(query=="NA"){
    acc2[i,6]<-query2
  }
  if(query2=="NA" & query=="NA"){
    acc2[i, 6]<-query3
  }
}

seq <- read.csv("seq.csv", header=F)
fix <- as.data.frame(sapply(seq, function(x) gsub(" ", "", x)))
seqs <- cbind(acc, fix) %>% select(accession, V1)
chars <- as.vector(seqs$V1)
count <- data.frame(chars = nchar(chars, type = "chars"))
fin <- cbind(seqs, count)   #fin <- cbind(seqs, count)
head(fin)
write.csv(fin, "seqs.csv", row.names=F)

seqs <- read.csv("seqs.csv") %>% select(!V1)
accs <- full_join(acc2, seqs)
accs[] <- lapply(accs, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})
accs2 <- accs %>% filter(!final=="NA")
write.csv(accs, file="accessions_genes.csv", row.names=F)

#################################################################
library(rentrez)
library(dplyr)
library(reshape2)
#####FIX ROSE: #####
rose <- read.csv("RoseS1.csv") # From Rose et al. (2016), table S1
r1 <- melt(rose, id.vars = "Organism", variable.name = "gene",value.name = "accession") %>% filter(!is.na(accession)) %>% select(accession)
r1 <- r1 %>% mutate(accession = as.character(accession))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "APNCPRBCLB", "D43696"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "APU87063", "U87063"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "BTHCPRBCL", "L13475"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "MGACPRBCL", "L13481"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "PCU87087", "U87087"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "TPU87091", "U87091"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "AY608581", "AH013789"))
r1 <- r1 %>% mutate(accession = replace(accession, accession ==  "AY608579", "AH013788"))
r <- r %>% mutate(accession = replace(accession, accession ==  "NC_013765 2", "NC_013765"))

##########get info ##################
accession <- read.csv("accessions_genes.csv") %>% select(accession)
acc <- rbind(accession, r1) %>% distinct()
acc$species <- NA
acc$product <- NA
acc$gene <- NA
acc$note <- NA


for (i in 1:nrow(acc)){
  print(i)
  query<-paste(acc$accession[i])
  see <- entrez_fetch(db="nuccore", id=query, rettype="gb")
  species="/organism=(.*?)\n"
  resultSp <- regmatches(see,regexec(species,see))
  dfsp <- data.frame(authors=resultSp)
  acc[i,2]<-paste(dfsp[2,1])
  product="/product=(.*?)\n"
  resultproduct <- regmatches(see,regexec(product,see))
  dfproduct <- data.frame(authors=resultproduct)
  acc[i,3]<-paste(dfproduct[2,1])
  gene="/gene=(.*?)\n"
  resultgene <- regmatches(see,regexec(gene,see))
  dfgene <- data.frame(authors=resultgene)
  acc[i,4]<-paste(dfgene[2,1])
  note="/note=(.*?)\n"
  resultnote <- regmatches(see,regexec(note,see))
  dfnote <- data.frame(authors=resultnote)
  acc[i,5]<-paste(dfnote[2,1])
  print(paste("values for accession",i,"calculated",sep=" "))
}

acc2 <- as.data.frame(sapply(acc, function(x) gsub("\"", "", x)))
acc2[] <- lapply(acc2, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})

acc2$final <- NA
for (i in 1:nrow(acc2)){
  query<-paste(acc2$gene[i])
  query2<-paste(acc2$product[i])
  query3<-paste(acc2$note[i])
  if(!query=="NA"){
    acc2[i, 6]<-query
  }
  if(query=="NA"){
    acc2[i,6]<-query2
  }
  if(query2=="NA" & query=="NA"){
    acc2[i, 6]<-query3
  }
}

#get rose seqs based on accessions: #####
#rose <- read.csv("RoseS1.csv")
#r1 <- melt(rose, id.vars = "Organism", variable.name = "gene",value.name = "accession") %>% filter(!is.na(accession)) %>% select(accession)
recs <- NULL
for(i in 1:nrow(r1)){
  print(i)
  query<-paste(r1$accession[i])
  b <- entrez_fetch(db="nuccore", id=query, rettype = "fasta")
  #b <- gsub(pattern = ">", replacement = paste(">",query, sep = ""), b)
  recs <- c(recs,b)
}
as.character(recs) -> char_recs
write(char_recs, "rose.fasta")

################get counts from seq#############
seq <- read.csv("seq.csv", header=F)
fix <- as.data.frame(sapply(seq, function(x) gsub(" ", "", x)))
chars <- as.vector(fix$V1)
count <- data.frame(chars = nchar(chars, type = "chars"))
counts <- cbind(fix, count)
nam <- read.csv("seqNam.csv", header=F)
colnames(nam)[colnames(nam)=="V1"] <- "accession"
nams <- as.data.frame(sapply(nam, function(x) gsub(">", "", x)))
done <- cbind(nams, counts)


#now rose:
Roseq <- read.csv("Roseq.csv", header=F)
Rofix <- as.data.frame(sapply(Roseq, function(x) gsub(" ", "", x)))
Rochars <- as.vector(Rofix$V1)
Rocount <- data.frame(Rochars = nchar(Rochars, type = "chars")) %>% rename(chars = Rochars)
Rocounts <- cbind(Rofix, Rocount)
Ronam <- read.csv("RoseqNam.csv", header=F)
colnames(Ronam)[colnames(Ronam)=="V1"] <- "accession"
Ronams <- as.data.frame(sapply(Ronam, function(x) gsub(">", "", x)))
Ronams <- Ronams %>% mutate(accession = as.character(accession))
Ronams <- as.data.frame(substr(Ronams$accession,1,nchar(Ronams$accession)-2))
colnames(Ronams)[colnames(Ronams)=="substr(Ronams$accession, 1, nchar(Ronams$accession) - 2)"] <- "accession"
Rodone <- cbind(Rocounts, Ronams)

#combine them:
bound <- rbind(Rodone, done) %>% select(accession, chars)
#continue:

write.csv(bound, "seqs.csv", row.names=F)

####################add counts to accs######################
bound <- read.csv("seqs.csv")
seqs <- left_join(acc2, bound)
accs <- seqs %>% distinct()
accs[] <- lapply(accs, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})
#check and fix!
accs2 <- accs %>% filter(!final=="NA")
write.csv(accs, file="accessions_genes.csv", row.names=F)


#####clean, RENAME##############
rm(list = ls())
library(stringr)
accs <- read.csv("accessions_genes.csv")
noP <-  accs %>% filter(!species == "Physcomitrella patens")
#accs <- accs %>% distinct()
checkGenes <- count(noP, final)
rose <- read.csv("roseAcc.csv")
roseGenes <- rose %>% select(gene) %>% distinct()
sapply(accs, class)
accs <- accs %>% mutate(final = as.character(final))
accs <- accs %>% mutate(final = replace(final, final ==  "[BankIt_ITS_wizard]; [rRNAITS_notfound];", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "26S ribosomal RNA", "X26S"))
accs <- accs %>% mutate(final = replace(final, final ==  "26S rRNA", "X26S"))
accs <- accs %>% mutate(final = replace(final, final ==  "26S large subunit ribosomal RNA", "X26S"))
accs <- accs %>% mutate(final = replace(final, final ==  "5.8S ribosomal RNA", "5.8s"))
accs <- accs %>% mutate(final = replace(final, final ==  "5.8S rRNA", "5.8s"))
accs <- accs %>% mutate(final = replace(final, final ==  "5.8S rRNA", "5.8s"))
accs <- accs %>% mutate(final = replace(final, final ==  "60S ribosomal protein L23a-like protein", "60S ribosomal protein L23A-like protein"))
accs <- accs %>% mutate(final = replace(final, final ==  "Aldh10A1", "ALDH10A1"))
accs <- accs %>% mutate(final = replace(final, final ==  "Aldh11A1", "ALDH11A1"))
accs <- accs %>% mutate(final = replace(final, final ==  "Aldh2B2", "ALDH2B2"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB-rbcL intergenic spacer region", "atpB.rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpI-atpH intergenic spacer region", "atpB.rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpI-atpH intergenic spacer, partial sequence", "atpI-atpH intergenic spacer"))
accs <- accs %>% mutate(final = replace(final, final ==  "authority: Calymperes bartramii Reese", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "authority: Calymperes rubiginosum (Mitt.) Reese", "trnL.trnF"))
accs$final[accs$accession == "KJ488179.1"] <- "ITS"
accs$final[accs$accession == "KJ488295.1"] <- "atpB"
accs <- accs %>% mutate(final = replace(final, final ==  "authority: Octoblepharum pulvinatum (Dozy & Molk.)", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "authority: Octoblepharum stramineum Mitt.", "trnL.trnF"))
accs$final[accs$accession == "MN794384.1"] <- "ITS"
accs$final[accs$accession == "MN794385.1"] <- "ITS"
accs$final[accs$accession == "MN794386.1"] <- "ITS"
accs$final[accs$accession == "MK371359.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371360.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371361.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371362.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371363.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371364.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371365.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371366.1"] <- "trnL.trnF"
accs$final[accs$accession == "MK371367.1"] <- "trnL.trnF"
accs$final[accs$accession == "KJ488180.1"] <- "ITS"
accs$final[accs$accession == "KJ488181.1"] <- "ITS"
accs$final[accs$accession == "KJ488182.1"] <- "ITS"
accs$final[accs$accession == "KJ488183.1"] <- "ITS"
accs$final[accs$accession == "KJ488184.1"] <- "ITS"
accs$final[accs$accession == "KJ488187.1"] <- "ITS"
accs$final[accs$accession == "KJ488188.1"] <- "ITS"
accs$final[accs$accession == "KJ488189.1"] <- "ITS"
accs$final[accs$accession == "KJ488190.1"] <- "ITS"
accs$final[accs$accession == "KJ488191.1"] <- "ITS"
accs$final[accs$accession == "KJ488192.1"] <- "ITS"
accs$final[accs$accession == "KJ488193.1"] <- "ITS"
accs$final[accs$accession == "KJ488195.1"] <- "ITS"
accs$final[accs$accession == "KJ488196.1"] <- "ITS"
accs$final[accs$accession == "KJ488197.1"] <- "ITS"
accs$final[accs$accession == "KJ488198.1"] <- "ITS"
accs$final[accs$accession == "KJ488199.1"] <- "ITS"
accs$final[accs$accession == "KJ488200.1"] <- "ITS"
accs$final[accs$accession == "KJ488201.1"] <- "ITS"
accs$final[accs$accession == "KJ488202.1"] <- "ITS"
accs$final[accs$accession == "KJ488203.1"] <- "ITS"
accs$final[accs$accession == "KJ488204.1"] <- "ITS"
accs$final[accs$accession == "KJ488205.1"] <- "ITS"
accs$final[accs$accession == "KJ488206.1"] <- "ITS"
accs$final[accs$accession == "KJ488207.1"] <- "ITS"
accs$final[accs$accession == "KJ488208.1"] <- "ITS"
accs$final[accs$accession == "KJ488209.1"] <- "ITS"
accs$final[accs$accession == "KJ488210.1"] <- "ITS"
accs$final[accs$accession == "KJ488211.1"] <- "ITS"
accs$final[accs$accession == "KJ488212.1"] <- "ITS"
accs$final[accs$accession == "KJ488213.1"] <- "ITS"
accs$final[accs$accession == "KJ488214.1"] <- "ITS"
accs$final[accs$accession == "KJ488185.1"] <- "ITS"
accs$final[accs$accession == "KJ488186.1"] <- "ITS"
accs$final[accs$accession == "KJ488215.1"] <- "ITS"
accs$final[accs$accession == "KJ488216.1"] <- "ITS"
accs$final[accs$accession == "KJ488217.1"] <- "ITS"
accs$final[accs$accession == "KJ488218.1"] <- "ITS"
accs$final[accs$accession == "KJ488219.1"] <- "ITS"
accs$final[accs$accession == "KJ488220.1"] <- "ITS"
accs$final[accs$accession == "KJ488221.1"] <- "ITS"
accs$final[accs$accession == "KJ488222.1"] <- "ITS"
accs$final[accs$accession == "KJ488223.1"] <- "ITS"
accs$final[accs$accession == "KJ488224.1"] <- "ITS"
accs$final[accs$accession == "KJ488225.1"] <- "ITS"
accs$final[accs$accession == "KJ488226.1"] <- "ITS"
accs$final[accs$accession == "KJ488227.1"] <- "ITS"
accs$final[accs$accession == "KJ488228.1"] <- "ITS"
accs$final[accs$accession == "KJ488229.1"] <- "ITS"
accs$final[accs$accession == "KJ488230.1"] <- "ITS"
accs$final[accs$accession == "KJ488231.1"] <- "ITS"
accs$final[accs$accession == "KJ488232.1"] <- "ITS"
accs$final[accs$accession == "KJ488233.1"] <- "ITS"
accs$final[accs$accession == "KJ488234.1"] <- "ITS"
accs$final[accs$accession == "KJ488235.1"] <- "ITS"
accs$final[accs$accession == "KJ488236.1"] <- "ITS"
accs$final[accs$accession == "KJ488237.1"] <- "ITS"
accs$final[accs$accession == "KJ488238.1"] <- "ITS"
accs$final[accs$accession == "KJ488239.1"] <- "ITS"
accs$final[accs$accession == "KJ488240.1"] <- "ITS"
accs$final[accs$accession == "KJ488241.1"] <- "ITS"
accs$final[accs$accession == "KJ488242.1"] <- "ITS"
accs$final[accs$accession == "KJ488243.1"] <- "ITS"
accs$final[accs$accession == "KJ488244.1"] <- "ITS"
accs$final[accs$accession == "KJ488245.1"] <- "ITS"
accs$final[accs$accession == "KJ488246.1"] <- "ITS"
accs$final[accs$accession == "KJ488247.1"] <- "ITS"
accs$final[accs$accession == "KJ488248.1"] <- "ITS"
accs$final[accs$accession == "KJ488249.1"] <- "ITS"
accs$final[accs$accession == "KJ488250.1"] <- "ITS"
accs$final[accs$accession == "KJ488251.1"] <- "ITS"
accs$final[accs$accession == "KJ488252.1"] <- "ITS"
accs$final[accs$accession == "KJ488253.1"] <- "ITS"
accs$final[accs$accession == "KJ488254.1"] <- "ITS"
accs$final[accs$accession == "KJ488255.1"] <- "ITS"
accs$final[accs$accession == "KJ488256.1"] <- "ITS"
accs$final[accs$accession == "KJ488257.1"] <- "ITS"
accs$final[accs$accession == "KJ488258.1"] <- "ITS"
accs$final[accs$accession == "KJ488259.1"] <- "ITS"
accs$final[accs$accession == "KJ488260.1"] <- "ITS"
accs$final[accs$accession == "KJ488261.1"] <- "ITS"
accs$final[accs$accession == "KJ488262.1"] <- "ITS"
accs$final[accs$accession == "KJ488263.1"] <- "ITS"
accs$final[accs$accession == "KJ488264.1"] <- "ITS"
accs$final[accs$accession == "KJ488265.1"] <- "ITS"
accs$final[accs$accession == "KJ488266.1"] <- "ITS"
accs$final[accs$accession == "KJ488267.1"] <- "ITS"
accs$final[accs$accession == "KJ488268.1"] <- "ITS"
accs$final[accs$accession == "KJ488269.1"] <- "ITS"
accs$final[accs$accession == "KJ488270.1"] <- "ITS"
accs$final[accs$accession == "KJ488271.1"] <- "ITS"
accs$final[accs$accession == "KJ488272.1"] <- "ITS"
accs$final[accs$accession == "KJ488273.1"] <- "ITS"
accs$final[accs$accession == "KJ488274.1"] <- "ITS"
accs$final[accs$accession == "KJ488275.1"] <- "ITS"
accs$final[accs$accession == "KJ488276.1"] <- "ITS"
accs$final[accs$accession == "KJ488277.1"] <- "ITS"
accs$final[accs$accession == "KJ488278.1"] <- "ITS"
accs$final[accs$accession == "KJ488279.1"] <- "ITS"
accs$final[accs$accession == "KJ488280.1"] <- "ITS"
accs$final[accs$accession == "KJ488281.1"] <- "ITS"
accs$final[accs$accession == "KJ488282.1"] <- "ITS"
accs$final[accs$accession == "KJ488283.1"] <- "ITS"
accs$final[accs$accession == "KJ488284.1"] <- "ITS"
accs$final[accs$accession == "KJ488285.1"] <- "ITS"
accs$final[accs$accession == "KJ488286.1"] <- "ITS"
accs$final[accs$accession == "KJ488287.1"] <- "ITS"
accs$final[accs$accession == "KJ488288.1"] <- "ITS"
accs$final[accs$accession == "KJ488289.1"] <- "ITS"
accs$final[accs$accession == "KJ488290.1"] <- "ITS"
accs$final[accs$accession == "KJ488291.1"] <- "ITS"
accs$final[accs$accession == "KJ488292.1"] <- "ITS"
accs$final[accs$accession == "KJ488293.1"] <- "ITS"
accs <- accs %>% mutate(final = replace(final, final ==  "conatins tRNA-Leu (trnL), trnL-trnF intergenic", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contain trnL, trnL-F IGS and trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains 18S ribosomal RNA, internal transcribed", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains 5.8S ribosomal RNA, and internal", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains 5.8S ribosomal RNA, internal transcribed", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains ATP synthase CF0 subunit III (atpH),", "atpH"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains ATP synthase CF0 subunit IV (atpI),", "atpI"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains atpB gene, atpB-rbcL intergenic spacer,", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains internal transcribed spacer 1 and 5.8S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains internal transcribed spacer 1, 5.8S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains psbD-trnT intergenic spacer and tRNA-Thr", "trnT.psbD"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains rps4 and rps4-trnS intergenic spacer", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains rps4 gene, rps4-trnT intergenic spacer,", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains small subunit ribosomal RNA and internal", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains small subunit ribosomal RNA, internal", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu (trnL) and trnL-trnF intergenic", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu (trnL) gene, trnL-trnF intergenic", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu (trnL), trnL-trnF intergenic", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu, trnL-trnF intergenic spacer, and", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnL, trnL-trnF intergenic spacer and", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnL,trnL-trnF intergenic spacer and trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "cotains trnL gene, trnL-trnF intergenic spacer and", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing NADH dehydrogenase subunit", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing partial nad5 gene for NADH", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing photosystem II subunit PsbT", "psbT"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing ribulose-1,5-bisphosphate", "rbcl"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing RNA polymerase C (rpoC1)", "rpoC1"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing rps4 gene + trnS-rps4 spacer", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing rps4 region;", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing rps4-trnL-F region;", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing rps4-trnT region;", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing rps4-trnT-trnL-trnF region;", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing tRNA-Gly (trnG) gene,", "trnG"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing trnT-trnL region;", "trnT.trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing trnT-trnL-trnF region;", "trnT"))
accs <- accs %>% mutate(final = replace(final, final ==  "genotype: ITS1&2", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "genotype: trnL", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "internal transcribed spacer 1", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "internal transcribed spacer 2", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "NAD5", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "Nad5", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "name Nad5", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "phosphoadenosine phosphosulfate reductase", "cysH"))
accs <- accs %>% mutate(final = replace(final, final ==  "phosphoadenosine phosphosulphate reductase", "cysH"))
accs <- accs %>% mutate(final = replace(final, final ==  "phosphoadenosine phosphosulfate reductase", "cysH"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbA-trnH intergenic spacer region; may also", "psbA.trnH"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbD-trnT intergenic spacer", "trnT.psbD"))
accs <- accs %>% mutate(final = replace(final, final ==  "rbcl", "rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "rbcL-atpB intergenic spacer", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "rbcL-atpB intergenic spacer", "rbcL.atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribulose-1,5-bisphosphate carboxylase/oxygenase", "rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "rpl16 intron", "rpl16"))
accs <- accs %>% mutate(final = replace(final, final ==  "sequence contains 18S rRNA gene, ITS1, 5.8S rRNA", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "sequence contains ITS1, 5.8S rRNA gene, ITS2 and", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "sequence contains ITS1, 5.8S rRNA gene, ITS2, 28S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "small ribosomal protein 4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "small subunit ribosomal RNA", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "small subunit ribosomal RNA", "ITS"))
accs$final[accs$accession == "MG584705.1"] <- "trnL.trnF"
accs <- accs %>% mutate(final = replace(final, final ==  "tRNA-His", "trnH"))
accs <- accs %>% mutate(final = replace(final, final ==  "tRNA-Leu", "trnL"))
accs$final[accs$accession == "KX396249.1"] <- "trnL.trnF"
accs <- accs %>% mutate(final = replace(final, final ==  "trnF-GAA", "trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnH-psbA intergenic spacer", "trnH.psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnI (CAU)", "trnI"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnI(cau)", "trnI"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnI-CAU", "trnI"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnI(CAU)", "trnI"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnK-UUU", "trnK"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL-trnF intergenic spacer", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL-trnF intergenic spacer region", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL-UAA", "trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnM (CAT)", "trnM"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnM-CAT", "trnM"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnM-CAU", "trnM"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnS-GGA", "trnS"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnT-psbD intergenic spacer", "trnT.psbD"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnT-trnE intergenic spacer", "trnE"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnV-UAC", "trnV"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpI-atpH intergenic spacer", "atpH.atpI"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpH-atpI intergenic spacer", "atpH.atpI"))
accs <- accs %>% mutate(final = replace(final, final ==  "18S rRNA", "18S"))
accs <- accs %>% mutate(final = replace(final, final ==  "rps7-atp6 intergenic spacer region", "rps7.atp6"))
accs <- accs %>% mutate(final = replace(final, final ==  "18S ribosomal RNA", "18S"))
accs <- accs %>% mutate(final = replace(final, final ==  "23S large subunit ribosomal RNA", "23S"))
accs <- accs %>% mutate(final = replace(final, final ==  "23S ribosomal RNA", "23S"))
accs$final[accs$accession == "AY607999"] <- "psbT"
accs$final[accs$accession == "AY608120"] <- "trnL.trnF"
accs$final[accs$accession == "JN089198"] <- "atpB.rbcL"
accs$final[accs$accession == "AY857611"] <- "ITS"
accs$final[accs$accession == "AY857613"] <- "ITS"
accs$final[accs$accession == "JN089244"] <- "rps7.atp6"
accs$final[accs$accession == "AY950399"] <- "trnL.trnF"
accs$final[accs$accession == "AF403632"] <- "ITS"
accs$final[accs$accession == "AY608002"] <- "psbT"
accs$final[accs$accession == "AY950401"] <- "trnL.trnF"
accs$final[accs$accession == "AY950402"] <- "trnL.trnF"
accs$final[accs$accession == "GU552296"] <- "trnL.trnF"
accs$final[accs$accession == "GU552284"] <- "trnL.trnF"
accs$final[accs$accession == "GU552295"] <- "trnL.trnF"
accs$final[accs$accession == "GU563727"] <- "trnL.trnF"
accs$final[accs$accession == "GU552282"] <- "trnL.trnF"
accs$final[accs$accession == "GU552279"] <- "trnL.trnF"
accs$final[accs$accession == "AY857596"] <- "ITS"
accs$final[accs$accession == "AY857553"] <- "trnL.trnF"
accs$final[accs$accession == "AY857597"] <- "ITS"
accs$final[accs$accession == "AY857554"] <- "trnL.trnF"
accs$final[accs$accession == "AY857598"] <- "ITS"
accs$final[accs$accession == "AY857555"] <- "trnL.trnF"
accs$final[accs$accession == "AY857599"] <- "ITS"
accs$final[accs$accession == "AY857556"] <- "trnL.trnF"
accs$final[accs$accession == "AY857600"] <- "ITS"
accs$final[accs$accession == "AY857558"] <- "trnL.trnF"
accs$final[accs$accession == "AY857601"] <- "ITS"
accs$final[accs$accession == "AY857602"] <- "ITS"
accs$final[accs$accession == "AY857560"] <- "trnL.trnF"
accs$final[accs$accession == "AY857603"] <- "ITS"
accs$final[accs$accession == "AY857604"] <- "ITS"
accs$final[accs$accession == "AY857562"] <- "trnL.trnF"
accs$final[accs$accession == "AY857605"] <- "ITS"
accs$final[accs$accession == "AY857606"] <- "ITS"
accs$final[accs$accession == "AY857564"] <- "trnL.trnF"
accs$final[accs$accession == "AY857608"] <- "ITS"
accs$final[accs$accession == "AY857566"] <- "trnL.trnF"
accs$final[accs$accession == "AY857609"] <- "ITS"
accs$final[accs$accession == "AY857567"] <- "trnL.trnF"
accs$final[accs$accession == "JN089244"] <- "rps7.atp6"
accs$final[accs$accession == "AY950403"] <- "trnL.trnF"
accs$final[accs$accession == "GU552288"] <- "trnL.trnF"
accs$final[accs$accession == "AY857615"] <- "ITS"
accs$final[accs$accession == "AY857572"] <- "trnL.trnF"
accs$final[accs$accession == "AY857610"] <- "ITS"
accs$final[accs$accession == "AY857573"] <- "trnL.trnF"
accs$final[accs$accession == "AY854393"] <- "ITS"
accs$final[accs$accession == "AY950406"] <- "trnL.trnF"
accs$final[accs$accession == "AY950418"] <- "trnL.trnF"
accs$final[accs$accession == "AY854409"] <- "ITS"
accs$final[accs$accession == "AY950420"] <- "trnL.trnF"
accs$final[accs$accession == "GU552281"] <- "trnL.trnF"
accs$final[accs$accession == "GU552291"] <- "trnL.trnF"
accs$final[accs$accession == "GU552289"] <- "trnL.trnF"
accs$final[accs$accession == "GU552292"] <- "trnL.trnF"
accs$final[accs$accession == "GU552287"] <- "trnL.trnF"
accs$final[accs$accession == "GU552280"] <- "trnL.trnF"
accs$final[accs$accession == "AF509843"] <- "ITS"
accs$final[accs$accession == "AY854410"] <- "ITS"
accs$final[accs$accession == "AY950422"] <- "trnL.trnF"
accs$final[accs$accession == "AY854412"] <- "ITS"
accs$final[accs$accession == "AY950423"] <- "trnL.trnF"
accs$final[accs$accession == "AY854413"] <- "ITS"
accs$final[accs$accession == "AY950424"] <- "trnL.trnF"
accs$final[accs$accession == "AY854415"] <- "ITS"
accs$final[accs$accession == "AY950426"] <- "trnL.trnF"
accs$final[accs$accession == "AY950429"] <- "trnL.trnF"
accs$final[accs$accession == "AY854420"] <- "ITS"
accs$final[accs$accession == "AY950430"] <- "trnL.trnF"
accs$final[accs$accession == "AY854421"] <- "ITS"
accs$final[accs$accession == "AY950431"] <- "trnL.trnF"
accs$final[accs$accession == "AY854423"] <- "ITS"
accs$final[accs$accession == "AY950435"] <- "trnL.trnF"
accs$final[accs$accession == "AY854425"] <- "ITS"
accs$final[accs$accession == "AY950436"] <- "trnL.trnF"
accs$final[accs$accession == "AY854426"] <- "ITS"
accs$final[accs$accession == "AY950437"] <- "trnL.trnF"
accs$final[accs$accession == "AY950438"] <- "trnL.trnF"
accs$final[accs$accession == "AY950440"] <- "trnL.trnF"
accs$final[accs$accession == "AY950428"] <- "trnL.trnF"
accs$final[accs$accession == "AY950350"] <- "atpB.rbcL"
accs$final[accs$accession == "AY950441"] <- "trnL.trnF"
accs$final[accs$accession == "AY854431"] <- "ITS"
accs$final[accs$accession == "AY950442"] <- "trnL.trnF"
accs$final[accs$accession == "AY854433"] <- "ITS"
accs$final[accs$accession == "AY950443"] <- "trnL.trnF"
accs <- accs %>% mutate(final = replace(final, final ==  "5S rRNA", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "ATP synthase beta", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "AtpB", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB.", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB-rbcL intergenic spacer", "atpB.rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB-rbcL intergenic spacer; may also contain rbcL", "atpB.rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB-rbcL non-coding spacer", "atpB.rbcL"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpF-atpH intergenic spacer; may contain AtpF and", "atpF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains 18S ribosomal RNA and internal", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains 5.8S ribosomal RNA and internal", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains AtpB and atpB-rbcL intergenic spacer", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains AtpB, intergenic spacer and RbcL", "atpB"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains photosystem II protein T (psbT), psbT-psbN", "psbT"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains PsbA (psbA) gene, psbA-trnH intergenic", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains PsbA and psbA-trnK intergenic spacer", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA_Leu (trnL) gene and trnL-trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu (trnL) and tRNA-Phe (trnF) genes,", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu (trnL) gene and trnL-trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu and tRNA-Phe", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu and trnL-trnF intergenic spacer", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu gene and trnL-trnF intergenic", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu, trnL-trnF intergenic spacer and", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Leu, trnL-trnF intergenic spacer,", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains tRNA-Thr (trnT) and tRNA-Leu (trnL) genes,", "trnT"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnG", "trnG"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnL intron, trnL-trnF intergenic spacer", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnL, intergenic spacer and trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "contains trnL, intergenic spacer and trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing chloroplast small ribosomal", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "genomic DNA containing NADH dehydrogenase subunit 5", "nad5"))
accs <- accs %>% mutate(final = replace(final, final ==  "internal transcribed spacer 1, ITS1", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "internal transcribed spacer 1; may also contain 18S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "internal transcribed spacer 2, ITS2", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "4	internal transcribed spacer I", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "Musci", "ITS"))
#accs <- accs %>% mutate(final = replace(final, final ==  "may contain small ribosomal protein subunit 4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "Nad4", "nad4"))
accs$final[accs$accession == "AF545013"] <- "trnL.trnF"
accs <- accs %>% mutate(final = replace(final, final ==  "PCR_primers=fwd_name: ITS1, rev_name: ITS4", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "photosystem II protein", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbA gene", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbA-trnH intergenic spacer", "psbA.trnH"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbA-trnH intergenic spacer; may also contain psbA", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbB operon", "psbB"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein L16", "rpl16"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein L16 G2 intron", "rpl16"))
accs <- accs %>% mutate(final = replace(final, final ==  "rpL16", "rpl16"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein S4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein subunit 4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein system 4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "Rps4", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "rps4, partial", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "tRNA-Ser", "trnS"))
accs <- accs %>% mutate(final = replace(final, final ==  "tRNA-Glu", "trnE"))
accs <- accs %>% mutate(final = replace(final, final ==  "tRNA-Thr", "trnT"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnE-UUC", "trnE"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnGUCC G2 intron", "trnG"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL (CAA)", "trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL (UAA)", "trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL-trnF region", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL-UUA", "trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL(uaa)", "trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnT-trnL intergenic spacer", "trnT.trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnT(UGU)-trnL(UAA) intergenic spacer", "trnT.trnL"))
accs <- accs %>% mutate(final = replace(final, final ==  "ribosomal protein L32", "rpl23"))

noP <-  accs %>% filter(!species == "Physcomitrella patens")
checkGenes <- count(noP, final)
#okay


#write.csv(accs, file="accessions_genes_cleaned.csv", row.names=F)
#more edits AND FIND WHOLE PLASTID/MITO GENOMES! ################
sapply(accs, class)
accs <- accs %>% mutate(chars = as.numeric(chars))
accs[79758,7] <- 105340 #NC_007945
accs[82521,7] <- 100725 #NC_024518
accs[82533,7] <- 109586 #NC_024523
accs[86147,7] <- 161162 #NC_004543
accs[86122,7] <- 122890 #AP005672
accs[86117,7] <- 120546 #NC_019628
accs[86113,7] <- 153208 #NC_020259
#accs <- accs[-c(79811, 79757, 86122, 81773, 82514, 82518, 82526),]
check <- accs %>% filter(is.na(chars))
accs <- accs %>% filter(!is.na(chars))

goddamnspaces <- accs %>% select(accession)
accs <- accs %>% mutate(accession = as.character(accession))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_007945 ", "NC_007945"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_013765 ", "NC_013765"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "AP005672 ", "AP005672"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_004543 ", "NC_004543"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_019628 ", "NC_019628"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_020259 ", "NC_020259"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_024523 ", "NC_024523"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "NC_024518 ", "NC_024518"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "FJ262471.1", "FJ262471.2"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "FJ262470.1", "FJ262470.2"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "AF413556.1", "AF413556.2"))
accs <- accs %>% mutate(accession = replace(accession, accession ==  "AF413547.1", "AF413547.2"))
accs <- accs %>% distinct()


genomes <- accs %>% filter(chars > 90000) %>% select(accession)
genomes$definition <- NA
for (i in 1:nrow(genomes)){
  print(i)
  query<-paste(genomes$accession[i])
  see <- entrez_fetch(db="nuccore", id=query, rettype="gb")
  pattern="DEFINITION(.*?)\nACCESSION"
  result <- regmatches(see,regexec(pattern,see))
  df <- data.frame(authors=result)
  genomes[i,2]<-paste(df[2,1])
  print(paste("values for accession",i,"calculated",sep=" "))
}

check <- entrez_fetch(db="nuccore", id="MN056355.1", rettype="gb")
check
pattern="DEFINITION(.*?)\nACCESSION"
regmatches(check,regexec(pattern,check))



genomes2 <- genomes
genomes2$mito <- NA
genomes2$chloro <- NA
genomes2$plastid <- NA
genomes2$other <- NA
genomes2$other2 <- NA
genomes2$final <- NA
for (i in 1:nrow(genomes2)){
  value1 <- "mitochondrion"
  value2 <- "chloroplast"
  value3 <- "plastid"
  value4 <- "genome"
  value5 <- "genomic"
  chars <- paste(genomes2$definition[i])
  genomes2[i,3] <- paste(grepl(value1, chars, fixed = TRUE))
  genomes2[i,4] <- paste(grepl(value2, chars, fixed = TRUE))
  genomes2[i,5] <- paste(grepl(value3, chars, fixed = TRUE))
  genomes2[i,6] <- paste(grepl(value4, chars, fixed = TRUE))
  genomes2[i,7] <- paste(grepl(value5, chars, fixed = TRUE))
  if(genomes2$mito[i]==TRUE){
    genomes2[i,8] <- paste("Mitochondrion")
  }
  if(genomes2$chloro[i]==TRUE){
    genomes2[i,8] <- paste("Chloroplast")
  }
  if(genomes2$plastid[i]==TRUE){
    genomes2[i,8] <- paste("Chloroplast")
  }
  if(genomes2$other[i]==TRUE){
    genomes2[i,8] <- paste("whole genome or genomic DNA")
  }
  if(genomes2$other2[i]==TRUE){
    genomes2[i,8] <- paste("whole genome or genomic DNA")
  }
}


for(i in 1:nrow(accs)){
  query<- accs$accession[i]
  if(query %in%  genomes2$accession)  {
    temp<- genomes2 %>% filter(accession==query)
    buck<-paste(temp$final)
    accs[i,6]<-paste(buck)
  }
}


accs$add<-"YES"
accs[82648,8]<-"NO"
accs[82649,8]<-"NO"
accs$accession<-as.character(accs$accession)

for(i in 79024:nrow(accs)){
  acc <- accs$add[i]
  bacc<-accs$accession[i]
  if (acc=="YES"){
    cacc<-paste(bacc, ".1", sep="")
    print(cacc)
    accs[i, 1] <- paste(as.character(cacc))
  }
}
accs$accession[79024] 

accs <- accs %>% select(!add)
#write.csv(accs, file="accessions_genes_cleaned.csv", row.names=F)

###########make sure that they  match with seqs!!!!
rm(list = ls())
library(dplyr)
accs <- read.csv("accessions_genes_cleaned.csv")
x <- read.csv("x3.csv")
accs$accession<-as.character(accs$accession)
x$X1 <- as.character(x$X1)
x$X2 <- as.character(x$X2)

accs2<-accs
for(i in 1:nrow(accs)){
  query<- accs$accession[i]
  if(query %in%  x$X1)  {
    temp<- x %>% filter(X1==query)
    buck<-paste(temp$X2)
    accs[i,1]<-paste(buck)
  }
}

write.csv(accs, file="accessions_genes_cleaned.csv", row.names=F)

##########################################################



#getSequences######################
###### make CSV's of each gene FROM ROSE!!! i.e. not incl. genes not in rose #######
rm(list = ls())
library(dplyr)
library(reshape2)
accs <- read.csv("accessions_genes_cleaned.csv")
accs <- accs %>% mutate(final = as.character(final))
accs <- accs %>% mutate(final = replace(final, final ==  "18S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "5.8s", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "nad4", "nad5.nad4"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnL", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnF", "trnL.trnF"))
accs <- accs %>% mutate(final = replace(final, final ==  "atpB-rbcL intergenic spacer", "atpB.rbcL")) #should have done this wayyy back but oops!
accs <- accs %>% mutate(final = replace(final, final ==  "atpB.rbcL", "atpB")) #see rose alignments!
accs <- accs %>% mutate(final = replace(final, final ==  "rps7", "atp6"))
accs <- accs %>% mutate(final = replace(final, final ==  "rps7.atp6", "atp6"))
accs <- accs %>% mutate(final = replace(final, final ==  "nad5", "nad5.nad4"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbA.trnH", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnH", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnH.psbA", "psbA"))
accs <- accs %>% mutate(final = replace(final, final ==  "psbT", "psbB"))
accs <- accs %>% mutate(final = replace(final, final ==  "rpl23", "rpl2"))
#rose fucked up and combined the two trnT's.... (trnT-UGU and trnT-GGU, which is fucking stupid. here' i'm calling these rps4 and psbD, respectively)
accs <- accs %>% mutate(final = replace(final, final ==  "trnT", "rps4"))
accs <- accs %>% mutate(final = replace(final, final ==  "trnT.psbD", "psbD"))
#get rid of species that are just "sp."
test <- accs[grep(" sp. ", accs$species),]
accs <- accs[! grepl(" sp. ", accs$species),]
#get rid of species that are "cf." a
#ctually jk.....
#test <- accs[grep(" cf. ", accs$species),]
#accs <- accs[! grepl(" cf. ", accs$species),]
accs <- accs[! grepl("anonymous", accs$final),]



noP <-  accs %>% filter(!species == "Physcomitrella patens")
checkGenes <- count(noP, final)
#roseGenes <- read.csv("roseAcc.csv") %>% select(gene) %>% distinct()
#thoseGenes <- checkGenes[ checkGenes$final %in% roseGenes$gene, ] 
#genes <- thoseGenes %>% select(final) %>% as.matrix()

#find species that are only represented by a single gene:
sp <- accs %>% group_by(species) %>% tally() %>% filter(n==1)
species <- left_join(sp, accs)
spGenes <- species %>% select(final) %>% group_by(final) %>% tally()


genes <- checkGenes %>% filter(n>38) %>% select(final) %>% filter(!final=="whole genome or genomic DNA") %>% filter(!final=="Chloroplast") %>% as.matrix()
listDF <- list()
for(i in genes) {
  df.tmp <- accs %>% filter(final == i)
  max <- df.tmp %>% group_by(species) %>% summarise(chars = max(chars))
  df.final <- inner_join(max, df.tmp)%>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
  listDF[[i]] <- df.final
}

rm(df.tmp, max, df.final)
setwd("~/seqs")
for(x in 1:length(listDF)){
  write.csv(listDF[[x]],file = paste(names(listDF[x]),'.csv',sep=""), row.names = F)
}

rm(list = ls())
temp = list.files(pattern="*.csv")

list <- lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
               read.csv)


######################
for(i in 1:length(list)){
  df.tmp <- list[[i]]
  if(nrow(df.tmp)<100){
    df.tmpA <- df.tmp[1:nrow(df.tmp),]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    write.dna(seqsA, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>100 & nrow(df.tmp)<200 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:nrow(df.tmp),]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    seqs <- c(seqsA, seqsB)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>300 & nrow(df.tmp)<400 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:nrow(df.tmp),]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>500 & nrow(df.tmp)<600 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:nrow(df.tmp),]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>600 & nrow(df.tmp)<700 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:nrow(df.tmp),]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>700 & nrow(df.tmp)<800 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:nrow(df.tmp),]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>1000 & nrow(df.tmp)<1100 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:799,]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    df.tmpI <- df.tmp[800:899,]
    df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
    queryI<-paste(df.tmpI$accession)
    seqsI <- read.GenBank(queryI, as.character = T)
    spp_namesI <- as.vector(df.tmpI2$species) 
    names(seqsI) <- spp_namesI
    rm(df.tmpI, df.tmpI2, queryI, spp_namesI)
    
    df.tmpJ <- df.tmp[900:999,]
    df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
    queryJ<-paste(df.tmpJ$accession)
    seqsJ <- read.GenBank(queryJ, as.character = T)
    spp_namesJ <- as.vector(df.tmpJ2$species) 
    names(seqsJ) <- spp_namesJ
    rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)
    
    df.tmpK <- df.tmp[1000:nrow(df.tmp),]
    df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
    queryK<-paste(df.tmpK$accession)
    seqsK <- read.GenBank(queryK, as.character = T)
    spp_namesK <- as.vector(df.tmpK2$species) 
    names(seqsK) <- spp_namesK
    rm(df.tmpK, df.tmpK2, queryK, spp_namesK)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>1200 & nrow(df.tmp)<1300 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:799,]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    df.tmpI <- df.tmp[800:899,]
    df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
    queryI<-paste(df.tmpI$accession)
    seqsI <- read.GenBank(queryI, as.character = T)
    spp_namesI <- as.vector(df.tmpI2$species) 
    names(seqsI) <- spp_namesI
    rm(df.tmpI, df.tmpI2, queryI, spp_namesI)
    
    df.tmpJ <- df.tmp[900:999,]
    df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
    queryJ<-paste(df.tmpJ$accession)
    seqsJ <- read.GenBank(queryJ, as.character = T)
    spp_namesJ <- as.vector(df.tmpJ2$species) 
    names(seqsJ) <- spp_namesJ
    rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)
    
    df.tmpK <- df.tmp[1000:1099,]
    df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
    queryK<-paste(df.tmpK$accession)
    seqsK <- read.GenBank(queryK, as.character = T)
    spp_namesK <- as.vector(df.tmpK2$species) 
    names(seqsK) <- spp_namesK
    rm(df.tmpK, df.tmpK2, queryK, spp_namesK)
    
    df.tmpL <- df.tmp[1100:1199,]
    df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
    queryL<-paste(df.tmpL$accession)
    seqsL <- read.GenBank(queryL, as.character = T)
    spp_namesL <- as.vector(df.tmpL2$species) 
    names(seqsL) <- spp_namesL
    rm(df.tmpL, df.tmpL2, queryL, spp_namesL)
    
    df.tmpM <- df.tmp[1200:nrow(df.tmp),]
    df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
    queryM<-paste(df.tmpM$accession)
    seqsM <- read.GenBank(queryM, as.character = T)
    spp_namesM <- as.vector(df.tmpM2$species) 
    names(seqsM) <- spp_namesM
    rm(df.tmpM, df.tmpM2, queryM, spp_namesM)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>1400 & nrow(df.tmp)<1500 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:799,]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    df.tmpI <- df.tmp[800:899,]
    df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
    queryI<-paste(df.tmpI$accession)
    seqsI <- read.GenBank(queryI, as.character = T)
    spp_namesI <- as.vector(df.tmpI2$species) 
    names(seqsI) <- spp_namesI
    rm(df.tmpI, df.tmpI2, queryI, spp_namesI)
    
    df.tmpJ <- df.tmp[900:999,]
    df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
    queryJ<-paste(df.tmpJ$accession)
    seqsJ <- read.GenBank(queryJ, as.character = T)
    spp_namesJ <- as.vector(df.tmpJ2$species) 
    names(seqsJ) <- spp_namesJ
    rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)
    
    df.tmpK <- df.tmp[1000:1099,]
    df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
    queryK<-paste(df.tmpK$accession)
    seqsK <- read.GenBank(queryK, as.character = T)
    spp_namesK <- as.vector(df.tmpK2$species) 
    names(seqsK) <- spp_namesK
    rm(df.tmpK, df.tmpK2, queryK, spp_namesK)
    
    df.tmpL <- df.tmp[1100:1199,]
    df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
    queryL<-paste(df.tmpL$accession)
    seqsL <- read.GenBank(queryL, as.character = T)
    spp_namesL <- as.vector(df.tmpL2$species) 
    names(seqsL) <- spp_namesL
    rm(df.tmpL, df.tmpL2, queryL, spp_namesL)
    
    df.tmpM <- df.tmp[1200:1299,]
    df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
    queryM<-paste(df.tmpM$accession)
    seqsM <- read.GenBank(queryM, as.character = T)
    spp_namesM <- as.vector(df.tmpM2$species) 
    names(seqsM) <- spp_namesM
    rm(df.tmpM, df.tmpM2, queryM, spp_namesM)
    
    df.tmpN <- df.tmp[1300:1399,]
    df.tmpN2 <- as.data.frame(sapply(df.tmpN, function(x) gsub(" ", "_", x)))
    queryN<-paste(df.tmpN$accession)
    seqsN <- read.GenBank(queryN, as.character = T)
    spp_namesN <- as.vector(df.tmpN2$species) 
    names(seqsN) <- spp_namesN
    rm(df.tmpN, df.tmpN2, queryN, spp_namesN)
    
    df.tmpO <- df.tmp[1400:nrow(df.tmp),]
    df.tmpO2 <- as.data.frame(sapply(df.tmpO, function(x) gsub(" ", "_", x)))
    queryO<-paste(df.tmpO$accession)
    seqsO <- read.GenBank(queryO, as.character = T)
    spp_namesO <- as.vector(df.tmpO2$species) 
    names(seqsO) <- spp_namesO
    rm(df.tmpO, df.tmpO2, queryO, spp_namesO)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM, seqsN, seqsO)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>2300 & nrow(df.tmp)<2400 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:799,]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    df.tmpI <- df.tmp[800:899,]
    df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
    queryI<-paste(df.tmpI$accession)
    seqsI <- read.GenBank(queryI, as.character = T)
    spp_namesI <- as.vector(df.tmpI2$species) 
    names(seqsI) <- spp_namesI
    rm(df.tmpI, df.tmpI2, queryI, spp_namesI)
    
    df.tmpJ <- df.tmp[900:999,]
    df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
    queryJ<-paste(df.tmpJ$accession)
    seqsJ <- read.GenBank(queryJ, as.character = T)
    spp_namesJ <- as.vector(df.tmpJ2$species) 
    names(seqsJ) <- spp_namesJ
    rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)
    
    df.tmpK <- df.tmp[1000:1099,]
    df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
    queryK<-paste(df.tmpK$accession)
    seqsK <- read.GenBank(queryK, as.character = T)
    spp_namesK <- as.vector(df.tmpK2$species) 
    names(seqsK) <- spp_namesK
    rm(df.tmpK, df.tmpK2, queryK, spp_namesK)
    
    df.tmpL <- df.tmp[1100:1199,]
    df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
    queryL<-paste(df.tmpL$accession)
    seqsL <- read.GenBank(queryL, as.character = T)
    spp_namesL <- as.vector(df.tmpL2$species) 
    names(seqsL) <- spp_namesL
    rm(df.tmpL, df.tmpL2, queryL, spp_namesL)
    
    df.tmpM <- df.tmp[1200:1299,]
    df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
    queryM<-paste(df.tmpM$accession)
    seqsM <- read.GenBank(queryM, as.character = T)
    spp_namesM <- as.vector(df.tmpM2$species) 
    names(seqsM) <- spp_namesM
    rm(df.tmpM, df.tmpM2, queryM, spp_namesM)
    
    df.tmpN <- df.tmp[1300:1399,]
    df.tmpN2 <- as.data.frame(sapply(df.tmpN, function(x) gsub(" ", "_", x)))
    queryN<-paste(df.tmpN$accession)
    seqsN <- read.GenBank(queryN, as.character = T)
    spp_namesN <- as.vector(df.tmpN2$species) 
    names(seqsN) <- spp_namesN
    rm(df.tmpN, df.tmpN2, queryN, spp_namesN)
    
    df.tmpO <- df.tmp[1400:1499,]
    df.tmpO2 <- as.data.frame(sapply(df.tmpO, function(x) gsub(" ", "_", x)))
    queryO<-paste(df.tmpO$accession)
    seqsO <- read.GenBank(queryO, as.character = T)
    spp_namesO <- as.vector(df.tmpO2$species) 
    names(seqsO) <- spp_namesO
    rm(df.tmpO, df.tmpO2, queryO, spp_namesO)
    
    df.tmpP <- df.tmp[1500:1599,]
    df.tmpP2 <- as.data.frame(sapply(df.tmpP, function(x) gsub(" ", "_", x)))
    queryP<-paste(df.tmpP$accession)
    seqsP <- read.GenBank(queryP, as.character = T)
    spp_namesP <- as.vector(df.tmpP2$species) 
    names(seqsP) <- spp_namesP
    rm(df.tmpP, df.tmpP2, queryP, spp_namesP)
    
    df.tmpQ <- df.tmp[1600:1699,]
    df.tmpQ2 <- as.data.frame(sapply(df.tmpQ, function(x) gsub(" ", "_", x)))
    queryQ<-paste(df.tmpQ$accession)
    seqsQ <- read.GenBank(queryQ, as.character = T)
    spp_namesQ <- as.vector(df.tmpQ2$species) 
    names(seqsQ) <- spp_namesQ
    rm(df.tmpQ, df.tmpQ2, queryQ, spp_namesQ)
    
    df.tmpR <- df.tmp[1700:1799,]
    df.tmpR2 <- as.data.frame(sapply(df.tmpR, function(x) gsub(" ", "_", x)))
    queryR<-paste(df.tmpR$accession)
    seqsR <- read.GenBank(queryR, as.character = T)
    spp_namesR <- as.vector(df.tmpR2$species) 
    names(seqsR) <- spp_namesR
    rm(df.tmpR, df.tmpR2, queryR, spp_namesR)
    
    df.tmpS <- df.tmp[1800:1899,]
    df.tmpS2 <- as.data.frame(sapply(df.tmpS, function(x) gsub(" ", "_", x)))
    queryS<-paste(df.tmpS$accession)
    seqsS <- read.GenBank(queryS, as.character = T)
    spp_namesS <- as.vector(df.tmpS2$species) 
    names(seqsS) <- spp_namesS
    rm(df.tmpS, df.tmpS2, queryS, spp_namesS)
    
    df.tmpT <- df.tmp[1900:1999,]
    df.tmpT2 <- as.data.frame(sapply(df.tmpT, function(x) gsub(" ", "_", x)))
    queryT<-paste(df.tmpT$accession)
    seqsT <- read.GenBank(queryT, as.character = T)
    spp_namesT <- as.vector(df.tmpT2$species) 
    names(seqsT) <- spp_namesT
    rm(df.tmpT, df.tmpT2, queryT, spp_namesT)
    
    df.tmpU <- df.tmp[2000:2099,]
    df.tmpU2 <- as.data.frame(sapply(df.tmpU, function(x) gsub(" ", "_", x)))
    queryU<-paste(df.tmpU$accession)
    seqsU <- read.GenBank(queryU, as.character = T)
    spp_namesU <- as.vector(df.tmpU2$species) 
    names(seqsU) <- spp_namesU
    rm(df.tmpU, df.tmpU2, queryU, spp_namesU)
    
    df.tmpV <- df.tmp[2100:2199,]
    df.tmpV2 <- as.data.frame(sapply(df.tmpV, function(x) gsub(" ", "_", x)))
    queryV<-paste(df.tmpV$accession)
    seqsV <- read.GenBank(queryV, as.character = T)
    spp_namesV <- as.vector(df.tmpV2$species) 
    names(seqsV) <- spp_namesV
    rm(df.tmpV, df.tmpV2, queryV, spp_namesV)
    
    df.tmpW <- df.tmp[2200:2299,]
    df.tmpW2 <- as.data.frame(sapply(df.tmpW, function(x) gsub(" ", "_", x)))
    queryW<-paste(df.tmpW$accession)
    seqsW <- read.GenBank(queryW, as.character = T)
    spp_namesW <- as.vector(df.tmpW2$species) 
    names(seqsW) <- spp_namesW
    rm(df.tmpW, df.tmpW2, queryW, spp_namesW)
    
    df.tmpX <- df.tmp[2300:nrow(df.tmp),]
    df.tmpX2 <- as.data.frame(sapply(df.tmpX, function(x) gsub(" ", "_", x)))
    queryX<-paste(df.tmpX$accession)
    seqsX <- read.GenBank(queryX, as.character = T)
    spp_namesX <- as.vector(df.tmpX2$species) 
    names(seqsX) <- spp_namesX
    rm(df.tmpX, df.tmpX2, queryX, spp_namesX)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM, seqsN, seqsO, seqsP, seqsQ, seqsR, seqsS, seqsT, seqsU, seqsV, seqsW, seqsX)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
  if(nrow(df.tmp)>2400 & nrow(df.tmp)<2500 ){
    df.tmpA <- df.tmp[1:99,]
    df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
    queryA<-paste(df.tmpA$accession)
    seqsA <- read.GenBank(queryA, as.character = T)
    spp_namesA <- as.vector(df.tmpA2$species) 
    names(seqsA) <- spp_namesA
    rm(df.tmpA, df.tmpA2, queryA, spp_namesA)
    
    df.tmpB <- df.tmp[100:199,]
    df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
    queryB<-paste(df.tmpB$accession)
    seqsB <- read.GenBank(queryB, as.character = T)
    spp_namesB <- as.vector(df.tmpB2$species) 
    names(seqsB) <- spp_namesB
    rm(df.tmpB, df.tmpB2, queryB, spp_namesB)
    
    df.tmpC <- df.tmp[200:299,]
    df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
    queryC<-paste(df.tmpC$accession)
    seqsC <- read.GenBank(queryC, as.character = T)
    spp_namesC <- as.vector(df.tmpC2$species) 
    names(seqsC) <- spp_namesC
    rm(df.tmpC, df.tmpC2, queryC, spp_namesC)
    
    df.tmpD <- df.tmp[300:399,]
    df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
    queryD<-paste(df.tmpD$accession)
    seqsD <- read.GenBank(queryD, as.character = T)
    spp_namesD <- as.vector(df.tmpD2$species) 
    names(seqsD) <- spp_namesD
    rm(df.tmpD, df.tmpD2, queryD, spp_namesD)
    
    df.tmpE <- df.tmp[400:499,]
    df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
    queryE<-paste(df.tmpE$accession)
    seqsE <- read.GenBank(queryE, as.character = T)
    spp_namesE <- as.vector(df.tmpE2$species) 
    names(seqsE) <- spp_namesE
    rm(df.tmpE, df.tmpE2, queryE, spp_namesE)
    
    df.tmpF <- df.tmp[500:599,]
    df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
    queryF<-paste(df.tmpF$accession)
    seqsF <- read.GenBank(queryF, as.character = T)
    spp_namesF <- as.vector(df.tmpF2$species) 
    names(seqsF) <- spp_namesF
    rm(df.tmpF, df.tmpF2, queryF, spp_namesF)
    
    df.tmpG <- df.tmp[600:699,]
    df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
    queryG<-paste(df.tmpG$accession)
    seqsG <- read.GenBank(queryG, as.character = T)
    spp_namesG <- as.vector(df.tmpG2$species) 
    names(seqsG) <- spp_namesG
    rm(df.tmpG, df.tmpG2, queryG, spp_namesG)
    
    df.tmpH <- df.tmp[700:799,]
    df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
    queryH<-paste(df.tmpH$accession)
    seqsH <- read.GenBank(queryH, as.character = T)
    spp_namesH <- as.vector(df.tmpH2$species) 
    names(seqsH) <- spp_namesH
    rm(df.tmpH, df.tmpH2, queryH, spp_namesH)
    
    df.tmpI <- df.tmp[800:899,]
    df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
    queryI<-paste(df.tmpI$accession)
    seqsI <- read.GenBank(queryI, as.character = T)
    spp_namesI <- as.vector(df.tmpI2$species) 
    names(seqsI) <- spp_namesI
    rm(df.tmpI, df.tmpI2, queryI, spp_namesI)
    
    df.tmpJ <- df.tmp[900:999,]
    df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
    queryJ<-paste(df.tmpJ$accession)
    seqsJ <- read.GenBank(queryJ, as.character = T)
    spp_namesJ <- as.vector(df.tmpJ2$species) 
    names(seqsJ) <- spp_namesJ
    rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)
    
    df.tmpK <- df.tmp[1000:1099,]
    df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
    queryK<-paste(df.tmpK$accession)
    seqsK <- read.GenBank(queryK, as.character = T)
    spp_namesK <- as.vector(df.tmpK2$species) 
    names(seqsK) <- spp_namesK
    rm(df.tmpK, df.tmpK2, queryK, spp_namesK)
    
    df.tmpL <- df.tmp[1100:1199,]
    df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
    queryL<-paste(df.tmpL$accession)
    seqsL <- read.GenBank(queryL, as.character = T)
    spp_namesL <- as.vector(df.tmpL2$species) 
    names(seqsL) <- spp_namesL
    rm(df.tmpL, df.tmpL2, queryL, spp_namesL)
    
    df.tmpM <- df.tmp[1200:1299,]
    df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
    queryM<-paste(df.tmpM$accession)
    seqsM <- read.GenBank(queryM, as.character = T)
    spp_namesM <- as.vector(df.tmpM2$species) 
    names(seqsM) <- spp_namesM
    rm(df.tmpM, df.tmpM2, queryM, spp_namesM)
    
    df.tmpN <- df.tmp[1300:1399,]
    df.tmpN2 <- as.data.frame(sapply(df.tmpN, function(x) gsub(" ", "_", x)))
    queryN<-paste(df.tmpN$accession)
    seqsN <- read.GenBank(queryN, as.character = T)
    spp_namesN <- as.vector(df.tmpN2$species) 
    names(seqsN) <- spp_namesN
    rm(df.tmpN, df.tmpN2, queryN, spp_namesN)
    
    df.tmpO <- df.tmp[1400:1499,]
    df.tmpO2 <- as.data.frame(sapply(df.tmpO, function(x) gsub(" ", "_", x)))
    queryO<-paste(df.tmpO$accession)
    seqsO <- read.GenBank(queryO, as.character = T)
    spp_namesO <- as.vector(df.tmpO2$species) 
    names(seqsO) <- spp_namesO
    rm(df.tmpO, df.tmpO2, queryO, spp_namesO)
    
    df.tmpP <- df.tmp[1500:1599,]
    df.tmpP2 <- as.data.frame(sapply(df.tmpP, function(x) gsub(" ", "_", x)))
    queryP<-paste(df.tmpP$accession)
    seqsP <- read.GenBank(queryP, as.character = T)
    spp_namesP <- as.vector(df.tmpP2$species) 
    names(seqsP) <- spp_namesP
    rm(df.tmpP, df.tmpP2, queryP, spp_namesP)
    
    df.tmpQ <- df.tmp[1600:1699,]
    df.tmpQ2 <- as.data.frame(sapply(df.tmpQ, function(x) gsub(" ", "_", x)))
    queryQ<-paste(df.tmpQ$accession)
    seqsQ <- read.GenBank(queryQ, as.character = T)
    spp_namesQ <- as.vector(df.tmpQ2$species) 
    names(seqsQ) <- spp_namesQ
    rm(df.tmpQ, df.tmpQ2, queryQ, spp_namesQ)
    
    df.tmpR <- df.tmp[1700:1799,]
    df.tmpR2 <- as.data.frame(sapply(df.tmpR, function(x) gsub(" ", "_", x)))
    queryR<-paste(df.tmpR$accession)
    seqsR <- read.GenBank(queryR, as.character = T)
    spp_namesR <- as.vector(df.tmpR2$species) 
    names(seqsR) <- spp_namesR
    rm(df.tmpR, df.tmpR2, queryR, spp_namesR)
    
    df.tmpS <- df.tmp[1800:1899,]
    df.tmpS2 <- as.data.frame(sapply(df.tmpS, function(x) gsub(" ", "_", x)))
    queryS<-paste(df.tmpS$accession)
    seqsS <- read.GenBank(queryS, as.character = T)
    spp_namesS <- as.vector(df.tmpS2$species) 
    names(seqsS) <- spp_namesS
    rm(df.tmpS, df.tmpS2, queryS, spp_namesS)
    
    df.tmpT <- df.tmp[1900:1999,]
    df.tmpT2 <- as.data.frame(sapply(df.tmpT, function(x) gsub(" ", "_", x)))
    queryT<-paste(df.tmpT$accession)
    seqsT <- read.GenBank(queryT, as.character = T)
    spp_namesT <- as.vector(df.tmpT2$species) 
    names(seqsT) <- spp_namesT
    rm(df.tmpT, df.tmpT2, queryT, spp_namesT)
    
    df.tmpU <- df.tmp[2000:2099,]
    df.tmpU2 <- as.data.frame(sapply(df.tmpU, function(x) gsub(" ", "_", x)))
    queryU<-paste(df.tmpU$accession)
    seqsU <- read.GenBank(queryU, as.character = T)
    spp_namesU <- as.vector(df.tmpU2$species) 
    names(seqsU) <- spp_namesU
    rm(df.tmpU, df.tmpU2, queryU, spp_namesU)
    
    df.tmpV <- df.tmp[2100:2199,]
    df.tmpV2 <- as.data.frame(sapply(df.tmpV, function(x) gsub(" ", "_", x)))
    queryV<-paste(df.tmpV$accession)
    seqsV <- read.GenBank(queryV, as.character = T)
    spp_namesV <- as.vector(df.tmpV2$species) 
    names(seqsV) <- spp_namesV
    rm(df.tmpV, df.tmpV2, queryV, spp_namesV)
    
    df.tmpW <- df.tmp[2200:2299,]
    df.tmpW2 <- as.data.frame(sapply(df.tmpW, function(x) gsub(" ", "_", x)))
    queryW<-paste(df.tmpW$accession)
    seqsW <- read.GenBank(queryW, as.character = T)
    spp_namesW <- as.vector(df.tmpW2$species) 
    names(seqsW) <- spp_namesW
    rm(df.tmpW, df.tmpW2, queryW, spp_namesW)
    
    df.tmpX <- df.tmp[2300:2399,]
    df.tmpX2 <- as.data.frame(sapply(df.tmpX, function(x) gsub(" ", "_", x)))
    queryX<-paste(df.tmpX$accession)
    seqsX <- read.GenBank(queryX, as.character = T)
    spp_namesX <- as.vector(df.tmpX2$species) 
    names(seqsX) <- spp_namesX
    rm(df.tmpX, df.tmpX2, queryX, spp_namesX)
    
    df.tmpY <- df.tmp[2400:nrow(df.tmp),]
    df.tmpY2 <- as.data.frame(sapply(df.tmpY, function(x) gsub(" ", "_", x)))
    queryY<-paste(df.tmpY$accession)
    seqsY <- read.GenBank(queryY, as.character = T)
    spp_namesY <- as.vector(df.tmpY2$species) 
    names(seqsY) <- spp_namesY
    rm(df.tmpY, df.tmpY2, queryY, spp_namesY)
    
    seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM, seqsN, seqsO, seqsP, seqsQ, seqsR, seqsS, seqsT, seqsU, seqsV, seqsW, seqsX, seqsY)
    write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
  }
}




###
for(i in 1:length(list)){
  print(nrow(list[[i]])) 
}

c <- c(53, 44, 1069, 59, 2422, 93, 1297, 65, 532, 394, 41, 1410, 623, 92, 104, 2449, 89, 663, 106, 2318, 41, 151, 39, 789)
sort(c)


#GET MITOS AND CHLOROS! ########
accs <- read.csv("accessions_genes_cleaned.csv") %>% filter(final=="Mitochondrion" | final== "Chloroplast")
genes <- c("Chloroplast", "Mitochondrion")
listDF <- list()
for(i in genes) {
  df.tmp <- accs %>% filter(final == i)
  max <- df.tmp %>% group_by(species) %>% summarise(chars = max(chars))
  df.final <- inner_join(max, df.tmp)%>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
  listDF[[i]] <- df.final
}
setwd("~/seqs/genomes")
for(x in 1:length(listDF)){
  write.csv(listDF[[x]],file = paste(names(listDF[x]),'.csv',sep=""), row.names = F)
}

rm(list = ls())
m <- read.csv("Mitochondrion.csv")
m2 <- as.data.frame(sapply(m, function(x) gsub(" ", "_", x)))
seqs <- read.GenBank(m$accession, as.character = T)
spp_names <- as.vector(m2$species) 
names(seqs) <- spp_names
write.dna(seqs, "Mitochondrion.fasta", format="fasta")
mito <- read.dna("Mitochondrion.fasta", format="fasta")
stringset <- mito %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(stringset)
bin <- as.DNAbin(seqs)
write.dna(seqs, "Mitochondrion.fasta", format="fasta")


c <- read.csv("Chloroplast.csv")
c2 <- as.data.frame(sapply(c, function(x) gsub(" ", "_", x)))
seqs <- read.GenBank(c$accession, as.character = T)
spp_names <- as.vector(c2$species) 
names(seqs) <- spp_names
write.dna(seqs, "Chloroplast.fasta", format="fasta")
chloro <- read.dna("Chloroplast.fasta", format="fasta")
stringset <- chloro %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(stringset)
bin <- as.DNAbin(seqs)
write.dna(seqs, "Chloroplast.fasta", format="fasta")

#try <- read.GenBank(query, as.character = T)
#write.dna(try, "try.fasta", format="fasta")

#REDOING ITS: ##########
setwd("~/seqs")
rm(list = ls())
accs <- read.csv("accessions_genes_cleaned.csv")
accs <- accs %>% mutate(final = as.character(final))
accs <- accs[! grepl(" sp. ", accs$species),]
#ITS proper: ######
ITS.tmp <- accs %>% filter(final=="ITS")
max <- ITS.tmp %>% group_by(species) %>% summarise(chars = max(chars))
ITS <- inner_join(max, ITS.tmp) %>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
write.csv(ITS, "ITS.csv", row.names=F)
rm(max, ITS.tmp)

df.tmp <- ITS
df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:399,]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

df.tmpE <- df.tmp[400:499,]
df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
queryE<-paste(df.tmpE$accession)
seqsE <- read.GenBank(queryE, as.character = T)
spp_namesE <- as.vector(df.tmpE2$species) 
names(seqsE) <- spp_namesE
rm(df.tmpE, df.tmpE2, queryE, spp_namesE)

df.tmpF <- df.tmp[500:599,]
df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
queryF<-paste(df.tmpF$accession)
seqsF <- read.GenBank(queryF, as.character = T)
spp_namesF <- as.vector(df.tmpF2$species) 
names(seqsF) <- spp_namesF
rm(df.tmpF, df.tmpF2, queryF, spp_namesF)

df.tmpG <- df.tmp[600:699,]
df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
queryG<-paste(df.tmpG$accession)
seqsG <- read.GenBank(queryG, as.character = T)
spp_namesG <- as.vector(df.tmpG2$species) 
names(seqsG) <- spp_namesG
rm(df.tmpG, df.tmpG2, queryG, spp_namesG)

df.tmpH <- df.tmp[700:799,]
df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
queryH<-paste(df.tmpH$accession)
seqsH <- read.GenBank(queryH, as.character = T)
spp_namesH <- as.vector(df.tmpH2$species) 
names(seqsH) <- spp_namesH
rm(df.tmpH, df.tmpH2, queryH, spp_namesH)

df.tmpI <- df.tmp[800:899,]
df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
queryI<-paste(df.tmpI$accession)
seqsI <- read.GenBank(queryI, as.character = T)
spp_namesI <- as.vector(df.tmpI2$species) 
names(seqsI) <- spp_namesI
rm(df.tmpI, df.tmpI2, queryI, spp_namesI)

df.tmpJ <- df.tmp[900:999,]
df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
queryJ<-paste(df.tmpJ$accession)
seqsJ <- read.GenBank(queryJ, as.character = T)
spp_namesJ <- as.vector(df.tmpJ2$species) 
names(seqsJ) <- spp_namesJ
rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)

df.tmpK <- df.tmp[1000:1099,]
df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
queryK<-paste(df.tmpK$accession)
seqsK <- read.GenBank(queryK, as.character = T)
spp_namesK <- as.vector(df.tmpK2$species) 
names(seqsK) <- spp_namesK
rm(df.tmpK, df.tmpK2, queryK, spp_namesK)

df.tmpL <- df.tmp[1100:1199,]
df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
queryL<-paste(df.tmpL$accession)
seqsL <- read.GenBank(queryL, as.character = T)
spp_namesL <- as.vector(df.tmpL2$species) 
names(seqsL) <- spp_namesL
rm(df.tmpL, df.tmpL2, queryL, spp_namesL)

df.tmpM <- df.tmp[1200:1299,]
df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
queryM<-paste(df.tmpM$accession)
seqsM <- read.GenBank(queryM, as.character = T)
spp_namesM <- as.vector(df.tmpM2$species) 
names(seqsM) <- spp_namesM
rm(df.tmpM, df.tmpM2, queryM, spp_namesM)

df.tmpN <- df.tmp[1300:1399,]
df.tmpN2 <- as.data.frame(sapply(df.tmpN, function(x) gsub(" ", "_", x)))
queryN<-paste(df.tmpN$accession)
seqsN <- read.GenBank(queryN, as.character = T)
spp_namesN <- as.vector(df.tmpN2$species) 
names(seqsN) <- spp_namesN
rm(df.tmpN, df.tmpN2, queryN, spp_namesN)

df.tmpO <- df.tmp[1400:1499,]
df.tmpO2 <- as.data.frame(sapply(df.tmpO, function(x) gsub(" ", "_", x)))
queryO<-paste(df.tmpO$accession)
seqsO <- read.GenBank(queryO, as.character = T)
spp_namesO <- as.vector(df.tmpO2$species) 
names(seqsO) <- spp_namesO
rm(df.tmpO, df.tmpO2, queryO, spp_namesO)

df.tmpP <- df.tmp[1500:nrow(df.tmp),]
df.tmpP2 <- as.data.frame(sapply(df.tmpP, function(x) gsub(" ", "_", x)))
queryP<-paste(df.tmpP$accession)
seqsP <- read.GenBank(queryP, as.character = T)
spp_namesP <- as.vector(df.tmpP2$species) 
names(seqsP) <- spp_namesP
rm(df.tmpP, df.tmpP2, queryP, spp_namesP)

seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM, seqsN, seqsO, seqsP)
write.dna(seqs, file= "~/seqs/fastas/its.n.friends/ITS.fasta", format="fasta")


setwd("~/seqs/fastas/its.n.friends")
rm(list = ls())
ITSBIN <- read.dna("ITS.fasta", format="fasta") 
ITS <- ITSBIN %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(ITS)
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, file="ITS_aln.fasta")
#5.8S #######
x5.8s.tmp <- accs %>% filter(final=="5.8s")
max <- x5.8s.tmp %>% group_by(species) %>% summarise(chars = max(chars))
x5.8s <- inner_join(max, x5.8s.tmp) %>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
write.csv(x5.8s, "x5.8s.csv", row.names=F)
rm(max, x5.8s.tmp)

df.tmp <- x5.8s
df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:nrow(df.tmp),]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

seqs <- c(seqsA, seqsB, seqsC, seqsD)
write.dna(seqs, file= "~/seqs/fastas/its.n.friends/x5.8s.fasta", format="fasta")


rm(list = ls())
x5.8sBIN <- read.dna("x5.8s.fasta", format="fasta") 
x5.8s <- x5.8sBIN %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(x5.8s)
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, file="x5.8s_aln.fasta")

#26S ######

x26s <- read.csv("X26S.csv")
df.tmp <- x26s
df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:399,]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

df.tmpE <- df.tmp[400:499,]
df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
queryE<-paste(df.tmpE$accession)
seqsE <- read.GenBank(queryE, as.character = T)
spp_namesE <- as.vector(df.tmpE2$species) 
names(seqsE) <- spp_namesE
rm(df.tmpE, df.tmpE2, queryE, spp_namesE)

df.tmpF <- df.tmp[500:599,]
df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
queryF<-paste(df.tmpF$accession)
seqsF <- read.GenBank(queryF, as.character = T)
spp_namesF <- as.vector(df.tmpF2$species) 
names(seqsF) <- spp_namesF
rm(df.tmpF, df.tmpF2, queryF, spp_namesF)

df.tmpG <- df.tmp[600:699,]
df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
queryG<-paste(df.tmpG$accession)
seqsG <- read.GenBank(queryG, as.character = T)
spp_namesG <- as.vector(df.tmpG2$species) 
names(seqsG) <- spp_namesG
rm(df.tmpG, df.tmpG2, queryG, spp_namesG)

df.tmpH <- df.tmp[700:nrow(df.tmp),]
df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
queryH<-paste(df.tmpH$accession)
seqsH <- read.GenBank(queryH, as.character = T)
spp_namesH <- as.vector(df.tmpH2$species) 
names(seqsH) <- spp_namesH
rm(df.tmpH, df.tmpH2, queryH, spp_namesH)

seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH)
write.dna(seqs, file= "~/seqs/fastas/its.n.friends/x26s.fasta", format="fasta")


rm(list = ls())
x26sBIN <- read.dna("x26s.fasta", format="fasta") 
x26s <- x26sBIN %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(x26s)
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, file="x26s_aln.fasta")
bin <- as.DNAbin(aligned)
library(Biostrings)


pid(p, type="PID1")

one <- toString(x26s[1])
two <- toString(x26s[2])
p <- pairwiseAlignment(one,two)

#18S#######
x18s.tmp <- accs %>% filter(final=="18S")
max <- x18s.tmp %>% group_by(species) %>% summarise(chars = max(chars))
x18s <- inner_join(max, x18s.tmp) %>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
write.csv(x18s, "x18s.csv", row.names=F)
rm(max, x18s.tmp)
df.tmp <- x18s

df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:399,]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

df.tmpE <- df.tmp[400:499,]
df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
queryE<-paste(df.tmpE$accession)
seqsE <- read.GenBank(queryE, as.character = T)
spp_namesE <- as.vector(df.tmpE2$species) 
names(seqsE) <- spp_namesE
rm(df.tmpE, df.tmpE2, queryE, spp_namesE)

df.tmpF <- df.tmp[500:599,]
df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
queryF<-paste(df.tmpF$accession)
seqsF <- read.GenBank(queryF, as.character = T)
spp_namesF <- as.vector(df.tmpF2$species) 
names(seqsF) <- spp_namesF
rm(df.tmpF, df.tmpF2, queryF, spp_namesF)

df.tmpG <- df.tmp[600:699,]
df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
queryG<-paste(df.tmpG$accession)
seqsG <- read.GenBank(queryG, as.character = T)
spp_namesG <- as.vector(df.tmpG2$species) 
names(seqsG) <- spp_namesG
rm(df.tmpG, df.tmpG2, queryG, spp_namesG)

df.tmpH <- df.tmp[700:799,]
df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
queryH<-paste(df.tmpH$accession)
seqsH <- read.GenBank(queryH, as.character = T)
spp_namesH <- as.vector(df.tmpH2$species) 
names(seqsH) <- spp_namesH
rm(df.tmpH, df.tmpH2, queryH, spp_namesH)

df.tmpI <- df.tmp[800:nrow(df.tmp),]
df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
queryI<-paste(df.tmpI$accession)
seqsI <- read.GenBank(queryI, as.character = T)
spp_namesI <- as.vector(df.tmpI2$species) 
names(seqsI) <- spp_namesI
rm(df.tmpI, df.tmpI2, queryI, spp_namesI)

seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI)
write.dna(seqs, file= "~/seqs/fastas/its.n.friends/x18s.fasta", format="fasta")


rm(list = ls())
x18sBIN <- read.dna("x18s.fasta", format="fasta") 
x18s <- x18sBIN %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(x18s)
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, file="x18s_aln.fasta")
###########

# align #######
setwd("~/seqs/fastas")
rm(list = ls())
library(DECIPHER)
library(magrittr)
#get fastas into env
temp = list.files(pattern="*.fasta")
temp
x <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
            read.dna, format = "fasta")
list2env(
  lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
         read.dna, format = "fasta"), envir = .GlobalEnv)

setwd("~/seqs/fastas/align")
nam <- names(x)
for(i in nam) {
  df.tmp <- get(i)
  stringset <- df.tmp %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  seqs <- OrientNucleotides(stringset)
  aligned <- AlignSeqs(seqs)
  bin <- as.DNAbin(aligned)
  write.dna(bin, file= paste0(i, "_align.fasta", sep=""), format="fasta")
}


# actually i'm going to just make sure that they are in the same order..... #####
setwd("~/add.to.rose/seqs")
rm(list = ls())
library(DECIPHER)
library(magrittr)
#get fastas into env
temp = list.files(pattern="*.fasta")
x <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
            read.dna, format = "fasta")
list2env(
  lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
         read.dna, format = "fasta"), envir = .GlobalEnv)

nam <- names(x)
for(i in nam) {
  df.tmp <- get(i)
  stringset <- df.tmp %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
  seqs <- OrientNucleotides(stringset)
  bin <- as.DNAbin(seqs)
  write.dna(bin, file= paste0(i, ".fasta", sep=""), format="fasta")
}



###### write fastas!!!! this also RENAMES and################
setwd("~/seqs/fastas")
for(i in 1:length(list)){
  df.tmp <- list[[i]]
  df.tmp2 <- as.data.frame(sapply(df.tmp, function(x) gsub(" ", "_", x)))
  seqs <- seq[c(which(names(seq) %in% df.tmp$accession))]
  spp_names <- as.vector(df.tmp2$species) 
  names(seqs) <- spp_names
  write.dna(seqs, file= paste(names(list[i]), ".fasta", sep=""), format="fasta")
}

###FIXING RPL16! #####
rm(list = ls())
library(dplyr)
library(reshape2)
accs <- read.csv("accessions_genes_cleaned.csv")
rpl16 <- accs %>% filter(final=="rpl16")
rpl16Acc <- rpl16 %>% select(accession)
rpl16Acc$definition <- NA

for (i in 1274:nrow(rpl16Acc)){
  print(i)
  query<-paste(rpl16Acc$accession[i])
  see <- entrez_fetch(db="nuccore", id=query, rettype="gb")
  pattern="DEFINITION(.*?)\nACCESSION"
  result <- regmatches(see,regexec(pattern,see))
  df <- data.frame(authors=result)
  rpl16Acc[i,2]<-paste(df[2,1])
  print(paste("values for accession",i,"calculated",sep=" "))
}

for (i in 1:nrow(rpl16Acc)){
  print(i)
  query<-paste(rpl16Acc$accession[i])
  see <- entrez_fetch(db="nuccore", id=query, rettype="gb")
  organelle="/organelle=\"(.*?)\n"
  resultOrg <- regmatches(see,regexec(organelle,see))
  dforg <- data.frame(authors=resultOrg)
  rpl16Acc[i,7]<-paste(dforg[2,1])
  print(paste("values for accession",i,"calculated",sep=" "))
}



see <- entrez_fetch(db="nuccore", id="HG009165.1", rettype="gb")
see
pattern="/organelle=\"(.*?)\n"
regmatches(see,regexec(pattern,see))

safe <- rpl16Acc
rpl16Acc <- safe

rpl16Acc$mito <- NA
rpl16Acc$chloro <- NA
rpl16Acc$plastid <- NA
rpl16Acc$final <- NA
for (i in 1:nrow(rpl16Acc)){
  value1 <- "mitochondrion"
  value2 <- "chloroplast"
  value3 <- "plastid"
  chars <- paste(rpl16Acc$definition[i])
  rpl16Acc[i,3] <- paste(grepl(value1, chars, fixed = TRUE))
  rpl16Acc[i,4] <- paste(grepl(value2, chars, fixed = TRUE))
  rpl16Acc[i,5] <- paste(grepl(value3, chars, fixed = TRUE))
  if(rpl16Acc$mito[i]==TRUE){
    rpl16Acc[i,6] <- paste("Mitochondrion")
  }
  if(rpl16Acc$chloro[i]==TRUE){
    rpl16Acc[i,6] <- paste("Chloroplast")
  }
  if(rpl16Acc$plastid[i]==TRUE){
    rpl16Acc[i,6] <- paste("Chloroplast")
  }
}

rpl16 <- rpl16Acc %>% filter(!is.na(final)) %>% select(accession)
sapply(rpl16, class)
sapply(accs, class)
fin <- left_join(rpl16, accs)
head(fin)

max <- fin %>% group_by(species) %>% summarise(chars = max(chars))
df.final <- inner_join(max, fin)%>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
write.csv(df.final, "~/seqs/rpl16.csv", row.names=F)

df.tmp <- df.final
df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:399,]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

df.tmpE <- df.tmp[400:499,]
df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
queryE<-paste(df.tmpE$accession)
seqsE <- read.GenBank(queryE, as.character = T)
spp_namesE <- as.vector(df.tmpE2$species) 
names(seqsE) <- spp_namesE
rm(df.tmpE, df.tmpE2, queryE, spp_namesE)

df.tmpF <- df.tmp[500:nrow(df.tmp),]
df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
queryF<-paste(df.tmpF$accession)
seqsF <- read.GenBank(queryF, as.character = T)
spp_namesF <- as.vector(df.tmpF2$species) 
names(seqsF) <- spp_namesF
rm(df.tmpF, df.tmpF2, queryF, spp_namesF)

seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF)
write.dna(seqs, file="~/seqs/fastas/rpl16.fasta", format="fasta")

setwd("~/seqs/fastas/align")
rm(list = ls())
rpl16BIN <- read.dna("~/seqs/fastas/rpl16.fasta", format="fasta") 
rpl16 <- rpl16BIN %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
seqs <- OrientNucleotides(rpl16)
aligned <- AlignSeqs(seqs)
writeXStringSet(aligned, file="rpl16_align.fasta")


###REDO ITS ######
rm(list = ls())
library(dplyr)
library(reshape2)
accs <- read.csv("accessions_genes_cleaned.csv")
accs <- accs %>% mutate(final = as.character(final))
accs <- accs %>% mutate(final = replace(final, final ==  "18S", "ITS"))
accs <- accs %>% mutate(final = replace(final, final ==  "5.8s", "ITS"))
#get rid of species that are just "sp."
test <- accs[grep(" sp. ", accs$species),]
accs <- accs[! grepl(" sp. ", accs$species),]

#ITS proper: ######
ITS.tmp <- accs %>% filter(final=="ITS")
max <- ITS.tmp %>% group_by(species) %>% summarise(chars = max(chars))
ITS <- inner_join(max, ITS.tmp) %>% distinct(species, chars, final, .keep_all = T) %>% mutate(species = as.character(species))
write.csv(ITS, "ITS.csv", row.names=F)
rm(max, ITS.tmp)

df.tmp <- ITS
df.tmpA <- df.tmp[1:99,]
df.tmpA2 <- as.data.frame(sapply(df.tmpA, function(x) gsub(" ", "_", x)))
queryA<-paste(df.tmpA$accession)
seqsA <- read.GenBank(queryA, as.character = T)
spp_namesA <- as.vector(df.tmpA2$species) 
names(seqsA) <- spp_namesA
rm(df.tmpA, df.tmpA2, queryA, spp_namesA)

df.tmpB <- df.tmp[100:199,]
df.tmpB2 <- as.data.frame(sapply(df.tmpB, function(x) gsub(" ", "_", x)))
queryB<-paste(df.tmpB$accession)
seqsB <- read.GenBank(queryB, as.character = T)
spp_namesB <- as.vector(df.tmpB2$species) 
names(seqsB) <- spp_namesB
rm(df.tmpB, df.tmpB2, queryB, spp_namesB)

df.tmpC <- df.tmp[200:299,]
df.tmpC2 <- as.data.frame(sapply(df.tmpC, function(x) gsub(" ", "_", x)))
queryC<-paste(df.tmpC$accession)
seqsC <- read.GenBank(queryC, as.character = T)
spp_namesC <- as.vector(df.tmpC2$species) 
names(seqsC) <- spp_namesC
rm(df.tmpC, df.tmpC2, queryC, spp_namesC)

df.tmpD <- df.tmp[300:399,]
df.tmpD2 <- as.data.frame(sapply(df.tmpD, function(x) gsub(" ", "_", x)))
queryD<-paste(df.tmpD$accession)
seqsD <- read.GenBank(queryD, as.character = T)
spp_namesD <- as.vector(df.tmpD2$species) 
names(seqsD) <- spp_namesD
rm(df.tmpD, df.tmpD2, queryD, spp_namesD)

df.tmpE <- df.tmp[400:499,]
df.tmpE2 <- as.data.frame(sapply(df.tmpE, function(x) gsub(" ", "_", x)))
queryE<-paste(df.tmpE$accession)
seqsE <- read.GenBank(queryE, as.character = T)
spp_namesE <- as.vector(df.tmpE2$species) 
names(seqsE) <- spp_namesE
rm(df.tmpE, df.tmpE2, queryE, spp_namesE)

df.tmpF <- df.tmp[500:599,]
df.tmpF2 <- as.data.frame(sapply(df.tmpF, function(x) gsub(" ", "_", x)))
queryF<-paste(df.tmpF$accession)
seqsF <- read.GenBank(queryF, as.character = T)
spp_namesF <- as.vector(df.tmpF2$species) 
names(seqsF) <- spp_namesF
rm(df.tmpF, df.tmpF2, queryF, spp_namesF)

df.tmpG <- df.tmp[600:699,]
df.tmpG2 <- as.data.frame(sapply(df.tmpG, function(x) gsub(" ", "_", x)))
queryG<-paste(df.tmpG$accession)
seqsG <- read.GenBank(queryG, as.character = T)
spp_namesG <- as.vector(df.tmpG2$species) 
names(seqsG) <- spp_namesG
rm(df.tmpG, df.tmpG2, queryG, spp_namesG)

df.tmpH <- df.tmp[700:799,]
df.tmpH2 <- as.data.frame(sapply(df.tmpH, function(x) gsub(" ", "_", x)))
queryH<-paste(df.tmpH$accession)
seqsH <- read.GenBank(queryH, as.character = T)
spp_namesH <- as.vector(df.tmpH2$species) 
names(seqsH) <- spp_namesH
rm(df.tmpH, df.tmpH2, queryH, spp_namesH)

df.tmpI <- df.tmp[800:899,]
df.tmpI2 <- as.data.frame(sapply(df.tmpI, function(x) gsub(" ", "_", x)))
queryI<-paste(df.tmpI$accession)
seqsI <- read.GenBank(queryI, as.character = T)
spp_namesI <- as.vector(df.tmpI2$species) 
names(seqsI) <- spp_namesI
rm(df.tmpI, df.tmpI2, queryI, spp_namesI)

df.tmpJ <- df.tmp[900:999,]
df.tmpJ2 <- as.data.frame(sapply(df.tmpJ, function(x) gsub(" ", "_", x)))
queryJ<-paste(df.tmpJ$accession)
seqsJ <- read.GenBank(queryJ, as.character = T)
spp_namesJ <- as.vector(df.tmpJ2$species) 
names(seqsJ) <- spp_namesJ
rm(df.tmpJ, df.tmpJ2, queryJ, spp_namesJ)

df.tmpK <- df.tmp[1000:1099,]
df.tmpK2 <- as.data.frame(sapply(df.tmpK, function(x) gsub(" ", "_", x)))
queryK<-paste(df.tmpK$accession)
seqsK <- read.GenBank(queryK, as.character = T)
spp_namesK <- as.vector(df.tmpK2$species) 
names(seqsK) <- spp_namesK
rm(df.tmpK, df.tmpK2, queryK, spp_namesK)

df.tmpL <- df.tmp[1100:1199,]
df.tmpL2 <- as.data.frame(sapply(df.tmpL, function(x) gsub(" ", "_", x)))
queryL<-paste(df.tmpL$accession)
seqsL <- read.GenBank(queryL, as.character = T)
spp_namesL <- as.vector(df.tmpL2$species) 
names(seqsL) <- spp_namesL
rm(df.tmpL, df.tmpL2, queryL, spp_namesL)

df.tmpM <- df.tmp[1200:1299,]
df.tmpM2 <- as.data.frame(sapply(df.tmpM, function(x) gsub(" ", "_", x)))
queryM<-paste(df.tmpM$accession)
seqsM <- read.GenBank(queryM, as.character = T)
spp_namesM <- as.vector(df.tmpM2$species) 
names(seqsM) <- spp_namesM
rm(df.tmpM, df.tmpM2, queryM, spp_namesM)

df.tmpN <- df.tmp[1300:1399,]
df.tmpN2 <- as.data.frame(sapply(df.tmpN, function(x) gsub(" ", "_", x)))
queryN<-paste(df.tmpN$accession)
seqsN <- read.GenBank(queryN, as.character = T)
spp_namesN <- as.vector(df.tmpN2$species) 
names(seqsN) <- spp_namesN
rm(df.tmpN, df.tmpN2, queryN, spp_namesN)

df.tmpO <- df.tmp[1400:1499,]
df.tmpO2 <- as.data.frame(sapply(df.tmpO, function(x) gsub(" ", "_", x)))
queryO<-paste(df.tmpO$accession)
seqsO <- read.GenBank(queryO, as.character = T)
spp_namesO <- as.vector(df.tmpO2$species) 
names(seqsO) <- spp_namesO
rm(df.tmpO, df.tmpO2, queryO, spp_namesO)

df.tmpP <- df.tmp[1500:nrow(df.tmp),]
df.tmpP2 <- as.data.frame(sapply(df.tmpP, function(x) gsub(" ", "_", x)))
queryP<-paste(df.tmpP$accession)
seqsP <- read.GenBank(queryP, as.character = T)
spp_namesP <- as.vector(df.tmpP2$species) 
names(seqsP) <- spp_namesP
rm(df.tmpP, df.tmpP2, queryP, spp_namesP)

seqs <- c(seqsA, seqsB, seqsC, seqsD, seqsE, seqsF, seqsG, seqsH, seqsI, seqsJ, seqsK, seqsL, seqsM, seqsN, seqsO, seqsP)
write.dna(seqs, file= "~/seqs/fastas/its.n.friends/ITS.fasta", format="fasta")





rm(list = ls())
library(dplyr)

####MITO ######
setwd("~/seqs/mito")
temp = list.files(pattern="*.csv")
list <- lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
               read.csv)
setwd("~/seqs/done/mito")
temp = list.files(pattern="*.fasta")
list2 <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta")

mito.bad <- read.csv("~/seqs/Mitochondrion.csv") %>% mutate(species = as.character(species))  %>% select(species, accession) 
mito <- as.data.frame(sapply(mito.bad, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species))
mito$oldName <- mito.bad$species
sapply(mito, class)

df_list <- list()
names <- names(list)
for(i in 1:length(list)){
  df.tmp.bad <- list[[i]]  
  df.tmp <- as.data.frame(sapply(df.tmp.bad, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species)) %>% select(species, accession)
  df.tmp$oldName <- df.tmp.bad$species
  fasta.tmp <- list2[[i]]
  sp <- data.frame(species = labels(fasta.tmp))
  sp2 <- as.data.frame(sapply(sp, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species))
  species <- as.data.frame(sapply(sp2, function(x) gsub("_", " ", x))) %>% mutate(species = as.character(species))
  sapply(species, class)
  mito_in <- inner_join(mito, species, by = "species")
  not_in <- anti_join(species, mito_in, by = "species")
  not_accessions <- left_join(not_in, df.tmp, by = "species") 
  joined <- rbind(mito_in, not_accessions)
  joined$gene <- names[i]
  df_list[[i]] <- joined
}

big_mito = do.call(rbind, df_list)


####chloro ######
setwd("~/seqs/plastid")
temp = list.files(pattern="*.csv")
list <- lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
               read.csv)
setwd("~/seqs/done/plastid")
temp = list.files(pattern="*.fasta")
list2 <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta")

chloro.bad <- read.csv("~/seqs/Chloroplast.csv") %>% mutate(species = as.character(species))  %>% select(species, accession) 
chloro <- as.data.frame(sapply(chloro.bad, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species))
chloro$oldName <- chloro.bad$species
sapply(chloro, class)

df_list <- list()
names <- names(list)
for(i in 1:length(list)){
  df.tmp.bad <- list[[i]]  
  df.tmp.bad2 <- as.data.frame(sapply(df.tmp.bad, function(x) gsub("\\[", "", x)))
  df.tmp.bad3 <- as.data.frame(sapply(df.tmp.bad2, function(x) gsub("\\]", "", x)))
  df.tmp <- as.data.frame(sapply(df.tmp.bad3, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species)) %>% select(species, accession)
  df.tmp$oldName <- df.tmp.bad$species
  fasta.tmp <- list2[[i]]
  sp <- data.frame(species = labels(fasta.tmp))
  sp2 <- as.data.frame(sapply(sp, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species))
  species <- as.data.frame(sapply(sp2, function(x) gsub("_", " ", x))) %>% mutate(species = as.character(species))
  sapply(species, class)
  chloro_in <- inner_join(chloro, species, by = "species")
  not_in <- anti_join(species, chloro_in, by = "species")
  not_accessions <- left_join(not_in, df.tmp, by = "species") 
  joined <- rbind(chloro_in, not_accessions)
  joined$gene <- names[i]
  df_list[[i]] <- joined
}

big_chloro = do.call(rbind, df_list)

big_compartments <- rbind(big_mito, big_chloro)






####FINALLY, ITS ######
rm(list = ls())
library(ape)
library(dplyr)
setwd("~/seqs/ITS")
temp = list.files(pattern="*.csv")
list <- lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
               read.csv)
setwd("~/seqs/done/ITS")
temp = list.files(pattern="*.fasta")
list2 <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta")


ITS_whole.bad <- do.call(rbind, list)
ITS_whole <- as.data.frame(sapply(ITS_whole.bad, function(x) gsub("-", "", x))) %>% mutate(species_fasta = as.character(species_fasta)) %>% select(species, species_fasta, accession)
ITS_whole$oldName <- ITS_whole.bad$species
ITS_whole <- ITS_whole[-543, ]
sapply(ITS_whole, class)

df_list <- list()
names <- names(list2)
for(i in 1:length(list2)){
  fasta.tmp <- list2[[i]]
  sp <- data.frame(species_fasta = labels(fasta.tmp))
  species <- as.data.frame(sapply(sp, function(x) gsub("-", "", x))) %>% mutate(species_fasta = as.character(species_fasta))
  sapply(species, class)
  done <- inner_join(ITS_whole, species, by = "species_fasta")
  done$gene <- names[i]
  df_list[[i]] <- done
}

big_ITS = do.call(rbind, df_list)



####FIN!
big_compartments$species_fasta <- NA
bigDone <- rbind(big_compartments, big_ITS)

speciesList <- bigDone %>% select(species) %>% distinct()
write.csv(speciesList, "species.csv", row.names=F)
write.csv(bigDone, "myaccessions.csv", row.names=F)



##ITSTEST#####
#df.tmp.bad <- list[[1]]
#df.tmp <- as.data.frame(sapply(df.tmp.bad, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species)) %>% select(species_fasta, accession)
#df.tmp$oldName <- df.tmp.bad$species
fasta.tmp <- list2[[1]]
sp <- data.frame(species_fasta = labels(fasta.tmp))
species <- as.data.frame(sapply(sp, function(x) gsub("-", "", x))) %>% mutate(species_fasta = as.character(species_fasta))
sapply(species, class)
done <- inner_join(ITS_whole, species, by = "species_fasta")
#not_in <- anti_join(species, seq_in, by = "species_fasta")
#not_accessions <- left_join(not_in, df.tmp, by = "species_fasta") 
#joined <- rbind(mito_in, not_accessions)
names <- names(list2)
done$gene <- names[1]
df_list <- list()
df_list[[1]] <- done






#one try: #######
df.tmp.bad <- list[[3]]  
df.tmp <- as.data.frame(sapply(df.tmp.bad, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species)) %>% select(species, accession)
df.tmp$oldName <- df.tmp.bad$species
fasta.tmp <- list2[[3]]
sp <- data.frame(species = labels(fasta.tmp))
sp2 <- as.data.frame(sapply(sp, function(x) gsub("-", "", x))) %>% mutate(species = as.character(species))
species <- as.data.frame(sapply(sp2, function(x) gsub("_", " ", x))) %>% mutate(species = as.character(species))
sapply(species, class)
mito_in <- inner_join(mito, species, by = "species")
not_in <- anti_join(species, mito_in, by = "species")
not_accessions <- left_join(not_in, df.tmp, by = "species") 
joined <- rbind(mito_in, not_accessions)
names <- names(list)
joined$gene <- names[3]
df_list <- list()
df_list[[1]] <- joined



#rename seqs##############
rm(list = ls())
library(dplyr)
library(openxlsx)
library(phylotools)
library(ape)
rename_seqs <- read.xlsx("rename_seqs.xlsx") 
fixSeq <- rename_seqs %>% select(species_OG, speciesFinal_rename)
###OKAY FIRST FIX ACCESSIONS
accessions <- read.csv("myaccessions.csv")
accessions$species_OG <- gsub(" ", "_",accessions$species)
accessions <- left_join(accessions, fixSeq)
accessions$speciesFinal_rename[is.na(accessions$speciesFinal_rename)] <- paste("DELETE")
write.xlsx(accessions, "accessionsFin.xlsx", asTable = F)

check <- fixSeq[-grep("_", fixSeq$species_OG),]
check <- fixSeq[-grep("_", fixSeq$speciesFinal_rename),] #should be 5

###BEGIN#####
rm(list = ls())
setwd("~/seqs/done/mito")
temp = list.files(pattern="*.fasta")
list2env(lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta"), envir = .GlobalEnv)
setwd("~seqs/done/plastid")
temp = list.files(pattern="*.fasta")
list2env(lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta"), envir = .GlobalEnv)
setwd("~seqs/done/ITS")
temp = list.files(pattern="*.fasta")
list2env(lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta"), envir = .GlobalEnv)

concat <- cbind(atp1, atp6, atpB, atpF, ITS1, ITS2, nad2, nad5.nad4, ndhF, psbA, psbB, psbE, rbcL, rpl16, rpl2, rpl32, rpoC1, rps4, trnE, trnG, trnK, trnLF, trnM, trnS, x26s, x58s, fill.with.gaps=TRUE)
write.dna(concat, file = "~seqs/done/concat.fasta", format = "fasta")

##########
rm(list = ls())
rename_seqs <- read.xlsx("rename_seqs.xlsx") 
fixSeq <- rename_seqs %>% select(species_OG, speciesFinal_rename)
seq.tmp <- read.dna("concat.fasta", format = "fasta")
labelsDF <- data.frame(species_OG = labels(seq.tmp))
fixLabelsDF <- left_join(labelsDF, fixSeq)
fixLabelsDF$speciesFinal_rename[is.na(fixLabelsDF$speciesFinal_rename)] <- paste("DELETE")
newLabels <- fixLabelsDF$speciesFinal_rename
row.names(seq.tmp) <- newLabels
write.dna(seq.tmp, "concat_rename_raw.fasta", format = "fasta")


#fine do it without deleting
temp = list.files(pattern="*.fasta") #should be 26
list <- lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
               read.dna, format = "fasta")
rm(temp)


setwd("~seqs/renameMe/done")
for(i in 1:length(list)){
  seq.tmp <- list[[i]]
  labelsDF <- data.frame(species_OG = labels(seq.tmp))
  fixLabelsDF <- left_join(labelsDF, fixSeq)
  fixLabelsDF$speciesFinal_rename[is.na(fixLabelsDF$speciesFinal_rename)] <- paste("DELETE")
  newLabels <- fixLabelsDF$speciesFinal_rename
  row.names(seq.tmp) <- newLabels
  write.dna(seq.tmp, paste(names(list[i]), ".fasta", sep = ""), format = "fasta")
}

i = 3
seq.tmp <- list[[i]]
labelsDF <- data.frame(species_OG = labels(seq.tmp))
fixLabelsDF <- left_join(labelsDF, fixSeq)
fixLabelsDF$speciesFinal_rename[is.na(fixLabelsDF$speciesFinal_rename)] <- paste("DELETE")
newLabels <- fixLabelsDF$speciesFinal_rename
char.tmp <- as.character(seq.tmp)
names(char.tmp) <- newLabels
labels(char.tmp)
done <- as.DNAbin(char.tmp)
labels(done)
names(list[i])
write.dna(seq.tmp, paste(names(list[i]), ".fasta", sep = ""), format = "fasta")


row.names(seq.tmp) <- newLabels




fin <- as.DNAbin(try)
write.dna(fin, "fin.fasta", format = "fasta")

#####concatenate####################
setwd("~/seqs/done")
temp = list.files(pattern="*.fasta")
list2env(lapply(setNames(temp, make.names(gsub("*.fasta$", "", temp))), 
                read.dna, format = "fasta"), envir = .GlobalEnv)

concat <- cbind(atp1, atp6, atpB, atpF, ITS1, ITS2, nad2, nad5.nad4, ndhF, psbA, psbB, psbE, rbcL, rpl16, rpl2, rpl32, rpoC1, rps4, trnE, trnG, trnK, trnLF, trnM, trnS, x26s, x58s, fill.with.gaps=TRUE)
write.dna(concat, file = "~/seqs/done/concat.fasta", format = "fasta")

