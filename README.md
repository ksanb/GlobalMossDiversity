# GlobalMossDiversity
R scripts from Spatial phylogenetics of mosses (Bryopsida) at a Global Scale: Current Status and Future Directions

#readme

bryogeography.sh contains a script to parse out the accessions and sequences from the fasta file provided by Rose et al. (2016)

getSequences2.R to obtain all sequences used in Rose et al. (2016), obtain new ones, write fastas, and concatenate into a supermatrix.

Nomenclature.R to deal with the nomenclature of all species present in the phylogeny. Manual edits were made, thus it requires two CSV's: (species_ids_man.csv and species_nomenclatureMAN.csv)

treesFin.R to root and bootstrap the trees

getDistributionData.R to get distribution data from GBIF and iDigBio

combineAndCleanDistData.R to combine and clean the distribution data

compareGeffert.R to compare floristic and specimen-based datasets

calculatemetrics.R to calculate biodiversity metrics

phyloregion.R calculates phyloregions

