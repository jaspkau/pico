#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")
#setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
#setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

###ROOT OMF ANALYSIS......................................

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")

d
###subset root samples
decon.d
#decontamination
decon
####subsetetting for months and scaling of envt data
d_r
taxa_names(d_r) = gsub("otu", "denovo", taxa_names(d_r))

####Tul OTUs
d.tul = subset_taxa(d_r, Family == "f:Tulasnellaceae")
d.tul

library(seqinr)
fas = read.fasta(file = "data/chimera_filtered_rep_set.fasta", forceDNAtolower = FALSE, set.attributes = FALSE)
fas2 = fas[names(fas) %in% taxa_names(d.tul)]
write.fasta(sequences=fas2, names= names(fas2) ,file.out="Tul_OTUs.fasta")

####Cerato OTUs
d.cer = subset_taxa(d_r, Family == "f:Ceratobasidiaceae")
d.cer

fas2 = fas[names(fas) %in% taxa_names(d.cer)]
write.fasta(sequences=fas2, names= names(fas2) ,file.out="Cer_OTUs.fasta")
