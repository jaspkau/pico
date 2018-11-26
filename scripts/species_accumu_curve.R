setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")
setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

#########read phyloseq object
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))

source("scripts/make_phyloseq_object.R")

decon = subset_samples(d, Source == "root")
decon
decon = subset_samples(decon, Population2 != "X")
decon

####decontaminate phyloseq object based on frequency
source("scripts/decontaminate_phyloseq.R")

d_r = subset_samples(decon.d, Sample_or_Control == "Sample")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

###make pc soil phyloseq

decon = subset_samples(d, Source == "soil")
decon
decon = subset_samples(decon, Species == "P. cooperi")
decon
decon = prune_taxa(taxa_sums(decon) >= 1, decon)
decon

####decontaminate phyloseq object based on frequency
source("scripts/decontaminate_phyloseq.R")

d.pc = subset_samples(decon.d, Sample_or_Control == "Sample")
d.pc
d.pc = prune_taxa(taxa_sums(d.pc) >= 1, d.pc)
d.pc

###make pp soil phyloseq

d.pp = subset_samples(d, Source == "soil")
d.pp
d.pp = subset_samples(d.pp, Species == "P. praeclara")
d.pp
d.pp = prune_taxa(taxa_sums(d.pp) >= 1, d.pp)
d.pp

####Combine phyloseq objects

d.fin = merge_phyloseq(d_r, d.pc, d.pp)
d.fin

taxa_names(d.fin) = gsub("otu", "denovo", taxa_names(d.fin))

sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source,".",sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Source)
sample_data(d.fin)$int = paste(sample_data(d.fin)$Species,".",sample_data(d.fin)$Stage)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

##merge samples
dpop = merge_samples(d.fin, "int")
  
# Species accumulation curve -------------------------------------------------------

otu_rc = data.frame(otu_table(dpop)) ###rows should be samples
sac = specaccum(otu_rc, method = "collector")
rarecurve(otu_rc, col = "blue", cex = 0.6, label = FALSE, xlab = "No. of sequences", ylab = "Observed no. of OTUs")
raremax <- min(rowSums(otu_rc))
rarecurve(otu_rc, step = 50, sample = raremax, col = "blue", cex = 0.6, xlab = "No. of sequences", ylab = "Observed no. of OTUs")

###
library(iNEXT)
otu_rc = data.frame(t(otu_rc)) ####columns should be samples
m <- c(1, 5, 20, 50, 100, 200, 400, 1000, 1500, 2000, 3000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m)
ggiNEXT(out, type=1, facet.var="none", color.var="site")
