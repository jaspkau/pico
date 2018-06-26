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

###ROOT OMF ANALYSIS......................................

source("scripts/make_phyloseq.R")
d
d.r = subset_samples(d, Source == "R")
d.r = subset_samples(d.r, Month == "Feb" | Month == "Apr")
d.r = prune_taxa(taxa_sums(d.r) >= 1, d.r)
d.r
d.r.py = merge_samples(d.r, "pop.year")
  
# Species accumulation curve -------------------------------------------------------

otu_rc = data.frame(otu_table(d.r.py)) ###rows should be samples
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
