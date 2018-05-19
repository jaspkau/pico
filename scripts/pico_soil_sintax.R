setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")
source("scripts/soil_phyloseq.R")

d
###subset soil samples
decon.d
#decontamination
decon
####subsetetting for months and scaling of envt data
d_s

# Alpha diversity ---------------------------------------------------------

temp = estimate_richness(d_s)
temp = merge(met, temp, by = "row.names")

a = summary(aov(Shannon ~ Population + Year + Month + Pop_size + Demo, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad
a = summary(aov(Simpson ~ Population + Year + Month + Pop_size + Demo, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad

plot_richness(d_r, x= "Population", measures=c("Shannon", "Simpson"))

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

#permanova with significant factors. 
a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size + sample_data(d.bd)$Population + sample_data(d.bd)$Year + sample_data(d.bd)$Demo + sample_data(d.bd)$Stage, strata = "Pop_size", permutations = 999)
a
###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")
