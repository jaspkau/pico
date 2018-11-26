setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")

d.abun = merge_samples(d, "Source")
mytaxa = taxa_names(d_r)
d.abun2 = prune_taxa(mytaxa, d.abun)

otu3 = data.frame(otu_table(d.abun2))
otu3 = decostand(otu3, method = "hellinger")
otu3 = otu_table(otu3, taxa_are_rows = FALSE)
d.abun3 = merge_phyloseq(otu3, tax_table(d.abun2), sample_data(d.abun2))
test = data.frame(otu_table(d.abun3))
test2 = data.frame(t(test))
plot(test2)
label = row.names(test2)

p = ggplot(test2, aes(R, S)) + geom_point() + 
  geom_text(aes(label=label),size = 3, vjust = "inward", 
            hjust = "inward", check_overlap = TRUE)
p

#######for only cer and tul

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
source("scripts/soil_phyloseq.R")

d.comb = merge_phyloseq(d_r, d_s)
d.comb

sample_data(d.comb)$int = paste(sample_data(d.comb)$Source,".",sample_data(d.comb)$Population,".",sample_data(d.comb)$Year)

d.fin = subset_taxa(d.comb, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d.fin

otu3 = data.frame(otu_table(d.fin))
otu3 = decostand(otu3, method = "hellinger")
otu3 = otu_table(otu3, taxa_are_rows = TRUE)
d.abun3 = merge_phyloseq(otu3, tax_table(d.fin), sample_data(d.fin))
test = data.frame(otu_table(d.abun3))
test2 = data.frame(t(test))
plot(test2)
label = row.names(test2)

p = ggplot(test2, aes(R, S)) + geom_point() + 
  geom_text(aes(label=label),size = 3, vjust = "inward", 
            hjust = "inward", check_overlap = TRUE)
p

d.abun = merge_samples(d.fin, "Source")

otu3 = data.frame(otu_table(d.abun))
otu3 = decostand(otu3, method = "hellinger")
otu3 = otu_table(otu3, taxa_are_rows = FALSE)
d.abun3 = merge_phyloseq(otu3, tax_table(d.abun), sample_data(d.abun))
test = data.frame(otu_table(d.abun3))
test2 = data.frame(t(test))
plot(test2)
label = row.names(test2)

p = ggplot(test2, aes(R, S)) + geom_point() + 
  geom_text(aes(label=label),size = 3, vjust = "inward", 
            hjust = "inward", check_overlap = TRUE)
p


                  