#setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
d_r = subset_taxa(d_r, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d_r

d2 = merge_samples(d_r, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int = vegdist(otu3, method = "bray")
dist_w_int[is.na(dist_w_int)] <- 0

####make soil phyloseq

source("scripts/soil_phyloseq.R")

d_s = subset_samples(d_s, Population == "PLF"| Population == "PLE"|
                     Population == "SCW"| Population == "SCE"|
                       Population == "CH"| Population == "MX")

d_s = subset_taxa(d_s, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d_s

d2 = merge_samples(d_s, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int.s = vegdist(otu3, method = "bray")
dist_w_int.s[is.na(dist_w_int.s)] <- 0

# MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d_r)
names = c("P2", "MN", "SAND", "ZN")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d4)
names = c("pg.sm2", "dr.st", "pg.st")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

soil.otus = gen_f[1:6,] ###need to create soil phyloseq object first

library(ecodist)

MRM(dist_w_py ~ dist(my.env2) + dist(my.soil2) + dist(soil.otus), nperm=1000)
summary(lm(dist_w_py ~ dist(my.env2) + dist(my.soil2) + dist(soil.otus)))

mrm.env = lm(dist_w_py ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil = lm(dist_w_int ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

MRM(dist_w_int ~ dist_w_int.s, nperm=1000)
mrm.soil.otus = lm(dist_w_int ~ dist_w_int.s)
summary(mrm.soil.otus)$adj.r.squared
