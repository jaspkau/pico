#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(vegan)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(dplyr)

#source("scripts/make_phyloseq.R")
#source("scripts/root_phyloseq.R")

# ANCOM (Analysis of composition of microbiome)-------------------------------------------------------------------

ancom.otu = t(data.frame(otu_table(d.an))) ##columns = OTUs and should be counts
ancom.otu = merge(ancom.otu, sample_data(d.an), by = "row.names")
row.names(ancom.otu) = ancom.otu$Code
ancom.otu = ancom.otu[,-1]
names(ancom.otu)
##look for the grouping variable you want to use
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Pop_size", names(ancom.otu))]

anc = ANCOM(ancom.fin, multcorr = 1, sig = 0.05)
anc$detected
plot_ancom(anc)

#SIMPER (similaritypercentage) analyses ----------------------------------------------

otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)
d.bd = merge_phyloseq(tax2, otu2, sample_data(d))

ind.df = data.frame(otu2)##the taxa should be columns and this otu table is hellinger tranfromed

#Identification of species most responsible for differences among groups of samples
#SIMPER(similaritypercentage), Based on abundance, does not weigh occurrence frequency as indicator species analysis does.

sim = simper(ind.df, sample_data(d.bd)$Pop_size)
sim.sum = summary(sim)
sim.df.popsize = data.frame(sim.sum$S_L)

sim.popsize.otus = row.names(sim.df.popsize)[1:50]

mann.popsize.df = ind.df[,names(ind.df) %in% sim.popsize.otus]

mann.popsize.df2 = merge(mann.popsize.df, met3, by = "row.names")
mann.popsize.df2$Pop_size = as.factor(mann.popsize.df2$Pop_size)

sim.kw.popsize = c()
for(i in 2:51){
  column = names(mann.popsize.df2[i])
  k = kruskal.test(mann.popsize.df2[,i]~Pop_size, data = mann.popsize.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k)))
  sim.kw.popsize = rbind(sim.kw.popsize, results)
} 

sim.kw.popsize$p.ad = p.adjust(sim.kw.popsize$pval, method = "bonferroni")

############OTUs differenting between demopgraphic classes

sim = simper(ind.df, sample_data(d.bd)$Demo)
sim.sum = summary(sim)
sim.df.demo = data.frame(sim.sum$F_NF)

sim.demo.otus = row.names(sim.df.demo)[1:50]

library(dplyr)

mann.demo.df = ind.df[,names(ind.df) %in% sim.demo.otus]

mann.demo.df2 = merge(mann.demo.df, met3, by = "row.names")
mann.demo.df2$Demo = as.factor(mann.demo.df2$Demo)

##Do kruskal wallis test with the first OTUs from simper analyses

sim.kw.demo = c()
for(i in 2:51){
  column = names(mann.demo.df2[i])
  k.demo = kruskal.test(mann.demo.df2[,i]~Demo, data = mann.demo.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  sim.kw.demo = rbind(sim.kw.demo, results)
} 

sim.kw.demo$p.ad = p.adjust(sim.kw.demo$pval, method = "bonferroni")


