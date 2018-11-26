setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")
setwd("/Users/administrator/Desktop/pico/")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
source("scripts/soil_phyloseq.R")

d.comb = merge_phyloseq(d_r, d_s)
d.comb

sample_data(d.comb)$int = paste(sample_data(d.comb)$Population,".",sample_data(d.comb)$Year)
sample_data(d.comb)$int = gsub(" ","",sample_data(d.comb)$int)

d.fin = subset_taxa(d.comb, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d.fin

# NETWORK ANALYSIS -------------------------------------------------------
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)

##make root otu file
##"Remove taxa not seen more than 3 times in at least 5% of the samples. 
#This protects against an OTU with small mean & trivially large C.V.
#d.r.net = filter_taxa(d.fin, function(x) sum(x > 1) > (0.05*length(x)), TRUE)

d.r.net = subset_samples(d.fin, Source == "R")
d.r.net
d.r.net = prune_taxa(taxa_sums(d.r.net) >= 1, d.r.net)
d.r.net
d.r.net = merge_samples(d.r.net, "int")
otu.r.net = data.frame(otu_table(d.r.net)) ##it should be non-normalized
###slelect OTUs first from which you want to build network, rows shud be samples
#otu.r.net = otu.r.net[row.names(otu.r.net) %in% row.names(sim.kw.popsize),]
colnames(otu.r.net) = paste(gsub("otu", "r", colnames(otu.r.net)))
taxa_names(d.r.net) = paste(gsub("otu", "r", taxa_names(d.r.net)))

d.s = subset_samples(d.fin, Population == "PLF"|Population == "PLE"
                     |Population == "SCW"| Population == "SCE"
                     |Population == "MX"| Population == "CH")

##make soil otu file
d.s.net = subset_samples(d.s, Source == "S")
d.s.net = prune_taxa(taxa_sums(d.s.net) >= 1, d.s.net)
d.s.net
#d.s.net = filter_taxa(d.s.net, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
d.s.net = merge_samples(d.s.net, "int")
d.s.net
otu.s.net = data.frame(otu_table(d.s.net)) ##it should no non-normalized
colnames(otu.s.net) = paste(gsub("otu", "s", colnames(otu.s.net)))
taxa_names(d.s.net) = paste(gsub("otu", "s", taxa_names(d.s.net)))

d.net = merge_phyloseq(d.r.net, d.s.net)

###Sparcc
##merge root and soil otu files together

net = cbind(ottau.s.net, otu.r.net, by = "row.names")
net = net[,-ncol(net)]

net = t(net)
write.csv(net, file = "results/net.csv")
#convert csv to .txt file
#and then do analyses with spaccWrapper.sh script and then import the files with following codes for editing to use in Cytoscape

library(reshape2)
spec.cor = read.delim("results/net/sim_cor.txt", sep = "\t")
row.names(spec.cor) = spec.cor[,1]
spec.cor = melt(spec.cor)

pval = read.delim("results/net/pvals_two_sided.txt", sep = "\t")
row.names(pval) = pval[,1]
pval = melt(pval)
library(splitstackshape)
library(plyr)
pval = plyr::rename(pval, c("value" = "p"))
spec.cor$p = pval$p

fam.df = as.data.frame(tax_table(d.net))
fam.df$otu = row.names(fam.df)

library(dplyr)
spec.cor = merge(spec.cor, fam.df, by.x = "X", by.y = "otu")
spec.cor$sfam = spec.cor$Family
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
spec.cor = spec.cor[ , !(names(spec.cor) %in% drops)]

spec.cor = merge(spec.cor, fam.df, by.x = "variable", by.y = "otu")
spec.cor$tfam = spec.cor$Family
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
spec.cor = spec.cor[ , !(names(spec.cor) %in% drops)]
spec.cor$fam = paste(substr(spec.cor$sfam, 3, 5),"-", substr(spec.cor$tfam, 3, 5))

spec.cor.t.sel = spec.cor[grepl("s", spec.cor$X), ] #target selection
spec.cor.fin.sel = spec.cor.t.sel[grepl("r", spec.cor.t.sel$variable), ]
write.csv(spec.cor.fin.sel, file = "results/net/spec.cor.csv")

###select for co-abundance network
spec.cor.fin.sel1 = spec.cor.fin.sel[spec.cor.fin.sel$value > 0.6,]
spec.cor.fin.sel1 = spec.cor.fin.sel1[spec.cor.fin.sel1$p == 0,]
write.csv(spec.cor.fin.sel1, file = "results/net/spec.cor.fin.sel.csv")

###select for co-exclusion network
spec.cor.fin.sel2 = spec.cor.fin.sel[spec.cor.fin.sel$value < -0.6,]
spec.cor.fin.sel2 = spec.cor.fin.sel2[spec.cor.fin.sel2$p == 0,]
write.csv(spec.cor.fin.sel2, file = "results/net/spec.cor.fin.sel_co_exclu.csv")

###select pair of root and soil OTUs
net = read.delim("results/net/spec.cor.csv", sep = ",")
net$X = gsub("s", "", net$X)
net$variable = gsub("r", "", net$variable)
net$sel = ifelse((net$X == net$variable)==TRUE, paste("sel"), paste(""))
#or
net$sel = ifelse(net$X == net$variable, paste("sel"), paste(""))
net = subset(net, sel == "sel")
write.csv(net, file = "results/net/net_sel.csv")

####SpiecEasi
spiec.out=spiec.easi(d.net, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(d.net)))
write.graph(spiec.graph,file="results/spieceasi.ncol.txt",format="ncol") 
plot_network(spiec.graph, d.net, type='taxa', color = "Family", label="value")
