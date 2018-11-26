#setwd("/Users/administrator/Desktop/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")
#setwd("/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

######################################

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

droot = tax_glom(d_r, "Kingdom")
taxa_sums(droot)

##subsetting for OMF OTUs

####Tul OTUs
d.tul = subset_taxa(d_r, Family == "f:Tulasnellaceae"| Family == "f:Ceratobasidiaceae")
d.tul
d.tul = tax_glom(d.tul, "Family")
taxa_sums(d.tul)


aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

# Once again, effective numbers to the resuce!
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef = 1/(1-temp$Simpson)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha = exp(temp$Shannon)


aldiv = estimate_richness(d.tul, measures = c("Shannon", "Simpson"))
temp = merge(temp, aldiv, by = "row.names")
row.names(temp) = temp[,1]

# Once again, effective numbers to the resuce!
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef.y = 1/(1-temp$Simpson.y)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha.y = exp(temp$Shannon.y)

cor.test(temp$ef.sha, temp$ef.sha.y)

d.r2 = prune_taxa(!(taxa_names(d_r) %in% taxa_names(d.tul)), d_r)
d.r2

aldiv = estimate_richness(d.r2, measures = c("Shannon", "Simpson"))
temp = merge(temp, aldiv, by = "row.names")
row.names(temp) = temp[,1]

# Once again, effective numbers to the resuce!
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef.z = 1/(1-temp$Simpson)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha.z = exp(temp$Shannon)

cor.test(temp$ef.sha.y, temp$ef.sha.z)

######################################################################

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
source("scripts/soil_phyloseq.R")

d.comb = merge_phyloseq(d_r, d_s)
d.comb

d.comb = subset_samples(d.comb, Pop_size != "Con")

d.comb = prune_taxa(taxa_sums(d.comb) > 0, d.comb)

seq_dep = estimate_richness(d.comb, measures = "Observed")
seq_dep$sequences = sample_sums(d.comb)
seq_dep = merge(sample_data(d.comb), seq_dep, by = "row.names")

sample_data(d.comb)$int = paste(sample_data(d.comb)$Source,".",sample_data(d.comb)$Population,".",sample_data(d.comb)$Year)

d.fin = subset_taxa(d.comb, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d.fin

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d.fin))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,10,11,37,38,39)]
env = met2[,12:36]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d.fin = merge_phyloseq(tax_table(d.fin), otu_table(d.fin), sample_data(met3))

###identify the overlap between root and soil Otus

d.r = subset_samples(d.fin, Source == "R")
d.r
d.r = prune_taxa(taxa_sums(d.r) >= 1, d.r)
d.r
d.r2 = tax_glom(d.r, "Kingdom")
taxa_sums(d.r2)

d.s.x = subset_samples(d.fin, Source == "S")
d.s.x
d.s.x = prune_taxa(taxa_sums(d.s.x) >= 1, d.s.x)
d.s.x
d.r.s = prune_taxa((taxa_names(d.r) %in% taxa_names(d.s.x)), d.r)
d.r.s

d.r.s2 = tax_glom(d.r.s, "Kingdom")
taxa_sums(d.r.s2)

library(vegan)
otu2 = t(data.frame(otu_table(d.r.s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d.r.s))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")

sample_data(d.bd)$int = paste(sample_data(d.bd)$Population, ".", sample_data(d.bd)$Year)
mod = betadisper(dist_w, sample_data(d.bd)$Year)
mod = betadisper(dist_w, sample_data(d.bd)$Population)
mod = betadisper(dist_w, sample_data(d.bd)$int)

anova(mod)

###PERMANOVA

###Weighted distance

#permanova with significant factors. 

a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")


d.int = merge_samples(d.r.s, "Population")
otu3 = data.frame(otu_table(d.int))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")

otu3_tab = otu_table(as.matrix(otu3), taxa_are_rows = F)
d.hc = merge_phyloseq(tax2, otu_table(as.matrix(otu3_tab), 
                                      taxa_are_rows = F), sample_data(d.int))

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h) %>% set("labels_cex", 0.5)
ggd1 <- as.ggdend(dhc)

d.an = d.r.s

source("scripts/differential_abundance.R")

####prune the otu table according to simper analyses results
#to retain otus which need to be shown in heatmap

###detected through ANCOM
an.otus = anc$detected

otu.hm.an = otu3[,colnames(otu3) %in% an.otus]
write.csv2(otu.hm.an, file = "data/root_ancom_otus.csv")

otu.hm = merge(t(otu.hm.an), tax_table(d.hc), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
otu.hm$rank = paste(as.character(otu.hm$Row.names),"|",substr(otu.hm$Family, 3, 5))
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g1 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu.hm, notecex=0.4,
               notecol="white")
