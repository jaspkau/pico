setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
d_r

source("scripts/soil_phyloseq.R")
d_s

# ###SOIL OMF ANALYSIS WITH OTUs identified in ROOTS ----------------------------------------------------------------

dr = data.frame(otu_table(d_r))
dr = t(dr)
who = colnames(dr)

###create soil phyloseq object by using above info

otu_s = data.frame(otu_table(d_s))
otu_s2 = otu_s[who,]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_s = merge_phyloseq(tax2, otu_s4, sample_data(d_s))
d_s
d_s = subset_samples(d_s, Population == "PLF"| Population == "PLE"|
                       Population == "SCW"| Population == "SCE"| 
                       Population == "MX"| Population == "CH")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

# ###SOIL OMF ANALYSIS WITH the MOST abundant root OTUS----------------------------------------------------------------

d_f = merge_samples(d_r, "Population")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(gen_f$Row.names)
list = as.character(gen_f$rank)
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]

###create soil phyloseq object by using above info

otu_s = data.frame(otu_table(d))
otu_s2 = otu_s[who,]
otu_s2 = otu_s2[,colSums(otu_s2) > 0]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_s = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_s
d_s = subset_samples(d_s, Source == "S")
d_s
d_s = subset_samples(d_s, Month == "Feb"| Month == "Apr")
d_s
d_s = subset_samples(d_s, Population == "PLF"| Population == "PLE"|
                       Population == "SCW"| Population == "SCE"| 
                       Population == "MX"| Population == "CH")
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####scale envt data according to above sample selection
met4 = as.data.frame(sample_data(d_s))
s_env_met = met4[,cbind(1,2,3,4,5,6,7,8,9,36,37,38)]
s_env = met4[,10:35]
s_env = scale(s_env)
met5 = merge(s_env_met, s_env, by = "row.names")
row.names(met5) = met5$Row.names
met5$int =  paste(met5$Population,".",met5$Pop_size,".",met5$Demo,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
#site_list = colnames(otu)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d2))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = stepacross(dist_w)

###PERMANOVA

###Weighted distance

a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size + sample_data(d.bd)$Population + sample_data(d.bd)$Year + sample_data(d.bd)$Demo, strata = "Pop_size", permutations = 999)
a
###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")


###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d.int = merge_samples(d_s, "int")
otu3 = data.frame(otu_table(d.int))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = data.frame(otu3)
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")
dist_w_int = dist.zeroes(otu3, dist_w_int)

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(10, 3), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               density.info = c("none"))

