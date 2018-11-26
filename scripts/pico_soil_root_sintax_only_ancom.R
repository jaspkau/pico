#setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(vegan)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(dplyr)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
d.an = d_r
source("scripts/differential_abundance.R")
an.otus = anc$detected

source("scripts/make_phyloseq.R")
source("scripts/soil_phyloseq.R")
d_s
d_s.all = d_s

d_s = subset_samples(d_s, Pop_size != "Con")

d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####SOIL OMF ANALYSIS WITH OTUs identified important by ANCOM in ROOTS ----------------------------------------------------------------

d_s = prune_taxa(taxa_names(d_s) %in% an.otus, d_s)
d_s

d_s = subset_samples(d_s, Pop_size == "L"| Pop_size == "S")
d_s

d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

d_s.pcoop = d_s

###proportion of Pcoop specific OTUs in soil

d_s.all = tax_glom(d_s.all, "Kingdom")
d_s.pcoop = tax_glom(d_s.pcoop, "Kingdom")
dprop = data.frame(sample_sums(d_s.pcoop)/sample_sums(d_s.all))
summary(dprop)

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d_s, measures = c("Shannon", "Simpson"))
temp = merge(sample_data(d_s), aldiv, by = "row.names")
temp = temp[,-1]
row.names(temp) = temp[,1]

bp <- ggplot(temp, aes(x=Year, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$Simpson)

alpha.kw = c()
for(i in c(5, 10)){
  column = names(temp[i])
  k.demo = kruskal.test(Simpson ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")

# Soil OMF Beta diversity with bray ------------------------------------------------

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
#site_list = colnames(otu)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = stepacross(dist_w)

###PERMANOVA

###Weighted distance

a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size + sample_data(d.bd)$Population, permutations = 999)
a
###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

# Hierarchial clustering --------------------------------------------------

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d.int = merge_samples(d_s, "pop.year")
otu3 = data.frame(otu_table(d.int))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = data.frame(otu3)
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int = vegdist(otu3, method = "bray")
dist_w_int[is.na(dist_w_int)] <- 0

otu3_tab = otu_table(as.matrix(otu3), taxa_are_rows = F)
d.hc = merge_phyloseq(tax2, otu_table(as.matrix(otu3_tab), 
                                      taxa_are_rows = F), sample_data(d.int))

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h) %>% set("labels_cex", 0.5)
ggd1 <- as.ggdend(dhc)

p1 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p1

otu.hm = merge(t(otu3), tax_table(d.hc), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
otu.hm$rank = paste(as.character(otu.hm$Row.names),"-",substr(otu.hm$Family, 3, 5))
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g2 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100),
               trace = "none")
