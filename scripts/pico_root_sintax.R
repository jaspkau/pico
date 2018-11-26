#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")
#setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
#setwd("/Users/administrator/Desktop/pico")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)
library(ggdendro)
library(dendextend)

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

##subsetting for OMF OTUs
d.romf

droot = tax_glom(d_r, "Kingdom")
taxa_sums(droot)

####Tul OTUs
d.tul = subset_taxa(d_r, Family == "f:Tulasnellaceae")
d.tul
d.tul = tax_glom(d.tul, "Family")
taxa_sums(d.tul)

####Cerato OTUs
d.cer = subset_taxa(d_r, Family == "f:Ceratobasidiaceae")
d.cer
d.cer = tax_glom(d.cer, "Family")
taxa_sums(d.cer)

d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)

####the sequences of 50 most abundant OTUs
d_f = data.frame(otu_table(d_r))
gen_f = t(d_f)
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:50]
d_gf = prune_taxa(taxa_names(d_r) %in% who, d_r)
taxa_sums(d_gf)

####OTUs represented by <10 sequences

d_gf = prune_taxa(taxa_sums(d_r) < 10, d_r)
d_gf

# Rarefaction -------------------------------------------------------------

d.rf = merge_samples(d_r, "Population")
otu_rf = data.frame(otu_table(d.rf))

library(iNEXT)
otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g1 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g1

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met, aldiv, by = "row.names")
row.names(temp) = temp[,1]

# Once again, effective numbers to the resuce!
# The conversion of Simpson diversity to effective numbers is 1/1-D
temp$ef = 1/(1-temp$Simpson)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha = exp(temp$Shannon)

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$ef.sha)

alpha.kw = c()
for(i in c(6, 7, 11)){
  column = names(temp[i])
  k.demo = kruskal.test(ef.sha ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")
alpha.kw

avg = temp %>%
  group_by(Pop_size) %>%
  summarise(simp = mean(ef.sha))
avg


# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")

sample_data(d.bd)$int = paste(sample_data(d.bd)$Population, ".", sample_data(d.bd)$Year)
mod = betadisper(dist_w, sample_data(d.bd)$Year)
mod = betadisper(dist_w, sample_data(d.bd)$Population)
anova(mod)

###PERMANOVA

###Weighted distance

#permanova with significant factors. 

a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size + sample_data(d.bd)$Population + sample_data(d.bd)$Population + sample_data(d.bd)$Stage, permutations = 999)
a

###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

# Hierarchial clustering --------------------------------------------------

#compressing the phyloseq object at level which is significantly different

d.int = merge_samples(d_r, "int")
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

d.an = d_r

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

#Relative abundance plots ---------------------------------------------

d.ra = d_r

source("scripts/relative_abundance_plots.R")

###relative abundances at family level
p1 = p.fam

###relative abundances at OTU level
p2 = p.otus

ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2)


library(dplyr)
gavg <- temp %>% 
  group_by(Year) %>% 
  summarise(Dia = mean(Shannon))
gavg
