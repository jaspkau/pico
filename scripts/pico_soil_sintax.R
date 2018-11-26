#setwd("/Users/administrator/Desktop/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)
library(ggdendro)
library(dendextend)

source("scripts/make_phyloseq.R")
source("scripts/soil_phyloseq.R")

d
###subset soil samples
decon.d
#decontamination
decon
####subsetetting for months and scaling of envt data
d_s

d_s = subset_samples(d_s, Pop_size != "Con")

d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

d_s.all = d_s
d_s.all = tax_glom(d_s.all, "Kingdom")
taxa_sums(d_s.all)

###select OMF families (Dearnaley, 2012 and Waud et al. 2014)

d_s = subset_taxa(d_s, Family=="f:Serendipitaceae"| Family=="f:Sebacinaceae"| Family=="f:Thelephoraceae"| Family=="f:Tulasnellaceae"| 
                    Family=="f:Ceratobasidiaceae"| Family=="f:Pezizaceae"| 
                    Family=="f:Pyronemataceae"| Family=="f:Tuberaceae"| Family=="f:Agaricaceae"| 
                    Family=="f:Clavulinaceae"| Family=="f:Corticiaceae"| Family=="f:Inocybaceae"| 
                    Family=="f:Marasmiaceae"| Family=="f:Russulaceae"| Family=="f:Tricholomataceae"| 
                    Family=="f:Typhulaceae"| Family=="f:Physalacriaceae"| Family=="f:Psathyrellaceae"|
                    Family=="f:Hymenogastraceae"| Family=="f:Incertae sedis")
d_s
d_s.omf = d_s

###proportion of OMF OTUs in soil

d_s.omf2 = tax_glom(d_s.omf, "Kingdom")
taxa_sums(d_s.omf2)
dprop = data.frame(sample_sums(d_s.omf2)/sample_sums(d_s.all))
summary(dprop)

####Thelephoraceae OTUs
d.thel = subset_taxa(d_s, Family == "f:Thelephoraceae")
d.thel
d.thel = tax_glom(d.thel, "Family")
taxa_sums(d.thel)

####Tul OTUs
d.tul = subset_taxa(d_s, Family == "f:Tulasnellaceae")
d.tul
d.tul = tax_glom(d.tul, "Family")
taxa_sums(d.tul)

####Cerato OTUs
d.cer = subset_taxa(d_s, Family == "f:Ceratobasidiaceae")
d.cer
d.cer = tax_glom(d.cer, "Family")
taxa_sums(d.cer)

####Agaricaceae OTUs
d.aga = subset_taxa(d_s, Family == "f:Agaricaceae")
d.aga
d.aga = tax_glom(d.aga, "Family")
taxa_sums(d.aga)

d_s = subset_taxa(d_s, Family == "f:Ceratobasidiaceae"| 
                    Family == "f:Tulasnellaceae")
d_s

d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

d.an = d_s
source("scripts/differential_abundance.R")
an.otus = anc$detected

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d_s, measures = c("Shannon", "Simpson"))
temp = merge(sample_data(d_s), aldiv, by = "row.names")
temp = temp[,-1]
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
for(i in c(5, 10)){
  column = names(temp[i])
  k.demo = kruskal.test(ef.sha ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")
alpha.kw

avg = temp %>%
  group_by(Population) %>%
  summarise(simp = mean(ef.sha))
avg

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

#ajust P-values
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

otu3_tab = otu_table(as.matrix(otu3), taxa_are_rows = F)
d.hc = merge_phyloseq(tax2, otu_table(as.matrix(otu3_tab), 
                                      taxa_are_rows = F), sample_data(d.int))

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h) %>% set("labels_cex", 0.5)
ggd1 <- as.ggdend(dhc)

p1 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p1

d.an =  d_s 

source("scripts/differential_abundance.R")

####prune the otu table according detected through ANCOM results
#to retain otus which need to be shown in heatmap

an.otus = anc$detected
otu.hm.an = otu3[,colnames(otu3) %in% an.otus]

otu.hm = merge(t(otu.hm.an), tax_table(d.hc), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
otu.hm$rank = paste(as.character(otu.hm$Row.names),"|",substr(otu.hm$Family, 3, 5))
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g3 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100),
               trace = "none")

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_s.omf, taxrank = "Family")
d_f = merge_samples(d_f, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = as.character(gen_f$Family)
#gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
gen_f$rank = ifelse(gen_f$Phylum == "unidentified", paste(as.character(gen_f$Kingdom), as.character(gen_f$Phylum), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Phylum != "unidentified" &  gen_f$Class == "unidentified", paste(as.character(gen_f$Phylum), as.character(gen_f$Class), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Class != "unidentified" &  gen_f$Order == "unidentified", paste(as.character(gen_f$Class), as.character(gen_f$Order), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Order != "unidentified" &  gen_f$Family == "unidentified", paste(as.character(gen_f$Order), as.character(gen_f$Family), sep = ";"), gen_f$rank)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
names(gen_f) = list
gen_f = gen_f/rowSums(gen_f)
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:16]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 5, name = "Blues"),brewer.pal(n = 10, name = "Paired"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "gold1", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 11, name = "Spectral")))

library(scales)

p.fam = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) +
  theme_bw(base_size = 20) + state_col2 + xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 10, color = "black")) +
  theme(legend.text = element_text(face = "italic", size = 10)) + guides(fill = guide_legend(ncol = 2, reverse=T, keywidth = 0.5, keyheight = 0.4))+ scale_y_continuous(labels = percent_format())
p.fam$data$variable = factor(p.fam$data$variable, ordered = TRUE, levels = rev(who))
p.fam
#ggsave(file="jc.treatment.nms.jpg")