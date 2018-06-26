#setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")
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

d_s.all = d_s

# Rarefaction -------------------------------------------------------------

d_rf = merge_samples(d_s, "Population")
otu_rf = data.frame(otu_table(d_rf))

library(iNEXT)
otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000, 70000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g1 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g1

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

d_s.all = tax_glom(d_s.all, "Kingdom")
d_s.omf = tax_glom(d_s.omf, "Kingdom")
dprop = data.frame(sample_sums(d_s.omf)/sample_sums(d_s.all))
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

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_s, taxrank = "Family")
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
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:19]
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