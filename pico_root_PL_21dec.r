setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run")
#source("C:\\Users\\jaspkaur\\Google Drive\\r\\packages.r")
source("C:\\Users\\jaspkaur\\Google Drive\\r\\libraries.r")

setwd("C://Users//jaspr//Google Drive//Metagenomics//pprec_july2017")
source("C:\\Users\\jaspr\\Google Drive\\r\\libraries.r")

#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
source("/Users/administrator/Documents/jaspreet/r/libraries.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/uclust_97%")

library(dunn.test)
library(adespatial)
library(phyloseq)
library(metagenomeSeq)
library(mixOmics)
library(sjPlot)

# Make phyloseq object ----------------------------------------------------

otu <- read.delim(file = "otu_table_no_singletons_sintax.txt", 
                  sep = "\t", header = T)
otu = otu[,-ncol(otu)]
row.names(otu) = otu[,1]
otu = otu[,-1]
#Rarefy(otu, depth = min(rowSums(otu)))
otu = otu[,colSums(otu) > 0]
site_list = colnames(otu)
otu_tab = otu_table(as.matrix(otu), taxa_are_rows = T)

#Format the taxonomy table

tax = read.delim(file = "parallel_blast/chimera_filtered_rep_set_tax_assignments2.txt", 
                 sep = "\t", header = F)
tax = tax[,-c(3:4)]
tax[,3:9] = colsplit(tax$V2, ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
tax = tax[,-2]
row.names(tax) = tax[,1]
tax = tax[,-1]
tax2 = tax_table(as.matrix(tax))
#source("R_scripts//read_newick.R")
#tre = read.newick("tree.tre")
#tre$tip.label = gsub("_", ":", tre$tip.label)
#tre$tip.label %in% row.names(otu)

# SINTAX taxonomy table ---------------------------------------------------

library(reshape2)
#####
#Format the taxonomy table
#####

tax = read.delim(file = "tax.sintax", sep = "\t", header = F)
#split first column to separate otu names
#tax[,1:2] = colsplit(tax$V1, " ", c("otu", "seq"))
#row.names(tax) = tax[,1]
row.names(tax) = tax$V1
list = tax$V2
tax2 = colsplit(list, pattern ="\\(|\\),", c("Kingdom", "Kingdom_conf", "Phylum", "Phylum_conf", "Class", "Class_conf", "Order", "Order_conf", "Family", "Family_conf", "Genus", "Genus_conf", "Species", "Species_conf"))
tax2$Species_conf = gsub("\\)", "", tax2$Species_conf)
tax2$Species_conf = as.numeric(tax2$Species_conf)
tax2[is.na(tax2)] <- 0
row.names(tax2) = row.names(tax)

###WE NEED to only use taxonomic assignments that are confident
##here is the pseudo code for what we need to do
# is taxonomy assignment confidence > 95% ?, if so, keep taxonomy otherwise assign unknown

tax_assign = function(x, conf, rank){
  # x = data frame containig taxonomies and confidences
  # conf = column with confidences
  # rank = corresponding taxonomic ranks
  
  #make a list of confidences as some level
  list = x[[conf]]
  ##make a list of ranks at some level (correspoding to confidences)
  list2 = as.character(x[[rank]])
  ###see if the confidence is <95%, if so, change the rank to unknown
  list3 = ifelse(list < 0.90, "unidentified", list2) ###if list is <95 than assign unknown otherwise go with list2
  #replace the original rank with the new filtered rank based on threshold
  x[rank] = list3
  return(x)
}

level = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
for (i in 1:7){ ###i 1 = kingdom, i 2 = phylum etc
  
  tax2 = tax_assign(tax2, paste(level[i], "_conf", sep = ""), level[i])
  
}

tax2 = tax2[,-c(2,4,6,8,10,12,14)]

tax2$row = row.names(tax2)
tax2[,8:9] = colsplit(tax2$row, " ", c("otu", "seq"))
row.names(tax2) = tax2$row
tax2 = tax2[,-c(8:9)]
tax2 = tax_table(as.matrix(tax2))

#meta data
library(gdata)
met <- read.xls("met.xlsx", sheet = 1)
row.names(met) = met[,1]
met$int = paste(met$Population,".",met$Stage,".",gsub("20", "", met$Year))
met$int = gsub(" ", "", met$int)
met$pop.year = paste(met$Population, ".", gsub("20", "", met$Year))
met$pop.year = gsub(" ", "", met$pop.year)

#phyloseq object

d = merge_phyloseq(tax2, otu_tab, sample_data(met))
d
d = subset_taxa(d, Kingdom == "d:Fungi")
d
d_r = subset_samples(d, Source == "R")
d_r
d_r = subset_samples(d_r, Month == "Feb"| Month == "April")
d_r
#d_r = subset_samples(d_r, Population == "PLF"| Population == "PLE" | Population == "PNPS")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d_r))
env_met = met2[,cbind(1,2,3,4,5,6,7, 62,63)]
env = met2[,8:60]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d_r = merge_phyloseq(tax_table(d_r), otu_table(d_r), sample_data(met3))
d_r

# march samples -----------------------------------------------------------

d_r = subset_samples(d, sample_data(d)$Month == "Mar")
# Beta diversity with bray ------------------------------------------------
#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
library(vegan)
otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_uw = vegdist(rel_otu_code, method = "bray", binary = TRUE)

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d3)$Population, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Stage, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$pop.year, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$int), permutations = 999)
a

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different

d2 = merge_samples(d_r, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")

otu3 = otu_table(as.matrix(otu3), taxa_are_rows = F)
d4 = merge_phyloseq(tax2, otu_table(as.matrix(otu3), 
                                    taxa_are_rows = F), sample_data(d2))

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
          cellnote = otu3, notecex=1.0,
          notecol="white")

# PCoA with bray --------------------------------------------------------------------

state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

####Weighted

pc = capscale(dist_w ~ 1, comm = rel_otu_code) ###~ means function of nothing
pc$CA$eig
s = summary(pc)
cap = data.frame(pc$CA$u)
plot(cap[,1], cap[,2])
cap = merge(sample_data(d), cap, by = "row.names")
#cap$Population = cap$Row.names
label = cap$int

p = ggplot(cap, aes(x= MDS1, y= MDS2, label=label))+theme_bw(base_size = 15) +
  geom_point() + state_col_ord + geom_text(aes(label=label),size = 3, vjust = "inward", hjust = "inward", check_overlap = TRUE) +
  labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
       y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))
p
# RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_code, met3, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]

test=forward.sel(fwdsel_df[,2:540], #OTUS#
                 fwdsel_df[,551:575], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_code ~ MG_PCT + MN + SAND + NA_PCT + ZN, data=fwdsel_df) ###this works for anova.cca
t = summary(cc)
#ZN + FE + MG + P2 + OM + MN + CEC + B + K + SC + NO3_N + P1 +SM + CU + NA. + CLAY + SILT + S + CA + PH, data=fwdsel_df) ###this works for anova.cca
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
ccdata = merge(ccdata, fwdsel_df, by = "row.names")
Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year)), size=2) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_code, met3, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]

test=forward.sel(fwdsel_df[,2:540], #OTUS#
                 fwdsel_df[,576:593], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_code ~ pg.sm2 + dr.at + dr.sm1 + smd.at + dr.sm2, data=fwdsel_df) ###this works for anova.cca
t = summary(cc)

cc = rda(rel_otu_code ~ MG_PCT + MN + SAND + NA_PCT + ZN + pg.sm2 + dr.at + dr.sm1, data=fwdsel_df) ###this works for anova.cca
t = summary(cc)

#ZN + FE + MG + P2 + OM + MN + CEC + B + K + SC + NO3_N + P1 +SM + CU + NA. + CLAY + SILT + S + CA + PH, data=fwdsel_df) ###this works for anova.cca
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)
arrowdata
#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
ccdata = merge(ccdata, fwdsel_df, by = "row.names")

state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year)), size=2) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()
# MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d2)
names = c("MG_PCT", "MN", "SAND", "NA_PCT", "ZN")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d2)
names = c("pg.sm2", "dr.at", "dr.sm1")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

library(ecodist)

MRM(dist_w_int ~ dist(my.env2) + dist(my.soil2), nperm=1000)
summary(lm(dist_w_int ~ dist(my.env2) + dist(my.soil2)))

mrm.env = lm(dist_w_int ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil.otus = lm(dist_w_int ~ dist(my.soil2))
summary(mrm.soil.otus)$adj.r.squared

# Realtive abundance plots at OTU level ------------------------------------------------

d_f = merge_samples(d_r, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
#gen_f$rank = as.character(gen_f$Family)
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
list = as.character(gen_f$rank)
#list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
#m$State2 = as.factor(m$State2)
#m$Time = ordered(met$Time, levels = c("pre-part", "post-part_3", "post-part_9", "post-part_25"))
#m$variable = gsub("_.*", "", m$variable)
#who = gsub("_.*", "", who)
state_col2 = scale_fill_manual(name = "State3", values=rev(c(brewer.pal(8, "Accent"), "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                             "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                             "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                             "slategray4", "seagreen4" , "aquamarine",
                                                             "tomato2")))
#col = scale_fill_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(19 - 11), brewer_pal(palette = "Spectral")(11))),"#CCCCCC"), name = "Family")
library(scales)
p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(reverse=TRUE))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = who)
p
# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_r, taxrank = "Family")
d_f = merge_samples(d_f, "pop.year")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = as.character(gen_f$Family)
#gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
#m$State2 = as.factor(m$State2)
#m$Time = ordered(met$Time, levels = c("pre-part", "post-part_3", "post-part_9", "post-part_25"))
#m$variable = gsub("_.*", "", m$variable)
#who = gsub("_.*", "", who)
state_col2 = scale_fill_manual(name = "State3", values=rev(c(brewer.pal(8, "Accent"), "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                             "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                             "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                             "slategray4", "seagreen4" , "aquamarine",
                                                             "tomato2")))
#col = scale_fill_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(19 - 11), brewer_pal(palette = "Spectral")(11))),"#CCCCCC"), name = "Family")
library(scales)
p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(reverse=TRUE))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = who)
p

# Soil OMF ----------------------------------------------------------------

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
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:30]

###create soil phyloseq object by using above info

otu_s = data.frame(otu_table(d))
otu_s2 = otu_s[who,]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_s = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_s
d_s = subset_samples(d_s, Source == "S")
d_s
d_s = subset_samples(d_s, Month == "Feb"| Month == "April")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax_table(d_s), otu2, sample_data(d_s))
rel_otu_code = data.frame(otu_table(d3))

library(BiodiversityR)
dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = dist.zeroes(rel_otu_code, dist_w)
dist_uw = vegdist(rel_otu_code, method = "bray", binary = TRUE)

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d3)$Population, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Year), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$int), permutations = 999)
a
# Realtive abundance plots ------------------------------------------------
#d_f = tax_glom(d_r, taxrank = "Family")
d_f = merge_samples(d_s, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
#gen_f$rank = paste(gen_f$Family)
list = as.character(gen_f$rank)
#list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:26]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
#m$State2 = as.factor(m$State2)
#m$Time = ordered(met$Time, levels = c("pre-part", "post-part_3", "post-part_9", "post-part_25"))
#m$variable = gsub("_.*", "", m$variable)
#who = gsub("_.*", "", who)

state_col2 = scale_fill_manual(name = "State3", values=rev(c(brewer.pal(8, "Accent"), "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                             "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                             "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                             "slategray4", "seagreen4" , "aquamarine",
                                                             "tomato2")))
#col = scale_fill_manual(values = c(rev(c(brewer_pal(palette = "Dark2")(19 - 11), brewer_pal(palette = "Spectral")(11))),"#CCCCCC"), name = "Family")
library(scales)
p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(reverse=TRUE))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = who)
p


