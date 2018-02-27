setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

setwd("C://Users//jaspr//Google Drive//Metagenomics//pprec_july2017")

#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/uclust_97%")

library(dunn.test)
library(adespatial)
library(phyloseq)
library(metagenomeSeq)
library(mixOmics)

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

#Format blast taxonomy table

tax = read.delim(file = "parallel_blast/chimera_filtered_rep_set_tax_assignments2.txt", 
                 sep = "\t", header = F)
tax = tax[,-c(3:4)]
tax[,3:9] = colsplit(tax$V2, ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
tax = tax[,-2]
row.names(tax) = tax[,1]
tax = tax[,-1]
tax2 = tax_table(as.matrix(tax))

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
d_r = subset_samples(d_r, Month == "Feb"| Month == "Apr")
d_r
#d_r = subset_samples(d_r, Population == "PLF"| Population == "PLE" | Population == "PNPS")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d_r))
env_met = met2[,cbind(1,2,3,4,5,6,7,52,53)]
env = met2[,8:51]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d_r = merge_phyloseq(tax_table(d_r), otu_table(d_r), sample_data(met3))
d_r

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

####At pop.year level

d4 = merge_samples(d_r, "pop.year")
otu4 = data.frame(otu_table(d4))
otu4 = decostand(otu4, method = "hellinger")
rel_otu_pop.year = data.frame(otu_table(otu4, taxa_are_rows = F))

rowSums(rel_otu_pop.year)
#otu4 = round(otu4, 2)

dist_w_py = vegdist(otu4, method = "bray")
dist_w_py = dist.zeroes(otu4, dist_w_py)

#weighted distance analysis
h = hclust(dist_w_py, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

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

fwdsel_df = merge(sample_data(d4), rel_otu_pop.year, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,56:930], #OTUS#
                 fwdsel_df[,12:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_pop.year ~ CA_PCT + MN + ZN + P2 + K_PCT, data=fwdsel_df) ###this works for anova.cca
summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
library(plyr)
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Year"))
#ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)

Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "slategray4", "yellow4"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year)), size=2) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()
g

# RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_pop.year, sample_data(d4), by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)

fwdsel_df2 = subset(fwdsel_df, Year == 2 | Year == 3)

test=forward.sel(fwdsel_df[,56:930], #OTUS#
                 fwdsel_df[,36:55], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_pop.year ~ pg.sm1 + pg.st + pg.at, data=fwdsel_df) ###this works for anova.cca
summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")
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
rda_met <- read.xls("met.xlsx", sheet = 3,
                    verbose = TRUE, na.strings="N/A", perl="C:/Perl64/bin/perL")
row.names(rda_met) = rda_met[,1]
ccdata = merge(ccdata, rda_met, by = "row.names")

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
g

# MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("MG_PCT", "MN", "SAND", "ZN")
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

mrm.soil = lm(dist_w_py ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

mrm.soil.otus = lm(dist_w_py ~ dist(soil.otus))
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
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
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
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Soil OMF ----------------------------------------------------------------

d_f = merge_samples(d_r, "Population")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(gen_f$Row.names)
list = as.character(gen_f$rank)
gen_f = gen_f[,-c(1,8:15)]
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
d_s = subset_samples(d_s, Month == "Feb"| Month == "Apr")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####scale envt data according to above sample selection
met4 = as.data.frame(sample_data(d_s))
s_env_met = met4[,cbind(1,2,3,4,5,6,7,52,53)]
s_env = met4[,8:51]
s_env = scale(s_env)
met5 = merge(s_env_met, s_env, by = "row.names")
row.names(met5) = met5$Row.names
met5$int =  paste(met5$Population,".",met5$Month,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s

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
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d2 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d2))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = dist.zeroes(rel_otu_code, dist_w)

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d2)$Population, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Depth), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d2)$pop.year, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$int), permutations = 999)
a

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d3 = merge_samples(d_s, "pop.year")
otu3 = data.frame(otu_table(d3))
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
               cellnote = otu3, notecex=1.0,
               notecol="white")

# Soil OMF RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d3), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,56:84], #OTUS#
                 fwdsel_df[,12:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT + CLAY + B + MN + NO3_N + CU + MG + K + CA + CEC +
           NA_PCT + PH + S_SALTS + OM, data=fwdsel_df) ###this works for anova.cca

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT, data=fwdsel_df) ###this works for anova.cca

summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")

scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3" = "Year"))
Year = as.factor(ccdata$Year)

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "slategray4", "yellow", "springgreen2","violetred4"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_int, sample_data(d3), by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)

fwdsel_df2 = subset(fwdsel_df, Year == 2 | Year == 3)

test=forward.sel(fwdsel_df2[,2:27], #OTUS#
                 fwdsel_df2[,63:79], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(fwdsel_df2[,2:27] ~ pg.sm2 + dr.st + dr.sm2 + ar.sm2
         + pg.at, data=fwdsel_df2)
t = summary(cc)
t

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
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3"="Year"))
ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)
Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("CEC", "SAND", "B", "ZN", "NO3_N", "K")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d4)
names = c("pg.ppt", "smd.sm2", "smd.ppt", "pg.sm2", "smd.st", "pg.st")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

library(ecodist)

MRM(dist_w_py ~ dist(my.env2) + dist(my.soil2), nperm=1000)
summary(lm(dist_w_py ~ dist(my.env2) + dist(my.soil2)))

mrm.env = lm(dist_w_py ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil = lm(dist_w_py ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared


