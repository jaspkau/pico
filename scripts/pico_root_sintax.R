setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico (1)/")

#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

library(dunn.test)
library(adespatial)
library(phyloseq)
library(metagenomeSeq)
library(mixOmics)

###ROOT OMF ANALYSIS......................................
# Make phyloseq object ----------------------------------------------------

otu <- read.delim(file = "data/95%/otu_table_no_singletons_sintax.txt", 
                  sep = "\t", header = T)
otu = otu[,-ncol(otu)]
row.names(otu) = otu[,1]
otu = otu[,-1]
#Rarefy(otu, depth = min(rowSums(otu)))
otu = otu[,colSums(otu) > 0]
site_list = colnames(otu)
otu_tab = otu_table(as.matrix(otu), taxa_are_rows = T)

###Format SINTAX taxonomy table

library(reshape2)

tax = read.delim(file = "data/95%/tax.sintax", sep = "\t", header = F)
row.names(tax) = tax$V1
list = tax$V2
tax2 = colsplit(list, pattern ="\\(|\\),", c("Kingdom", "Kingdom_conf", "Phylum", "Phylum_conf", "Class", "Class_conf", "Order", "Order_conf", "Family", "Family_conf", "Genus", "Genus_conf", "Species", "Species_conf"))
tax2$Species_conf = gsub("\\)", "", tax2$Species_conf)
tax2$Species_conf = as.numeric(tax2$Species_conf)
tax2[is.na(tax2)] <- 0
row.names(tax2) = row.names(tax)

source("scripts/tax_func.R") #90% conf
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
library(readxl)
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))
row.names(met) = met$Code
met$int = paste(met$Population,".",met$Pop_size,".",met$Demo,".", gsub("20", "", met$Year))
met$int = gsub(" ", "", met$int)
met$int2 = paste(met$Pop_size,".",met$Demo)
met$int2 = gsub(" ", "", met$int2)
met$pop.year = paste(met$Population, ".", gsub("20", "", met$Year))
met$pop.year = gsub(" ", "", met$pop.year)

#phyloseq object

d = merge_phyloseq(tax2, otu_tab, sample_data(met))
d
d = subset_taxa(d, Kingdom == "d:Fungi")
d
#d_t = subset_taxa(d, Family == "f:Tulasnellaceae")
#t = data.frame(taxa_names(d_t))
#write.csv(t, file = "tul_otu.csv")
d_r = subset_samples(d, Source == "R")
d_r
d_r = subset_samples(d_r, Month == "Feb"| Month == "Apr")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d_r))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,40,41)]
env = met2[,10:39]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d_r = merge_phyloseq(tax_table(d_r), otu_table(d_r), sample_data(met3))
d_r

# Alpha diversity ---------------------------------------------------------

temp = estimate_richness(d_r)
temp = merge(met, temp, by = "row.names")
p =  ggplot(temp, aes(Population, Chao1))+ geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

a = summary(aov(Observed ~ Population, data = temp))
a
a = summary(aov(Simpson ~ Population, data = temp))
a
a = summary(aov(Shannon ~ Population, data = temp))
a
a = summary(aov(Observed ~ Year, data = temp))
a
a = summary(aov(Simpson ~ Year, data = temp))
a
a = summary(aov(Shannon ~ Year, data = temp))
a

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d3)$Population, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Stage, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Month, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Pop_size, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Demo, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Population*as.factor(sample_data(d3)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Population*sample_data(d3)$Pop_size, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Population*sample_data(d3)$Demo, permutations = 999)
a

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different

d2 = merge_samples(d_r, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
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

# Venn Diagrams -----------------------------------------------------------

install.packages("VennDiagram")
library(VennDiagram)
library(gplots)
plf.otus = subset_samples(d_r, Population == "PLF")
plf.otus = prune_taxa(taxa_sums(plf.otus) >= 1, plf.otus)
plf.otus = taxa_names(plf.otus)

ple.otus = subset_samples(d_r, Population == "PLE")
ple.otus = prune_taxa(taxa_sums(ple.otus) >= 1, ple.otus)
ple.otus = taxa_names(ple.otus)

sce.otus = subset_samples(d_r, Population == "SCE")
sce.otus = prune_taxa(taxa_sums(sce.otus) >= 1, sce.otus)
sce.otus = taxa_names(sce.otus)

scw.otus = subset_samples(d_r, Population == "SCW")
scw.otus = prune_taxa(taxa_sums(scw.otus) >= 1, scw.otus)
scw.otus = taxa_names(scw.otus)

ch.otus = subset_samples(d_r, Population == "CH")
ch.otus = prune_taxa(taxa_sums(ch.otus) >= 1, ch.otus)
ch.otus = taxa_names(ch.otus)

mx.otus = subset_samples(d_r, Population == "MX")
mx.otus = prune_taxa(taxa_sums(mx.otus) >= 1, mx.otus)
mx.otus = taxa_names(mx.otus)

venn(list(plf.otus, ple.otus, sce.otus, scw.otus, ch.otus, mx.otus))

venn(list(plf.otus, sce.otus))

venn(list(plf.otus, scw.otus))

venn(list(plf.otus, ple.otus))

venn(list(plf.otus, mx.otus))

# Indicator species analyses ----------------------------------------------

library(indicspecies)
ind.df = data.frame(otu2)##the taxa should be columns and this otu table is hellinger tranfromed

#Identification of species most responsible for differences among groups of samples
#SIMPER(similaritypercentage), Based on abundance, does not weigh occurrence frequency as indicator species analysis does.

sim = simper(ind.df, sample_data(d3)$Pop_size)
sim.sum = summary(sim)
sim.df.popsize = data.frame(sim.sum$S_L)

sim = simper(ind.df, sample_data(d3)$Demo)
sim.sum = summary(sim)
sim.df.demo = data.frame(sim.sum$F_NF)

l.mann.otus = unique(c(row.names(sim.df.popsize)[1:35], row.names(sim.df.demo)[1:35]))

library(dplyr)

mann.df = ind.df[,names(ind.df) %in% l.mann.otus]

mann.df2 = merge(mann.df, met3, by = "row.names")
mann.df2$Pop_size = as.factor(mann.df2$Pop_size)
mann.df2$Demo = as.factor(mann.df2$Demo)

##Do kruskal wallis test with the first OTUs from simper analyses

sim.kw = c()
for(i in 2:36){
  column = names(mann.df2[i])
  k = kruskal.test(mann.df2[,i]~Pop_size, data = mann.df2)$"p.value"
  k.demo = kruskal.test(mann.df2[,i]~Demo, data = mann.df2)$"p.value"
  results = data.frame(otu = paste(column), pvalue.popsize = paste(k), pvalue.demo = paste(k.demo))
  sim.kw = rbind(sim.kw, results)
} 

##fdr correction
##embed p values in relative abundance graph

# Network analyses --------------------------------------------------------

library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

d.otu = data.frame(otu_tab) ##it should no non-normalized
otu.net = d.otu[row.names(d.otu) %in% l.mann.otus,]
otu_tab.net = otu_table(as.matrix(otu.net), taxa_are_rows = T)
d.net = merge_phyloseq(otu_tab.net, tax2, sample_data(d3))

spiec.out=spiec.easi(d.net, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(d.net)))
write.graph(spiec.graph,file="spieceasi.ncol.txt",format="ncol") 
plot_network(spiec.graph, d.net, type='taxa', color = "Family", label=NULL)

clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]

betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 

otu.ids=colnames(spiec.out$data)
edges=E(spiec.graph)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(spiec.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"forestgreen")
  }else if(beta<0){
    edge.colors=append(edge.colors,"red")
  }
}
E(spiec.graph)$color=edge.colors

spiec.graph.b=spiec.graph
nodenames=V(spiec.graph.b)$name
V(spiec.graph.b)$name=getTaxonomy(nodenames, tax2, useRownames=TRUE)
E(spiec.graph.b)$arrow.size=5
V(spiec.graph.b)$color="white"
V(spiec.graph.b)$frame.color="black"
tkplot(spiec.graph.b)

# Realtive abundance plots at OTU level ------------------------------------------------

d_f = merge_samples(d_r, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = as.character(gen_f$Row.names)
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:35]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 3, name = "Pastel1"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 8, name = "Accent")))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_r, taxrank = "Family")
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
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 3, name = "Pastel1"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 8, name = "Accent")))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

#ggsave(file="jc.treatment.nms.jpg")

###......................................................

# Environmnetal data comparisons ------------------------------------------

require(partykit)

##final will store final results

final = c()
for(i in 10:33){
  column = names(met2[i])
  f = anova(aov(met2[,i]~ Pop_size, data=met2))$"F value"
  av = anova(aov(met2[,i]~ Pop_size, data=met2))$"Pr(>F)"
  results = data.frame(otu = paste(column), F.value = paste(f), Pr..F. = paste(av))
  final = rbind(final, results)
} 

write.csv(final, file='q1_table.csv')

##need to find variables which can dicriminate between groups



###environmental data

env <- as.data.frame(read_excel("data/met.xlsx", sheet = 3))

library(lme4)
install.packages("lmerTest")
library(lmerTest)

env$int = paste(env$Population,".",env$Year)
p = ggplot(env, aes(int, rf)) + geom_point(aes(color = month))
p 
model = lmer(atemp ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)

# RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d4), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,44:ncol(fwdsel_df)], #OTUS#
                 fwdsel_df[,14:37], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~ P2 + ZN + MN + CLAY + K, data=fwdsel_df) ###this works for anova.cca
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

