setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

setwd("C:/Users/jaspr/Google Drive/Metagenomics/pico_comb_run/pico/")

#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

install.packages("doParallel")
install.packages("DT")
install.packages("exactRankTests")
install.packages("foreach")
install.packages("ggplot2")
install.packages("Rcpp")
install.packages("shiny")
install.packages("coin")
install.packages("openxlsx")
install.packages("scripts/ancom.R.zip", repos = NULL, type="source")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

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

# Alpha diversity ---------------------------------------------------------

temp = estimate_richness(d_r)
temp = merge(met, temp, by = "row.names")

a = summary(aov(Shannon ~ Population + Year + Month + Stage + Pop_size + Demo, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad
a = summary(aov(Simpson ~ Population + Year + Month + Stage + Pop_size + Demo, data = temp))
a
p.ad = p.adjust(a[[1]]$`Pr(>F)`)
p.ad

plot_richness(d_r, x= "Population", measures=c("Shannon", "Simpson") )

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

#permanova with significant factors. 
a = adonis2(dist_w ~ sample_data(d.bd)$Pop_size + sample_data(d.bd)$Population + sample_data(d.bd)$Year + sample_data(d.bd)$Demo + sample_data(d.bd)$Stage, strata = "Pop_size", permutations = 999)
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
                                    taxa_are_rows = F), sample_data(d2))

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

####prune the otu table according to simper analyses results
#to retain otus which need to be shown in heatmap

hm.otus = unique(c(row.names(sim.df.popsize)[1:50], row.names(sim.df.demo)[1:50]))
otu.hm = otu3[,colnames(otu3) %in% hm.otus]
colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g1 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu.hm, notecex=0.7,
               notecol="white")

#Relative abundance plots ---------------------------------------------

d.ra = d_r

source("scripts/relative_abundance_plots.R")

###relative abundances at OTU level
p.otus

###relative abundances at phylum level
p.phy

###relative abundances at family level
p.fam


# MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("P2", "MN", "SAND", "ZN")
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

mrm.soil = lm(dist_w_int ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

mrm.soil.otus = lm(dist_w_py ~ dist(soil.otus))
summary(mrm.soil.otus)$adj.r.squared


# Association with soil OTUs ----------------------------------------------

net = read.delim("data/net.txt", sep = "\t", header = T)
row.names(net) = net[,1]
net = net[,-1]
ind = indpower(net) #columns shud be OTUs
ind.col = melt(ind)
ind.col.t.sel = ind.col[grepl("t.r", ind.col$Var2), ] #target selection
ind.fin.sel = ind.col.t.sel[grepl("i.s", ind.col.t.sel$Var1), ]
#ind.col.sel = ind.col[ind.col$value > 0.1,]
write.csv(ind.fin.sel, file = "ind.fin.sel.csv")
##first select roots target OTUs from Var2 column and then further subset
#dataframe for VAr1 to include only soil indicators
library(igraph)
x = graph_from_data_frame(ind.fin.sel) ##convert data frame to igraph object, first and second columns will be vertices and edges respectively, third and any other follwing columns will be edge attributes
write.graph(x, file="ind.fin.sel.txt",format="ncol") 


###select the specific vertix for plotting
neigh.nodes <- neighbors(x, V(x)$t.s15, mode="out")
vcol <- rep("grey40", vcount(x))
vcol[V(x)$t.s15] <- "gold"
vcol[neigh.nodes] <- "#ff9d00"
plot(x, vertex.color=vcol)
#Identify the edges going into or out of a vertex, for instance the WSJ. For a single node, use incident(), for multiple nodes use incident_edges()

