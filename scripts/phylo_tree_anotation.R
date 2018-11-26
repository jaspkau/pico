#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//oab/")
#setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#setwd("/Users/administrator/Documents/jaspreet/pico/pico/")
#setwd("/Users/administrator/Desktop/oab")

#https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/library("ape")
library("Biostrings")
library("ggplot2")
#devtools::install_github("GuangchuangYu/treeio")
#devtools::install_github("GuangchuangYu/ggtree")
library(treeio)
library(ape)
library(ggtree)

#https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html
#https://guangchuangyu.github.io/presentation/2016-ggtree-chinar/

x = read.mrbayes("results/phylo/root_soil/cer_raxml/RAxML_bipartitions.bootFinal")
#x = read.raxml("results/phylo/tul_raxml/RAxML_bipartitionsBranchLabels.bootFinal")

###for mrbayes
tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=prob, vjust=-.5, size=3))

tree$data$prob[is.na(tree$data$prob)] = ""

tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, size=1))

###for raxml 

x = read.raxml("results/phylo/root_soil/tul_raxml/RAxML_bipartitionsBranchLabels.bootFinal")
tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=bootstrap, vjust=-.5, size=3))

###remove the <50 bootstraps
tree$data$bootstrap = ifelse(tree$data$bootstrap < 50, paste(""), tree$data$bootstrap)

###import associated matrix

library(phyloseq)
library(reshape2)
library(readxl)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

#########read phyloseq object

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
source("scripts/soil_phyloseq.R")

d.comb = merge_phyloseq(d_r, d_s)
d.comb

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

taxa_names(d.fin) = gsub("otu", "denovo", taxa_names(d.fin))

sample_data(d.fin)$int = paste(sample_data(d.fin)$Source,".",sample_data(d.fin)$Pop_size)
sample_data(d.fin)$int = gsub(" ", "", sample_data(d.fin)$int)

##merge samples
dpop = merge_samples(d.fin, "int")

####for abudance based associated matrix
otu = data.frame(t(otu_table(dpop)))
row.names(otu) = gsub("otu", "denovo", row.names(otu))

otu2 = as.data.frame(t(otu))
otu2 = as.data.frame(t(decostand(otu2, method = "hellinger")))

##make binary matrix
otu2 = as.data.frame(ifelse(otu == 0, 0, 1))
###

tree = ggtree(x, color="black", size=0.7) + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = bootstrap), vjust=-.5, hjust = 0.5, size=3)

tree$data$bootstrap = ifelse(tree$data$bootstrap < 50, paste(""), tree$data$bootstrap)

g = gheatmap(tree, otu2, offset = 0.0, width=0.5, font.size=3, colnames_angle=-45, hjust = 0, 
             low = "grey97", high = "black")
g

g = gheatmap(tree, otu2, offset = 0.0, width=0.5, font.size=3, colnames_angle=-45, hjust = 0, 
             low = "red3", high = "limegreen")
g

##missleneaous

ggtree(x) +
  geom_label(mapping = aes(label = node), size = 2)
get
exp = ggtree(x, aes(color=branch.length), size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, hjust = 0.8, size=0.5)) +
  theme(legend.position="bottom")

cp = collapse(exp, node = 183)
cp + geom_point2(aes(subset=(node == 183)), size=5, shape=23, fill="steelblue")

####Collapse nodes

#Branch Length = consider an alignment of length 100, and RAxML 
#estimates a particular branch length to be 0.1,  
#the total number of substitutions on that branch should be 0.1 X 100
# collapse if prob < 0.6 and branch legnth is < 0.00

nod.claps = ifelse((tree$data$branch.length < 0.03)==TRUE, paste(tree$data$node), paste(""))
exp %>% collapse(nod.claps)

