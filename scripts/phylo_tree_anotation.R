#setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")
#setwd("C:/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico/")
#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
#setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

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

x = read.mrbayes("results/phylo/tul_bayes/tul_otus_aln_trim.nexusi.con.tre")
#x = read.raxml("results/phylo/tul_raxml/RAxML_bipartitionsBranchLabels.bootFinal")

###for mrbayes
tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=prob, vjust=-.5, size=3))

tree$data$prob[is.na(tree$data$prob)] = ""

tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, size=1))

###for raxml 

x = read.raxml("results/phylo/tul_raxml/RAxML_bipartitionsBranchLabels.bootFinal")
tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=bootstrap, vjust=-.5, size=3))

x = read.raxml("results/phylo/cer_raxml/RAxML_bipartitionsBranchLabels.bootFinal")
tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label=bootstrap, vjust=-.5, size=3))

###import associated matrix

library(phyloseq)
#devtools::install_github("benjjneb/decontam")
library(decontam)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")

d
###subset root samples
decon.d
#decontamination
decon
####subsetetting for months and scaling of envt data
d_r
##merge samples
dpop = merge_samples(d_r, "Population")

otu = data.frame(t(otu_table(dpop)))
row.names(otu) = gsub("otu", "denovo", row.names(otu))
##make binary matrix
otu2 = as.data.frame(ifelse(otu == 0, 0, 1))

tree = ggtree(x, color="black", size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, hjust = 0.5, size=0.5))

g = gheatmap(tree, otu2, offset = 0.5, width=0.5, font.size=3, colnames_angle=-45, hjust = 1, 
             low = "grey", high = "black")
g

##missleneaous

ggtree(x) +
  geom_label(mapping = aes(label = node), size = 2)

exp = ggtree(x, aes(color=branch.length), size=1, linetype="dotted") + geom_tiplab(size=3, color="black") +
  geom_nodelab(aes(x=branch, label = round(as.numeric(prob), 2), vjust=-.5, hjust = 0.5, size=0.5)) +
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

