# Network analyses --------------------------------------------------------

library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

#install.packages("rPython")
library(rPython)

setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

source("scripts/pico_root_phyloseq_object_sintax.R")

##make root otu file
##"Remove taxa not seen more than 3 times in at least 5% of the samples. 
#This protects against an OTU with small mean & trivially large C.V.
#d.r.net = filter_taxa(d_r, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
d.r.net = merge_samples(d_r, "int")
otu.r.net = data.frame(otu_table(d.r.net)) ##it should be non-normalized
###slelect OTUs first from which you want to build network, rows shud be samples
#otu.r.net = otu.r.net[row.names(otu.r.net) %in% row.names(sim.kw.popsize),]
taxa_names(d.r.net) = paste(gsub("o", "r", taxa_names(d.r.net)))
colnames(otu.r.net) = paste(gsub("o", "r", colnames(otu.r.net)))

d.s = subset_samples(d, Population == "PLF"|Population == "PLE"
                         |Population == "SCW"| Population == "SCE"
                         |Population == "MX"| Population == "CH")

##make soil otu file
d.s.net = subset_samples(d.s, Source == "S")
d.s.net = subset_samples(d.s.net, Month == "Feb"| Month == "Apr")
d.s.net = prune_taxa(taxa_sums(d.s.net) >= 1, d.s.net)
d.s.net
d.s.net = filter_taxa(d.s.net, function(x) sum(x > 1) > (0.05*length(x)), TRUE)
d.s.net = merge_samples(d.s.net, "int")
d.s.net
otu.s.net = data.frame(otu_table(d.s.net)) ##it should no non-normalized
taxa_names(d.s.net) = paste(gsub("o", "s", taxa_names(d.s.net)))
colnames(otu.s.net) = paste(gsub("o", "s", colnames(otu.s.net)))

d.net = merge_phyloseq(d.r.net, d.s.net)

###Sparcc
spar = sparcc(data.frame(otu_table(d.net)))
sparcc.graph <- abs(spar$Cor) >= 0.3
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
ig.sparcc <- adj2igraph(sparcc.graph)
write.graph(ig.sparcc,file="ig.sparcc.txt",format="ncol") 

vsize <- rowMeans(clr(data.frame(otu_table(d.net)), 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

plot(ig.sparcc, vertex.label=NA, main="sparcc")

####SpiecEasi

spiec.out=spiec.easi(d.net, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(d.net)))
write.graph(spiec.graph,file="spieceasi.ncol.txt",format="ncol") 
plot_network(spiec.graph, d.net, type='taxa', color = "Family", label="value")

clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]

#Compute matrix of regression coefficients
betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))

positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 

#SPIEC-EASI with positive and negative edge colors
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

#SPIEC-EASI clustering with positive edges only

otu.ids=colnames(spiec.out$data)
edges=E(spiec.graph)
filtered.edges=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(spiec.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta<0){
    filtered.edges=c(filtered.edges,edges[e.index])
  }
}
spiec.graph.pos=delete_edges(spiec.graph, filtered.edges)
write.graph(spiec.graph.pos,file="spieceasi.ncol.pos.txt",format="ncol") 
