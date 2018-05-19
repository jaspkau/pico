setwd("C:/Users/jaspr/Google Drive/Metagenomics/pico_comb_run/pico/")

#library(adespatial)  
library(phyloseq)
library(ancom.R)

###ROOT OMF ANALYSIS......................................
# Make phyloseq object ----------------------------------------------------

otu <- read.delim(file = "data/97%/ITS2/otu_table_no_singletons_sintax.txt", 
                  sep = "\t", header = T)
otu = otu[,-ncol(otu)]
row.names(otu) = paste(gsub("denovo", "o", otu[,1]))
otu = otu[,-1]
otu_tab = otu_table(as.matrix(otu), taxa_are_rows = T)

###Format SINTAX taxonomy table

library(reshape2)

tax = read.delim(file = "data/97%/ITS2/tax.sintax", sep = "\t", header = F)
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
row.names(tax2) = paste(gsub("denovo", "o", tax2$row))
tax2 = tax2[,-c(8:9)]
tax2 = tax_table(as.matrix(tax2))

#meta data
library(readxl)
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))
row.names(met) = met$Code
met$int = paste(met$Population,".",met$Pop_size,".",met$Demo,".", gsub("20", "", met$Year))
met$int = gsub(" ", "", met$int)
met$pop.year = paste(met$Population, ".", gsub("20", "", met$Year))
met$pop.year = gsub(" ", "", met$pop.year)
met$bip.int = paste(met$Population,".",met$Source,".",gsub("20", "", met$Year))
met$bip.int = gsub(" ", "", met$bip.int)

#phyloseq object

d = merge_phyloseq(tax2, otu_tab, sample_data(met))
d
d = subset_taxa(d, Kingdom == "d:Fungi")
d

####decontaminate phyloseq object based on frequency and prevelence
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

df <- as.data.frame(sample_data(d)) # Put sample_data into a ggplot-friendly d
df$LibrarySize <- sample_sums(d)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

###decontaminate based on combined prevelanec and frequency methods
sample_data(d)$is.neg <- sample_data(d)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(d, conc="DNA_conc", method="combined", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon <- prune_taxa(!contamdf.prev$contaminant, d)
decon

d.sel = subset_samples(decon, Month == "Feb"| Month == "Apr")
d.sel
d.sel2 = prune_taxa(taxa_sums(d.sel) >= 1, d.sel)
d.sel2

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d.sel2))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,10,11,38,39,40)]
env = met2[,12:37]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d.fin = merge_phyloseq(tax_table(d.sel2), otu_table(d.sel2), sample_data(met3))
d.fin

##"Remove taxa not seen more than 1 times in at least 2% of the samples. 
#This protects against an OTU with small mean & trivially large C.V.
d.fin2 = filter_taxa(d.fin, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
d.fin2

# make bipartite graph ----------------------------------------------------
library(bipartite)
d.bip = merge_samples(d.fin2, "bip.int")
otu.bip = data.frame(otu_table(d.bip))
plotweb(otu.bip)
visweb(otu.bip) 
networklevel(otu.bip)
specieslevel(otu.bip)

#Use circlize 
#https://stackoverflow.com/questions/27500041/r-make-circle-chord-diagram-with-circlize-from-dataframe?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
