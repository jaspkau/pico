decon.d = subset_samples(d, Source == "S")
decon.d

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

decon

d_s = subset_samples(decon, Month == "Feb"| Month == "Apr")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####scale envt data according to above sample selection

sample_data(d_s)$int =  paste(sample_data(d_s)$Population,".",gsub("20", "", sample_data(d_s)$Year))
sample_data(d_s)$int = gsub(" ", "", sample_data(d_s)$int)
