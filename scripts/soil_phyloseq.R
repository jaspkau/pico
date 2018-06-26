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
met4 = as.data.frame(sample_data(d_s))
s_env_met = met4[,cbind(1,2,3,4,5,6,7,8,9,10,11,38,39)]
s_env = met4[,12:37]
s_env = scale(s_env)
met5 = merge(s_env_met, s_env, by = "row.names")
row.names(met5) = met5$Row.names
met5$int =  paste(met5$Population,".",met5$Pop_size,".",met5$Demo,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s