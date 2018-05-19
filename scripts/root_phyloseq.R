decon.d = subset_samples(d, Source == "R")
decon.d

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

decon

d_r = subset_samples(decon, Month == "Feb"| Month == "Apr")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d_r))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,10,11,38,39)]
env = met2[,12:37]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d_r = merge_phyloseq(tax_table(d_r), otu_table(d_r), sample_data(met3))