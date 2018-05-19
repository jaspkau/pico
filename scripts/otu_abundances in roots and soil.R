setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")

source("scripts/make_phyloseq.R")

d.abun = merge_samples(d, "Source")
mytaxa = taxa_names(d_r)
d.abun2 = prune_taxa(mytaxa, d.abun)

otu3 = data.frame(otu_table(d.abun2))
otu3 = decostand(otu3, method = "hellinger")
otu3 = otu_table(otu3, taxa_are_rows = FALSE)
d.abun3 = merge_phyloseq(otu3, tax_table(d.abun2), sample_data(d.abun2))
test = data.frame(otu_table(d.abun3))
test2 = data.frame(t(test))
plot(test2)
label = row.names(test2)

p = ggplot(test2, aes(R, S)) + geom_point() + 
  geom_text(aes(label=label),size = 3, vjust = "inward", 
            hjust = "inward", check_overlap = TRUE)
p
                  