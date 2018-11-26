decon.d = subset_samples(d, Source == "R")
decon.d

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

decon

d_r = subset_samples(decon, Month == "Feb"| Month == "Apr")
d_r
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_r)
d_r

d.romf = subset_taxa(d_r, Family=="f:Serendipitaceae"| Family=="f:Sebacinaceae"| Family=="f:Thelephoraceae"| 
                    Family=="f:Ceratobasidiaceae"| Family=="f:Pezizaceae"| Family=="f:Tulasnellaceae"| 
                    Family=="f:Pyronemataceae"| Family=="f:Tuberaceae"| Family=="f:Agaricaceae"| 
                    Family=="f:Clavulinaceae"| Family=="f:Corticiaceae"| Family=="f:Inocybaceae"| 
                    Family=="f:Marasmiaceae"| Family=="f:Russulaceae"| Family=="f:Tricholomataceae"| 
                    Family=="f:Typhulaceae"| Family=="f:Physalacriaceae"| Family=="f:Psathyrellaceae"|
                    Family=="f:Hymenogastraceae"| Family=="f:Incertae sedis"| Family=="f:Hydnangiaceae"|
                      Family=="f:Hymenochaetaceae"| Family=="f:Atheliaceae"| Family=="f:Auriculariaceae"|
                      Family=="f:Cantharellaceae"| Family=="f:Sclerodermataceae"| Family=="f:Helvellaceae"|
                      Family=="f:Herpotrichiellaceae"| Family=="f:Nectriaceae")
