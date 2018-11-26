a = prune_taxa(taxa_names(d_r) %in% anc$detected, d_r)
b = subset_samples(a, Pop_size == "L")
c = merge_samples(b, "Pop_size")
e = tax_glom(c, "Family")

taxa_sums(e)

keep = c("otu9916")
a = prune_taxa(taxa_names(d_r) %in% keep, d_r)
b = merge_samples(a, "Population")
c = merge_samples(b, "Pop_size")
e = tax_glom(c, "Family")

sample_sums(b)

a = subset_taxa(d_s, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
e = tax_glom(a, "Family")
taxa_sums(e)