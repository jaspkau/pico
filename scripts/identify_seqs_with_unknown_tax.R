####first identify the most abundant OTUs in decontaminated phyloseq
##object

d_f = merge_samples(d_s, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
gen_f$rank = paste(as.character(gen_f$Row.names),"|",substr(gen_f$Family, 3, 5))
list = as.character(gen_f$rank)
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank", "otu", "pval", "p.ad")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:1000]
f = gen_f[,names(gen_f) %in% who]

####Identify the OTUs which have unidentified taxonomy at family level
f2 = f[, grepl("ide", names(f))]
names(f2) = gsub("\\|", "", names(f2))
names(f2) = gsub("ide", "", names(f2))
names(f2) = gsub("otu", "denovo", names(f2))
names(f2) = gsub(" ", "", names(f2))

############read fasta file
library(devtools)
#install_github("GuillemSalazar/FastaUtils")
library(FastaUtils)

fas = read.fasta(file = "data/chimera_filtered_rep_set.fasta", clean_name = FALSE)
fas2 = fas[fas$seq.name %in% names(f2),]

library(seqinr)

write.fasta(sequences=as.list(fas2$seq.text),names= fas2$seq.name,file.out="unidentified_OTUs.fasta")
