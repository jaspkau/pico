####decontaminate phyloseq object based on frequency and prevelence
df <- as.data.frame(sample_data(decon.d)) # Put sample_data into a ggplot-friendly d
df$LibrarySize <- sample_sums(decon.d)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

###combined method of decomtamination based on prevelanec and frequency
sample_data(decon.d)$is.neg <- sample_data(decon.d)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(decon.d, method="combined", neg="is.neg", conc="DNA_conc", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon <- prune_taxa(!contamdf.prev$contaminant, decon.d)
decon