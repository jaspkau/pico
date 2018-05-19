setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

#library(adespatial)  
library(phyloseq)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(adespatial)
library(vegan)
library(splitstackshape)
library(plyr)

###ROOT OMF ANALYSIS......................................

source("scripts/make_phyloseq.R")

decon.d = subset_samples(d, Source == "R")
decon.d

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

d_fin = subset_samples(decon, Month == "Feb"| Month == "Apr")
d_fin
d_fin = prune_taxa(taxa_sums(d_fin) >= 1, d_fin)
d_fin

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d_fin))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,10,11,38,39)]
env = met2[,12:37]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d_fin = merge_phyloseq(tax_table(d_fin), otu_table(d_fin), sample_data(met3))
d_fin

####most abundant OTUs in roots
d_r = subset_samples(d_fin, Source == "R")
d_r = prune_taxa(taxa_sums(d_r) >= 1, d_fin)
d_r
d_f = merge_samples(d_r, "Population")
r.otus = data.frame(otu_table(d_f))
r.otus = r.otus/rowSums(r.otus)
r.otus.sel = names(sort(colMeans(r.otus), decreasing = TRUE))[1:25]

###

otu_rda = data.frame(otu_table(d_fin))
otu_rda2 = otu_rda[r.otus.sel,]
otu_rda2 = otu_rda2[,colSums(otu_rda2) > 0]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_rda3 = otu_table(as.matrix(otu_rda2), taxa_are_rows = T)

d_fin2 = merge_phyloseq(tax_table(d_fin), otu_rda3, sample_data(met3))
d_fin2

# Hierarchial clustering --------------------------------------------------

#compressing the phyloseq object at level which is significantly different

d2 = merge_samples(d_fin2, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")

otu3_tab = otu_table(as.matrix(otu3), taxa_are_rows = F)
d4 = merge_phyloseq(tax2, otu_table(as.matrix(otu3_tab), 
                                    taxa_are_rows = F), sample_data(d2))

# RDA with soil ---------------------------------------------------------------------
fwdsel_df = merge(sample_data(d4), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,42:ncol(fwdsel_df)], #OTUS#
                 fwdsel_df[,16:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05
cc = rda(rel_otu_int ~ P2 + ZN + MN + SAND + PH + CA + OM, data=fwdsel_df) ###this works for anova.cca
summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
#install.packages("splitstackshape")
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Source", "site_2"="Population", "site_3"="Pop_size", "site_4"="Demo", "site_5"="Year"))
#ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)

Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "chartreuse3", "cyan2", "darkorange2", "coral4"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = Source, size = 3)) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()
g
