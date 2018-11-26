#setwd("/Users/administrator/Desktop/pico")
setwd("/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico")
#setwd("/Users/jas/Google Drive/Metagenomics/pico_comb_run/pico")

library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(vegan)

source("scripts/make_phyloseq.R")
source("scripts/root_phyloseq.R")
source("scripts/soil_phyloseq.R")

d.comb = merge_phyloseq(d_r, d_s)
d.comb

d.comb = subset_samples(d.comb, Pop_size != "Con")

d.comb = prune_taxa(taxa_sums(d.comb) > 0, d.comb)

seq_dep = estimate_richness(d.comb, measures = "Observed")
seq_dep$sequences = sample_sums(d.comb)
seq_dep = merge(sample_data(d.comb), seq_dep, by = "row.names")

sample_data(d.comb)$int = paste(sample_data(d.comb)$Source,".",sample_data(d.comb)$Population,".",sample_data(d.comb)$Year)

d.fin = subset_taxa(d.comb, Family == "f:Ceratobasidiaceae"| 
                      Family == "f:Tulasnellaceae")
d.fin

####scale envt data according to above sample selection
met2 = data.frame(sample_data(d.fin))
env_met = met2[,cbind(1,2,3,4,5,6,7,8,9,10,11,37,38,39)]
env = met2[,12:36]
env = scale(env)
met3 = merge(env_met, env, by = "row.names")
row.names(met3) = met3$Row.names

d.fin = merge_phyloseq(tax_table(d.fin), otu_table(d.fin), sample_data(met3))

###identify the overlap between root and soil Otus
d.x = subset_taxa(d.fin, Family == "f:Tulasnellaceae")
d.x
d.r = subset_samples(d.fin, Source == "R")
d.r
d.r = prune_taxa(taxa_sums(d.r) >= 1, d.r)
d.r
d.s.x = subset_samples(d.fin, Source == "S")
d.s.x
d.s.x = prune_taxa(taxa_sums(d.s.x) >= 1, d.s.x)
d.s.x
d.r.s = prune_taxa(!(taxa_names(d.s.x) %in% taxa_names(d.r)), d.s.x)
d.r.s

d.r.s2 = prune_taxa(taxa_names(d.r) %in% taxa_names(d.r.s), d.r)
d.r.s2
alltaxa = taxa_names(d.r)
d.r.s2 <- alltaxa[!(alltaxa %in% taxa_names(d.r.s))]
ex1 = prune_taxa(d.r.s2, d.r)
taxa_sums(d.r)
taxa_sums(ex1)
d.r.s3 = tax_glom(ex1, "Family")
summary(taxa_sums(d.r.s3))
summary(taxa_sums(ex1))

# Rarefaction -------------------------------------------------------------

sample_data(d.comb)$int = paste(sample_data(d.comb)$Source,".",sample_data(d.comb)$Stage)

d.rf = merge_samples(d_r, "Stage")
otu_rf = data.frame(otu_table(d.rf))

library(iNEXT)
otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(0, 3000, 10000, 20000, 50000, 100000, 110000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g1 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g1

# Alpha diversity ---------------------------------------------------------

aldiv = estimate_richness(d.fin, measures = c("Shannon", "Simpson"))
temp = merge(sample_data(d.fin), aldiv, by = "row.names")
temp = temp[,-1]
row.names(temp) = temp[,1]

temp$ef = 1/(1-temp$Simpson)

# The conversion of Shannon diversity to effective numbers is exp(H)
temp$ef.sha = exp(temp$Shannon)

#####Comparisons with catergorical variables

####Use shannon diversity index coz simpson is inflating diversity in samples with 0 seqs
shapiro.test(temp$ef.sha)

alpha.kw = c()
for(i in c(8)){
  column = names(temp[i])
  k.demo = kruskal.test(ef.sha ~ as.factor(temp[,i]), data = temp)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  alpha.kw = rbind(alpha.kw, results)
}

alpha.kw$p.ad = p.adjust(alpha.kw$pval, method = "bonferroni")
alpha.kw

avg = temp %>%
  group_by(Source) %>%
  summarise(simp = mean(ef.sha))
avg

bp <- ggplot(temp, aes(x=Source, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bp

# Soil OMF Beta diversity with bray ------------------------------------------------

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d.fin)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
#site_list = colnames(otu)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d.bd = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d.bd))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = stepacross(dist_w)

###PERMANOVA

###Weighted distance

a = adonis2(dist_w ~ sample_data(d.bd)$Source, permutations = 999)
a
###ajust P-values
p.adjust(a$`Pr(>F)`, method = "bonferroni")

# Hierarchial clustering --------------------------------------------------
#compressing the phyloseq object at level which is significantly different

d2 = merge_samples(d.fin, "int")
otu3 = data.frame(otu_table(d2))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = otu3
rowSums(otu3)
otu3 = round(otu3, 2)

dist_w_int = vegdist(otu3, method = "bray")
dist_w_int[is.na(dist_w_int)] <- 0

otu3_tab = otu_table(as.matrix(otu3), taxa_are_rows = F)
d4 = merge_phyloseq(tax2, otu_table(as.matrix(otu3_tab), 
                                    taxa_are_rows = F), sample_data(d2))

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h) %>% set("labels_cex", 0.5)
ggd1 <- as.ggdend(dhc)

p1 = ggplot(ggd1, horiz = TRUE,theme = theme_minimal())
p1

d.an = d_r

source("scripts/differential_abundance.R")

####prune the otu table according to simper analyses results
#to retain otus which need to be shown in heatmap

###detected through ANCOM
an.otus = anc$detected

otu.hm.an = otu3[,colnames(otu3) %in% an.otus]
write.csv2(otu.hm.an, file = "data/root_ancom_otus.csv")

otu.hm = merge(t(otu.hm.an), tax_table(d.hc), by = "row.names")
#gen_f = merge(gen_f, sim.kw.popsize, by.x = "Row.names", by.y = "otu")
otu.hm$rank = paste(as.character(otu.hm$Row.names),"|",substr(otu.hm$Family, 3, 5))
row.names(otu.hm) = otu.hm$rank
drops <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
otu.hm = otu.hm[ , !(names(otu.hm) %in% drops)]
otu.hm = data.frame(t(otu.hm))

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g1 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu.hm, notecex=0.4,
               notecol="white")

# PCoA between root and soil OMF ------------------------------------------
#####Weighted

pc = capscale(dist_w_int ~ 1, comm = rel_otu_int) ###~ means function of nothing
pc$CA$eig
s = summary(pc)
cap = data.frame(pc$CA$u)
plot(cap[,1], cap[,2])
cap$site = row.names(cap)
#install.packages("splitstackshape")
library(splitstackshape)
cap = cSplit(cap, "site", ".")
cap = plyr::rename(cap, c("site_1"= "Source", "site_2"= "Population", "site_3" = "Year"))
cap$popL = cap$Population=="SCE"|cap$Population=="PLF"
#cap$Population = cap$Row.names
label = cap$Population

state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "chartreuse3", "cyan2", "darkorange2", "coral4"))

g = ggplot(data = cap, aes(MDS1 , MDS2, color = Population)) + 
  geom_point(aes(color = Population, shape = Source, size = 3)) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") + stat_ellipse(aes(MDS1 , MDS2, group = popL)) +
  coord_cartesian() +
  labs(x = paste("Axis 1 (", round(s$cont$importance[2,1]*100, digits =2), "%)", sep = ''),
       y = paste("Axis 2 (", round(s$cont$importance[2,2]*100, digits =2), "%)", sep = ''))

# RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d4), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
library(adespatial)
test=forward.sel(fwdsel_df[,41:ncol(fwdsel_df)], #OTUS#
                 fwdsel_df[,13:31], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the

cc = rda(rel_otu_int ~ ZN + P2 + NA. + PH + OM, data= fwdsel_df) ###this works for anova.cca

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
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = plyr::rename(ccdata, c("site_1"="Source", "site_2"= "Population", "site_3" = "Year"))
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

# ANCOM -------------------------------------------------------------------
d.an = d.fin
ancom.otu = t(data.frame(otu_table(d.an))) ##columns = OTUs and should be counts
ancom.otu = merge(ancom.otu, sample_data(d.an), by = "row.names")
row.names(ancom.otu) = ancom.otu$Code
ancom.otu = ancom.otu[,-1]
names(ancom.otu)
##look for the grouping variable you want to use
ancom.fin = ancom.otu[, grepl("otu", names(ancom.otu))|grepl("Source", names(ancom.otu))]

anc = ANCOM(ancom.fin, sig = 0.05, multcorr = 1)
anc$detected
plot_ancom(anc)


