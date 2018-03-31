setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

source("scripts/pico_root_phyloseq_object_sintax.R")

# Soil: Analyses of total fungal community ----------------------------------------------------------------

###create soil phyloseq object

decon.s = subset_samples(d, Source == "S")
decon.s

####decontaminate phyloseq object based on frequency and prevelence

df <- as.data.frame(sample_data(decon.s)) # Put sample_data into a ggplot-friendly d
df$LibrarySize <- sample_sums(decon.s)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Samp_con)) + geom_point()

###prevelanec based
sample_data(decon.s)$is.neg <- sample_data(decon.s)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(decon.s, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
decon.s <- prune_taxa(!contamdf.prev$contaminant, decon.s)
decon.s

contamdf.freq <- isContaminant(decon.s, method="frequency", conc="DNA_conc")
table(contamdf.freq$contaminant)
which(contamdf.freq$contaminant)
decon.s <- prune_taxa(!contamdf.freq$contaminant, decon.s)
decon.s

d_s = subset_samples(decon.s, Month == "Feb"| Month == "Apr")
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
met5$int =  paste(met5$Population,".",met5$Month,".",met5$Pop_size,".",met5$Demo,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d2 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d2))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d2)$Population, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Depth), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Year), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Pop_size), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Demo), permutations = 999)
a

#SIMPER (similaritypercentage) analyses ----------------------------------------------

ind.df = data.frame(otu2)##the taxa should be columns and this otu table is hellinger tranfromed

#Identification of species most responsible for differences among groups of samples
#SIMPER(similaritypercentage), Based on abundance, does not weigh occurrence frequency as indicator species analysis does.

sim = simper(ind.df, sample_data(d2)$Pop_size)
sim.sum = summary(sim)
sim.df.popsize = data.frame(sim.sum$L_S)

sim.popsize.otus = row.names(sim.df.popsize)[1:1000]

library(dplyr)

mann.popsize.df = ind.df[,names(ind.df) %in% sim.popsize.otus]

mann.popsize.df2 = merge(mann.popsize.df, met5, by = "row.names")
mann.popsize.df2$Pop_size = as.factor(mann.popsize.df2$Pop_size)

sim.kw.popsize = c()
for(i in 2:1001){
  column = names(mann.popsize.df2[i])
  k = kruskal.test(mann.popsize.df2[,i]~Pop_size, data = mann.popsize.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k)))
  sim.kw.popsize = rbind(sim.kw.popsize, results)
} 

sim.kw.popsize$p.ad = p.adjust(sim.kw.popsize$pval, method = "bonferroni")

############OTUs differenting between demopgraphic classes

sim = simper(ind.df, sample_data(d2)$Demo)
sim.sum = summary(sim)
sim.df.demo = data.frame(sim.sum$F_NF)

sim.demo.otus = row.names(sim.df.demo)[1:1000]

library(dplyr)

mann.demo.df = ind.df[,names(ind.df) %in% sim.demo.otus]

mann.demo.df2 = merge(mann.demo.df, met5, by = "row.names")
mann.demo.df2$Demo = as.factor(mann.demo.df2$Demo)

##Do kruskal wallis test with the first OTUs from simper analyses

sim.kw.demo = c()
for(i in 2:1000){
  column = names(mann.demo.df2[i])
  k.demo = kruskal.test(mann.demo.df2[,i]~Demo, data = mann.demo.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  sim.kw.demo = rbind(sim.kw.demo, results)
} 

sim.kw.demo$p.ad = p.adjust(sim.kw.demo$pval, method = "bonferroni")

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d3 = merge_samples(d_s, "int")
otu3 = data.frame(otu_table(d3))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = data.frame(otu3)
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")
#dist_w_int = dist.zeroes(otu3, dist_w_int)

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(10, 3), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")

# Realtive abundance plots ------------------------------------------------
d_f = tax_glom(d_s, taxrank = "Family")
d_f = merge_samples(d_f, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
#gen_f$rank = paste(gen_f$Family)
gen_f$rank = as.character(gen_f$Family)
#gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
gen_f$rank = ifelse(gen_f$Phylum == "unidentified", paste(as.character(gen_f$Kingdom), as.character(gen_f$Phylum), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Phylum != "unidentified" &  gen_f$Class == "unidentified", paste(as.character(gen_f$Phylum), as.character(gen_f$Class), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Class != "unidentified" &  gen_f$Order == "unidentified", paste(as.character(gen_f$Class), as.character(gen_f$Order), sep = ";"), gen_f$rank)
gen_f$rank = ifelse(gen_f$Order != "unidentified" &  gen_f$Family == "unidentified", paste(as.character(gen_f$Order), as.character(gen_f$Family), sep = ";"), gen_f$rank)
list = as.character(gen_f$rank)
list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# ###SOIL OMF ANALYSIS WITH OTUs identified in ROOTS ----------------------------------------------------------------

dr = data.frame(otu_table(d_r))
dr = t(dr)
who = colnames(dr)

###create soil phyloseq object by using above info

otu_s = data.frame(otu_table(d))
otu_s2 = otu_s[who,]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_s = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_s
d_s = subset_samples(d_s, Source == "S")
d_s
d_s = subset_samples(d_s, Month == "Feb"| Month == "Apr")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####scale envt data according to above sample selection
met4 = as.data.frame(sample_data(d_s))
s_env_met = met4[,cbind(1,2,3,4,5,6,7,52,53)]
s_env = met4[,8:51]
s_env = scale(s_env)
met5 = merge(s_env_met, s_env, by = "row.names")
row.names(met5) = met5$Row.names
met5$int =  paste(met5$Population,".",met5$Month,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s

# Realtive abundance plots ------------------------------------------------
#d_f = tax_glom(d_r, taxrank = "Family")
d_f = merge_samples(d_s, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
#gen_f$rank = paste(gen_f$Family)
list = as.character(gen_f$rank)
#list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d2 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d2))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = stepacross(dist_w, path = "extended")

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d2)$Population, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Depth), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d2)$pop.year, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$int), permutations = 999)
a

###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d3 = merge_samples(d_s, "pop.year")
otu3 = data.frame(otu_table(d3))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = data.frame(otu3)
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")
dist_w_int = stepacross(dist_w_int)

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(10, 3), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")

# Soil OMF RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d3), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,56:84], #OTUS#
                 fwdsel_df[,12:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT + CLAY + B + MN + NO3_N + CU + MG + K + CA + CEC +
           NA_PCT + PH + S_SALTS + OM, data=fwdsel_df) ###this works for anova.cca

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT, data=fwdsel_df) ###this works for anova.cca

summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")

scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3" = "Year"))
Year = as.factor(ccdata$Year)

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "slategray4", "yellow", "springgreen2","violetred4"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_int, sample_data(d3), by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)

fwdsel_df2 = subset(fwdsel_df, Year == 2 | Year == 3)

test=forward.sel(fwdsel_df2[,2:27], #OTUS#
                 fwdsel_df2[,63:79], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(fwdsel_df2[,2:27] ~ pg.sm2 + dr.st + dr.sm2 + ar.sm2
         + pg.at, data=fwdsel_df2)
t = summary(cc)
t

#ZN + FE + MG + P2 + OM + MN + CEC + B + K + SC + NO3_N + P1 +SM + CU + NA. + CLAY + SILT + S + CA + PH, data=fwdsel_df) ###this works for anova.cca
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)
arrowdata
#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3"="Year"))
ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)
Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("CEC", "SAND", "B", "ZN", "NO3_N", "K")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d4)
names = c("pg.ppt", "smd.sm2", "smd.ppt", "pg.sm2", "smd.st", "pg.st")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

library(ecodist)

MRM(dist_w_py ~ dist(my.env2) + dist(my.soil2), nperm=1000)
summary(lm(dist_w_py ~ dist(my.env2) + dist(my.soil2)))

mrm.env = lm(dist_w_py ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil = lm(dist_w_py ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

# ###SOIL OMF ANALYSIS WITH 35 MOST abundant root OTUS----------------------------------------------------------------

d_f = merge_samples(d_r, "Population")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(gen_f$Row.names)
list = as.character(gen_f$rank)
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:35]

###create soil phyloseq object by using above info

otu_s = data.frame(otu_table(d))
otu_s2 = otu_s[who,]
otu_s2 = otu_s2[,colSums(otu_s2) > 0]
#otu_s3 = otu_s2[, colSums(otu_s2 > 0)]
otu_s4 = otu_table(as.matrix(otu_s2), taxa_are_rows = T)

d_s = merge_phyloseq(tax2, otu_s4, sample_data(met))
d_s
d_s = subset_samples(d_s, Source == "S")
d_s
d_s = subset_samples(d_s, Month == "Feb"| Month == "Apr")
d_s
d_s = prune_taxa(taxa_sums(d_s) >= 1, d_s)
d_s

####scale envt data according to above sample selection
met4 = as.data.frame(sample_data(d_s))
s_env_met = met4[,cbind(1,2,3,4,5,6,7,8,9,36,37,38)]
s_env = met4[,10:35]
s_env = scale(s_env)
met5 = merge(s_env_met, s_env, by = "row.names")
row.names(met5) = met5$Row.names
met5$int =  paste(met5$Population,".",met5$Month,".",met5$Pop_size,".",met5$Demo,".",gsub("20", "", met5$Year))
met5$int = gsub(" ", "", met5$int)

d_s = merge_phyloseq(tax_table(d_s), otu_table(d_s), sample_data(met5))
d_s

# Realtive abundance plots ------------------------------------------------
#d_f = tax_glom(d_r, taxrank = "Family")
d_f = merge_samples(d_s, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
gen_f$rank = paste(as.character(gen_f$Row.names), "_", gen_f$Family)
#gen_f$rank = paste(gen_f$Family)
list = as.character(gen_f$rank)
#list = paste(list, "_", rep(1:length(list)), sep = "")
gen_f = gen_f[,-1]
drops <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank")
gen_f = gen_f[ , !(names(gen_f) %in% drops)]
gen_f = data.frame(t(gen_f))
gen_f = gen_f/rowSums(gen_f)
names(gen_f) = list
#met$Sample = ordered(met$Sample, levels = c("A", "B", "C", "D", "E", "F", "G"))
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:25]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c("azure3", "burlywood1", "coral2", "wheat4", "violetred4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", "pink1", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2"))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Soil OMF Beta diversity with bray ------------------------------------------------

#New OTU table with relative abundances 

#First find the variable which is making the difference for weighted and unweighted data
otu2 = t(data.frame(otu_table(d_s)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu2[rowSums(otu2) > 0,]
#site_list = colnames(otu)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d2 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d2))

dist_w = vegdist(rel_otu_code, method = "bray")
dist_w = stepacross(dist_w)

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d2)$Population, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Depth), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Month), permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d2)$Year), permutations = 999)
a


###Do the hierarchial clustering by compressing the
#phyloseq object at level which is significantly different
library(BiodiversityR)

d3 = merge_samples(d_s, "int")
otu3 = data.frame(otu_table(d3))
otu3 = decostand(otu3, method = "hellinger")
rel_otu_int = data.frame(otu3)
rowSums(otu3)
otu3 = round(otu3, 2)

dist_uw_int = vegdist(otu3, method = "bray", binary = TRUE)
dist_w_int = vegdist(otu3, method = "bray")
dist_w_int = dist.zeroes(otu3, dist_w_int)

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)
g1 = heatmap.2(as.matrix(otu3), 
               Rowv = as.dendrogram(h), margins = c(10, 3), col = colfunc(50), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu3, notecex=1.0,
               notecol="white")

# Soil OMF RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d3), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,56:84], #OTUS#
                 fwdsel_df[,12:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT + CLAY + B + MN + NO3_N + CU + MG + K + CA + CEC +
           NA_PCT + PH + S_SALTS + OM, data=fwdsel_df) ###this works for anova.cca

cc = rda(rel_otu_int ~  P2 + FE + P1 + NA. + ZN + SILT + MG_PCT + S + 
           CA_PCT + K_PCT, data=fwdsel_df) ###this works for anova.cca

summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")

scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)

#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3" = "Year"))
Year = as.factor(ccdata$Year)

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "slategray4", "yellow", "springgreen2","violetred4"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_int, sample_data(d3), by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)

fwdsel_df2 = subset(fwdsel_df, Year == 2 | Year == 3)

test=forward.sel(fwdsel_df2[,2:27], #OTUS#
                 fwdsel_df2[,63:79], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(fwdsel_df2[,2:27] ~ pg.sm2 + dr.st + dr.sm2 + ar.sm2
         + pg.at, data=fwdsel_df2)
t = summary(cc)
t

#ZN + FE + MG + P2 + OM + MN + CEC + B + K + SC + NO3_N + P1 +SM + CU + NA. + CLAY + SILT + S + CA + PH, data=fwdsel_df) ###this works for anova.cca
scrs<-scores(cc,display="bp")
arrowdata<- data.frame(scrs)
arrowdata$variables <-rownames(arrowdata)
arrowdata
#mul<-vegan:::ordiArrowMul(scrs)
#mul
#arrowdata<- data.frame(scrs*mul)
#arrowdata$variables <-rownames(arrowdata)

ccdata = as.data.frame(scores(cc)$sites)
ccdata$site = row.names(ccdata)
library(splitstackshape)
ccdata = cSplit(ccdata, "site", ".")
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Month", "site_3"="Year"))
ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)
Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year), size = Month)) + state_col_ord + 
  scale_size_manual(values=c(2,4)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()

# Soil OMF MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("CEC", "SAND", "B", "ZN", "NO3_N", "K")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d4)
names = c("pg.ppt", "smd.sm2", "smd.ppt", "pg.sm2", "smd.st", "pg.st")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

library(ecodist)

MRM(dist_w_py ~ dist(my.env2) + dist(my.soil2), nperm=1000)
summary(lm(dist_w_py ~ dist(my.env2) + dist(my.soil2)))

mrm.env = lm(dist_w_py ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil = lm(dist_w_py ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

###......................................................

###......................................................


