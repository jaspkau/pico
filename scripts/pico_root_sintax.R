setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

setwd("C:/Users/jaspr/Google Drive/Metagenomics/pico_comb_run/pico/")

#source("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/packages.r")
setwd("/Users/administrator/Documents/jaspreet/pico/pico_comb_run/pico")

#library(adespatial)
library(phyloseq)

###ROOT OMF ANALYSIS......................................
# Make phyloseq object ----------------------------------------------------

otu <- read.delim(file = "data/97%/ITS2/otu_table_no_singletons_sintax.txt", 
                  sep = "\t", header = T)
otu = otu[,-ncol(otu)]
row.names(otu) = paste(gsub("denovo", "o", otu[,1]))
otu = otu[,-1]
#Rarefy(otu, depth = min(rowSums(otu)))
#otu = otu[,colSums(otu) > 0]
#site_list = colnames(otu)
otu_tab = otu_table(as.matrix(otu), taxa_are_rows = T)

###Format SINTAX taxonomy table

library(reshape2)

tax = read.delim(file = "data/97%/ITS2/tax.sintax", sep = "\t", header = F)
row.names(tax) = tax$V1
list = tax$V2
tax2 = colsplit(list, pattern ="\\(|\\),", c("Kingdom", "Kingdom_conf", "Phylum", "Phylum_conf", "Class", "Class_conf", "Order", "Order_conf", "Family", "Family_conf", "Genus", "Genus_conf", "Species", "Species_conf"))
tax2$Species_conf = gsub("\\)", "", tax2$Species_conf)
tax2$Species_conf = as.numeric(tax2$Species_conf)
tax2[is.na(tax2)] <- 0
row.names(tax2) = row.names(tax)

source("scripts/tax_func.R") #90% conf
level = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
for (i in 1:7){ ###i 1 = kingdom, i 2 = phylum etc
  tax2 = tax_assign(tax2, paste(level[i], "_conf", sep = ""), level[i])
}

tax2 = tax2[,-c(2,4,6,8,10,12,14)]

tax2$row = row.names(tax2)
tax2[,8:9] = colsplit(tax2$row, " ", c("otu", "seq"))
row.names(tax2) = paste(gsub("denovo", "o", tax2$row))
tax2 = tax2[,-c(8:9)]
tax2 = tax_table(as.matrix(tax2))

#meta data
library(readxl)
met <- as.data.frame(read_excel("data/met.xlsx", sheet = 1))
row.names(met) = met$Code
met$int = paste(met$Population,".",met$Pop_size,".",met$Demo,".", gsub("20", "", met$Year))
met$int = gsub(" ", "", met$int)
met$pop.year = paste(met$Population, ".", gsub("20", "", met$Year))
met$pop.year = gsub(" ", "", met$pop.year)

#phyloseq object

d = merge_phyloseq(tax2, otu_tab, sample_data(met))
d
d = subset_taxa(d, Kingdom == "d:Fungi")
d
decon.r = subset_samples(d, Source == "R")
decon.r

####decontaminate phyloseq object based on frequency and prevelence
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

df <- as.data.frame(sample_data(decon.r)) # Put sample_data into a ggplot-friendly d
df$LibrarySize <- sample_sums(decon.r)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

###prevelanec based
sample_data(decon.r)$is.neg <- sample_data(decon.r)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(decon.r, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
decon.r <- prune_taxa(!contamdf.prev$contaminant, decon.r)
decon.r

contamdf.freq <- isContaminant(decon.r, method="frequency", conc="DNA_conc")
table(contamdf.freq$contaminant)
which(contamdf.freq$contaminant == "TRUE")
decon.r <- prune_taxa(!contamdf.freq$contaminant, decon.r)
decon.r

d_r = subset_samples(decon.r, Month == "Feb"| Month == "Apr")
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
d_r

# Alpha diversity ---------------------------------------------------------

plot_richness(d_r, x= "Population", measures=c("Observed", "Shannon", "Simpson") )

plot_richness(d_r, x= "Year", measures=c("Observed", "Shannon", "Simpson"))

plot_richness(d_r, x= "Month", measures=c("Observed", "Shannon", "Simpson"))

plot_richness(d_r, x= "Stage", measures=c("Observed", "Shannon", "Simpson"))

plot_richness(d_r, x= "Pop_size", measures=c("Observed", "Shannon", "Simpson"))

plot_richness(d_r, x= "Demo", measures=c("Observed", "Shannon", "Simpson"))

temp = estimate_richness(d_r)
temp = merge(met, temp, by = "row.names")

a = summary(aov(Observed ~ Population, data = temp))
a
a = summary(aov(Simpson ~ Population, data = temp))
a
a = summary(aov(Shannon ~ Population, data = temp))
a
a = summary(aov(Observed ~ Year, data = temp))
a
a = summary(aov(Simpson ~ Year, data = temp))
a
a = summary(aov(Shannon ~ Year, data = temp))
a
a = summary(aov(Observed ~ Month, data = temp))
a
a = summary(aov(Simpson ~ Month, data = temp))
a
a = summary(aov(Shannon ~ Month, data = temp))
a
a = summary(aov(Observed ~ Stage, data = temp))
a
a = summary(aov(Simpson ~ Stage, data = temp))
a
a = summary(aov(Shannon ~ Stage, data = temp))
a
a = summary(aov(Observed ~ Pop_size, data = temp))
a
a = summary(aov(Simpson ~ Pop_size, data = temp))
a
a = summary(aov(Shannon ~ Pop_size, data = temp))
a
a = summary(aov(Observed ~ Demo, data = temp))
a
a = summary(aov(Simpson ~ Demo, data = temp))
a
a = summary(aov(Shannon ~ Demo, data = temp))
a

# Beta diversity with bray ------------------------------------------------
library(vegan)
otu2 = t(data.frame(otu_table(d_r)))
otu2 = decostand(otu2, method = "hellinger")
rowSums(otu2)
otu2 = otu_table(as.matrix(otu2), taxa_are_rows = F)

d3 = merge_phyloseq(tax2, otu2, sample_data(d))
rel_otu_code = data.frame(otu_table(d3))

dist_w = vegdist(rel_otu_code, method = "bray")

###PERMANOVA

###Weighted distance

a = adonis(dist_w ~ sample_data(d3)$Population, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Stage, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Month, permutations = 999)
a
a = adonis(dist_w ~ as.factor(sample_data(d3)$Year), permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Pop_size, permutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Demo, pe?indrmutations = 999)
a
a = adonis(dist_w ~ sample_data(d3)$Population*as.factor(sample_data(d3)$Year), permutations = 999)
a

#SIMPER (similaritypercentage) analyses ----------------------------------------------

ind.df = data.frame(otu2)##the taxa should be columns and this otu table is hellinger tranfromed

#Identification of species most responsible for differences among groups of samples
#SIMPER(similaritypercentage), Based on abundance, does not weigh occurrence frequency as indicator species analysis does.

sim = simper(ind.df, sample_data(d3)$Pop_size)
sim.sum = summary(sim)
sim.df.popsize = data.frame(sim.sum$S_L)

sim.popsize.otus = row.names(sim.df.popsize)[1:50]

library(dplyr)

mann.popsize.df = ind.df[,names(ind.df) %in% sim.popsize.otus]

mann.popsize.df2 = merge(mann.popsize.df, met3, by = "row.names")
mann.popsize.df2$Pop_size = as.factor(mann.popsize.df2$Pop_size)

sim.kw.popsize = c()
for(i in 2:51){
  column = names(mann.popsize.df2[i])
  k = kruskal.test(mann.popsize.df2[,i]~Pop_size, data = mann.popsize.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k)))
  sim.kw.popsize = rbind(sim.kw.popsize, results)
} 

sim.kw.popsize$p.ad = p.adjust(sim.kw.popsize$pval, method = "bonferroni")

############OTUs differenting between demopgraphic classes

sim = simper(ind.df, sample_data(d3)$Demo)
sim.sum = summary(sim)
sim.df.demo = data.frame(sim.sum$F_NF)

sim.demo.otus = row.names(sim.df.demo)[1:50]

library(dplyr)

mann.demo.df = ind.df[,names(ind.df) %in% sim.demo.otus]

mann.demo.df2 = merge(mann.demo.df, met3, by = "row.names")
mann.demo.df2$Demo = as.factor(mann.demo.df2$Demo)

##Do kruskal wallis test with the first OTUs from simper analyses

sim.kw.demo = c()
for(i in 2:51){
  column = names(mann.demo.df2[i])
  k.demo = kruskal.test(mann.demo.df2[,i]~Demo, data = mann.demo.df2)$"p.value"
  results = data.frame(otu = paste(column), pval = as.numeric(paste(k.demo)))
  sim.kw.demo = rbind(sim.kw.demo, results)
} 

sim.kw.demo$p.ad = p.adjust(sim.kw.demo$pval, method = "bonferroni")

# Hierarchial clustering --------------------------------------------------

#compressing the phyloseq object at level which is significantly different

d2 = merge_samples(d_r, "int")
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

#weighted distance analysis
h = hclust(dist_w_int, method = "average")
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

####prune the otu table according to simper analyses results
#to retain otus which need to be shown in heatmap

hm.otus = unique(c(row.names(sim.df.popsize)[1:50], row.names(sim.df.demo)[1:50]))
otu.hm = otu3[,colnames(otu3) %in% hm.otus]
colfunc <- colorRampPalette(c("grey", "black"))
library(gplots)

g1 = heatmap.2(as.matrix(otu.hm), 
               Rowv = as.dendrogram(h), margins = c(10, 10), col = colfunc(100), 
               xlab = "Weighted Bray Curtis dissimilarity distances",
               trace = "none",
               cellnote = otu.hm, notecex=0.7,
               notecol="white")

# Realtive abundance plots at OTU level ------------------------------------------------

d_f = merge_samples(d_r, "int")
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
who = names(sort(colMeans(gen_f), decreasing = TRUE))[1:50]
f = gen_f[,names(gen_f) %in% who]
f$Other = 1-rowSums(f)
who = c(who, "Other")
dd = f
dd$sl = row.names(dd)
m = melt(dd, id.vars = c("sl"), measure.vars = who)
library(RColorBrewer)
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 5, name = "Blues"),brewer.pal(n = 10, name = "Paired"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "gold1", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 11, name = "Spectral")))

library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

# Realtive abundance plots at Family level ------------------------------------------------

d_f = tax_glom(d_r, taxrank = "Family")
d_f = merge_samples(d_f, "int")
gen_f = data.frame(otu_table(d_f))
gen_f = t(gen_f)
gen_f = merge(gen_f, tax_table(d_f), by = "row.names")
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
state_col2 = scale_fill_manual(name = "State3", values=c(brewer.pal(n = 3, name = "Pastel1"), "azure3", "burlywood1", "cornflowerblue", "wheat4", "cyan4", "turquoise3", "hotpink", "tan2", 
                                                         "springgreen2", "slateblue2", "red3", "navyblue", 
                                                         "magenta", "olivedrab1", "blue2", "black", "yellow1",
                                                         "dodgerblue1", "orangered4", "yellow4", "deeppink4", 
                                                         "slategray4", "seagreen4" , "aquamarine",
                                                         "tomato2", brewer.pal(n = 8, name = "Accent")))
library(scales)

p = ggplot(m, aes(sl, fill = variable)) + geom_bar(aes(weight = value)) + 
  theme_bw(base_size = 20) + state_col2 + theme(axis.text.x = element_text(angle = 0, hjust=.5, size = 12)) +
  xlab("Sample") + ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  theme(legend.text = element_text(face = "italic")) + guides(fill = guide_legend(ncol = 1, reverse=T))+ scale_y_continuous(labels = percent_format())
p$data$variable = factor(p$data$variable, ordered = TRUE, levels = rev(who))
p

#ggsave(file="jc.treatment.nms.jpg")

###......................................................

# Environmnetal data comparisons ------------------------------------------

require(partykit)

##final will store final results

final = c()
for(i in 10:29){
  column = names(met2[i])
  f = anova(aov(met2[,i]~ Population, data=met2))$"F value"
  av = anova(aov(met2[,i]~ Population, data=met2))$"Pr(>F)"
  results = data.frame(otu = paste(column), F.value = paste(f), Pr..F. = paste(av))
  final = rbind(final, results)
} 

write.csv(final, file='soil_aov.csv')

##need to find variables which can dicriminate between groups

soil <- as.data.frame(read_excel("data/met.xlsx", sheet = 2))

library(rpart)
con = rpart.control(cp = 0.01, 
                    maxcompete = 6)
fit <- rpart(Population ~ OM + P1 + P2 + PH +	K +	MG + CA +	
                               CEC + NO3_N + S + ZN	+ MN + FE +	CU + B +	S__SALTS
                          + SAND + SILT + CLAY, method="class", control = con, data = soil)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

#install.packages("rpart.plot")
library(rpart.plot)

rpart.plot(fit)

###environmental data

env <- as.data.frame(read_excel("data/met.xlsx", sheet = 3))

library(lme4)
install.packages("lmerTest")
library(lmerTest)

env$int = paste(env$Population,".",env$Year)
p = ggplot(env, aes(int, rf)) + geom_point(aes(color = month))
p 
model = lmer(atemp ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Population + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Pop_size + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ Demo + (1|month),
             data=env,
             REML=TRUE)
anova(model)

model = lmer(atemp ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(stemp ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rf ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)
model = lmer(rh ~ as.factor(Year) + (1|month),
             data=env,
             REML=TRUE)
anova(model)

# RDA with soil ---------------------------------------------------------------------

fwdsel_df = merge(sample_data(d4), rel_otu_int, by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)
test=forward.sel(fwdsel_df[,41:ncol(fwdsel_df)], #OTUS#
                 fwdsel_df[,15:34], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~ P2 + ZN + MN + SAND + CLAY + S + S_SALTS + FE, data=fwdsel_df) ###this works for anova.cca
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
library(plyr)
ccdata = rename(ccdata, c("site_1"="Population", "site_2"="Pop_size", "site_3"="Demo", "site_4"="Year"))
#ccdata = rename(ccdata, Population = site_1, Month = site_2, Year = site_3)

Year = as.factor(ccdata$Year)
state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red", "blue", "magenta",
                                            "slategray4", "yellow1"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year)), size= Pop_size) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()
g

# RDA with environment ---------------------------------------------------------------------

fwdsel_df = merge(rel_otu_int, sample_data(d4), by = "row.names")
row.names(fwdsel_df) = fwdsel_df[,1]
names(fwdsel_df)

fwdsel_df2 = subset(fwdsel_df, Year == 2 | Year == 3)

test=forward.sel(fwdsel_df[,56:930], #OTUS#
                 fwdsel_df[,36:55], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the
#AdjR2Cim and whose p value id <0.05

cc = rda(rel_otu_int ~ pg. + pg.st + pg.at, data=fwdsel_df) ###this works for anova.cca
summary(cc)
anova.cca(cc)
anova.cca(cc, by = "axis")
anova.cca(cc, by = "terms")
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
rda_met <- read.xls("met.xlsx", sheet = 3,
                    verbose = TRUE, na.strings="N/A", perl="C:/Perl64/bin/perL")
row.names(rda_met) = rda_met[,1]
ccdata = merge(ccdata, rda_met, by = "row.names")

state_col_ord = scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                            "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                            "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                            "hotpink", "yellow1", "tan2", "red3", "pink1"))

state_col_ord = scale_color_manual(values=c("black", "red"))

g = ggplot(data = ccdata, aes(RDA1 , RDA2)) + 
  geom_point(aes(color = Population, shape = as.factor(Year)), size=2) + state_col_ord + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian()

g = g + geom_segment(data=arrowdata,aes(x=0,xend=RDA1,y=0,yend=RDA2),
                     arrow = arrow(length = unit(0.01,"npc")), colour="black") + 
  geom_text(data=arrowdata,aes(x=RDA1,y=RDA2,label= variables),size=3)+
  coord_cartesian() + theme_bw()
g

# MRM and varition partioning ---------------------------------------------

my.soil = sample_data(d4)
names = c("P2", "MN", "SAND", "ZN")
my.soil2 = my.soil[,names] ####order of rows should be same as community data matrix

my.env = sample_data(d4)
names = c("pg.sm2", "dr.st", "pg.st")
my.env2 = my.env[,names] ####order of rows should be same as community data matrix

soil.otus = gen_f[1:6,] ###need to create soil phyloseq object first

library(ecodist)

MRM(dist_w_py ~ dist(my.env2) + dist(my.soil2) + dist(soil.otus), nperm=1000)
summary(lm(dist_w_py ~ dist(my.env2) + dist(my.soil2) + dist(soil.otus)))

mrm.env = lm(dist_w_py ~ dist(my.env2))
summary(mrm.env)$adj.r.squared

mrm.soil = lm(dist_w_int ~ dist(my.soil2))
summary(mrm.soil)$adj.r.squared

mrm.soil.otus = lm(dist_w_py ~ dist(soil.otus))
summary(mrm.soil.otus)$adj.r.squared


# Association with soil OTUs ----------------------------------------------
#Use indpower function
ind = indpower(mann.popsize.df)
ind.col = melt(ind)
ind.col.sel = ind.col[ind.col$value > 0.4,]

##then select Var1 only for important OTUs from roots to select for soil OTUs