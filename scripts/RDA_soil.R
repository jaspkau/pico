setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico//")

#library(adespatial)  
library(phyloseq)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)
library(adespatial)
library(vegan)
library(splitstackshape)
library(plyr)

source("scripts/make_phyloseq.R")

d = subset_samples(d, Month == "Feb"| Month == "Apr")
d

sample_data(d)$int = paste(sample_data(d)$Source,".",sample_data(d)$Population,".",sample_data(d)$Year)
  
d.fin = subset_taxa(d, Family == "f:Ceratobasidiaceae"| 
                  Family == "f:Tulasnellaceae")

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
dhc <- as.dendrogram(h)
nodePar <- list(lab.cex = 1, pch = c(NA, 19), cex = 0.7, col = "blue")
p = plot(dhc,  xlab = "Weighted Bray-Curtis distance", nodePar = nodePar, horiz = TRUE)
p

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
test=forward.sel(fwdsel_df[,42:ncol(fwdsel_df)], #OTUS#
                 fwdsel_df[,17:35], #environmental variables#
                 nperm = 999, R2thresh = 0.9, adjR2thresh = 9, alpha = 1)
print(test) ###look at the results and select variables which are incrasing the

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

cc = rda(rel_otu_int ~ P1 + P2 + ZN + CA + S_SALTS +
           CU + CEC + PH + OM, data= fwdsel_df) ###this works for anova.cca

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