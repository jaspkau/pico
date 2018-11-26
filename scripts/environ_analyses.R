setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

library(readxl)
library(vegan)
library(lubridate)#workwithdates
library(dplyr)#datamanipulation(filter,summarize,mutate)
library(ggplot2)#graphics
library(gridExtra)#tileseveralplotsnexttoeachother
#install.packages("tseries")
library(tseries)
#install.packages("forecast")
library(forecast)
library(lmtest)
library(gdata)
library(dendextend)
library(ggdendro)

library(TSclust)
library(gdata)
library(imputeTS)
library(lsmeans)
library(multcomp)

# Environmnetal data comparisons ------------------------------------------

require(partykit)

##final will store final results

soil <- as.data.frame(read_excel("data/met.xlsx", sheet = 2))
soil = subset(soil, Month == "Feb"|Month == "April")

final = c()
for(i in 6:24){
  column = names(soil[i])
  k.year = kruskal.test(soil[,i]~ as.factor(soil$Year), data = soil)$"p.value"
  k.pop = kruskal.test(soil[,i]~ as.factor(soil$Population), data = soil)$"p.value"
  k.popsize = kruskal.test(soil[,i]~ as.factor(soil$Pop_size), data = soil)$"p.value"
  results = data.frame(otu = paste(column), pval.year = as.numeric(paste(k.year)), pval.pop = as.numeric(paste(k.pop)), pval.popsize = as.numeric(paste(k.popsize)))
  final = rbind(final, results)
} 

final$pad.year = p.adjust(final$pval.year, method = "bonferroni")
final$pad.pop = p.adjust(final$pval.pop, method = "bonferroni")
final$pad.popsize = p.adjust(final$pval.popsize, method = "bonferroni")

write.csv(final, file='results/soil_kw.csv')

library(FSA) ###for multiple pairwise comparisons after kw
DT = dunnTest(soil$SILT ~ as.factor(soil$Population), data=soil, method="bonferroni")

library(rcompanion) ##for compact letter display after dunn
DT = DT$res
p.ad = cldList(P.adj ~ Comparison, data = DT, threshold = 0.05)
p.ad

###environmental data

data <- read.xls("data/pico_comb_envt.xlsx",
                 sheet=4,verbose=TRUE,na.strings="N/A")

keep = c("Site", "pop.size","year", "month", "day2", "st")
#pg_env_avg[is.na(pg_env_avg)] <- " "
avg = data[ ,(names(data) %in% keep)]
avg = na.omit(avg)
avg$year = as.factor(avg$year)
avg$month = as.factor(avg$month)

#random = ~ 1 | rep/subject
#Indicates that each subject-within-rep unit will have its own intercept
model =  lmerTest::lmer(as.numeric(st) ~ pop.size + Site + year + (1|month/day2), data = avg)

a = anova(model)
a

a$p.ad = p.adjust(a$`Pr(>F)`, method = "bonferroni")
a$p.ad

#############################Year
posthoc <- glht(model, linfct = mcp(Site = "Tukey"))
summary(posthoc)

cld(posthoc,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")

#############################Year
posthoc <- glht(model, linfct = mcp(year = "Tukey"))
summary(posthoc)

cld(posthoc,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")

