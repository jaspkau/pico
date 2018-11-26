setwd("C:/Users/jaspkaur/Google Drive/Metagenomics/pico_comb_run/pico/")

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

##show table for site and year combination for all the variables and then write p values at the end and
#make plots only for the sgnificant interactions

data <- read.xls("data/pico_comb_envt.xlsx",
               sheet=4,verbose=TRUE,na.strings="N/A")


keep = c("Site", "pop.size","year", "month", "day2", "rf")
#pg_env_avg[is.na(pg_env_avg)] <- " "
avg = data[ ,(names(data) %in% keep)]
avg = na.omit(avg)
avg$year = as.factor(avg$year)
avg$month = as.factor(avg$month)

#random = ~ 1 | rep/subject
#Indicates that each subject-within-rep unit will have its own intercept
model = lmerTest::lmer(as.numeric(rf) ~ pop.size + Site + year + (1|month/day2), data = avg)

a = anova(model)
a

a$p.ad = p.adjust(a$`Pr(>F)`, method = "bonferroni")
a$p.ad

############################### Site
posthoc <- glht(model, linfct = mcp(pop.size = "Tukey"))
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
