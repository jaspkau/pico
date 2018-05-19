setwd("C://Users//jaspkaur//Google Drive//Metagenomics//pico_comb_run//pico/")

#library(adespatial)  
library(phyloseq)
library(ancom.R)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

###ROOT OMF ANALYSIS......................................

source("scripts/make_phyloseq.R")

decon.d = subset_samples(d, Source == "R")
decon.d

####decontaminate phyloseq object based on frequency and prevelence

source("scripts/decontaminate_phyloseq.R")

decon

d_r = subset_samples(decon, Month == "Feb"| Month == "Apr")
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

