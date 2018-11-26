d.rf = merge_samples(d_r, "Pop_size")
otu_rf = data.frame(otu_table(d.rf))

library(iNEXT)
otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000, 70000, 80000, 100000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g1 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g1


d.rf = merge_samples(d_r, "Population")
otu_rf = data.frame(otu_table(d.rf))

library(iNEXT)
otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000, 100000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g2 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g2

d.rf = merge_samples(d_r, "Year")
otu_rf = data.frame(otu_table(d.rf))

otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000, 100000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g3 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g3

d.rf = merge_samples(d_r, "Stage")
otu_rf = data.frame(otu_table(d.rf))

otu_rc = data.frame(t(otu_rf)) ####columns should be samples
m <- c(3000, 10000, 20000, 50000, 70000, 100000)
out = iNEXT(otu_rc, q=0, datatype="abundance", size=m, nboot = 100)
g = ggiNEXT(out, type=1, se = FALSE, facet.var="none")

g4 = g + scale_color_manual(values=c("wheat4", "violetred4", "turquoise3", "tomato2", "springgreen2",
                                     "slateblue2", "navyblue", "magenta", "blue2", "black", "seagreen4",
                                     "dodgerblue1", "orangered4", "yellow4", "slategray4", "olivedrab1","deeppink4", "aquamarine",
                                     "hotpink", "yellow1", "tan2", "red3", "pink1"))
g4

ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2)