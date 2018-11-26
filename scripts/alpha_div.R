aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met2, aldiv, by = "row.names")
row.names(temp) = temp[,1]

g1 <- ggplot(temp, aes(x=Pop_size, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g1

aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met2, aldiv, by = "row.names")
row.names(temp) = temp[,1]

g2 <- ggplot(temp, aes(x=Population, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g2


aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met2, aldiv, by = "row.names")
row.names(temp) = temp[,1]

g3 <- ggplot(temp, aes(x=Year, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g3

aldiv = estimate_richness(d_r, measures = c("Shannon", "Simpson"))
temp = merge(met2, aldiv, by = "row.names")
row.names(temp) = temp[,1]

g4 <- ggplot(temp, aes(x=Stage, y=Simpson)) + 
  geom_boxplot(aes(fill= "slategray4")) + 
  labs(x = paste("Site"), 
       y = paste("Simpson diversity index (H)")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g4

ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2)