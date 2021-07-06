library(tidyr)
library(dplyr)
phy <- read.table("Ramu et al.- Positions.txt", sep="\t", header = TRUE, row.names = 1) 
phy <- filter(phy, Chromosome==1) %>%
  slice(-grep("SB",phy$Locus.name)) %>%
  select(c("Locus.name","Physical.Map.Position"))
phy[,1] <- gsub("X","",phy[,1])
write.table(phy,"phy_r.txt",sep = "\t",quote = FALSE, row.names = TRUE)