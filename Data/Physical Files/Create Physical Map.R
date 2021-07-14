library(readr)
library(tidyr)
library(dplyr)
phy <- read_tsv("Ramu et al.- Positions.txt", col_names = TRUE) %>%
  filter(Chromosome==4)
phy <- slice(phy,-grep("SB",phy$`Locus name`)) %>%
  select(`Locus name`,`Physical Map Position`)
phy$`Locus name` <- gsub("X","",phy$`Locus name`) %>%
  tolower()
write.table(phy,"phy_r3.txt",sep = "\t",quote = FALSE, row.names = TRUE)