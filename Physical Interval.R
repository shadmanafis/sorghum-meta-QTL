library(readr)
library(tidyr)
library(dplyr)

# Load Genetic Data

phy <- read_tsv("Ramu et al.- Positions.txt", col_names = TRUE) %>%
  filter(Chromosome==1)
phy <- slice(phy,-grep("SB",phy$`Locus name`)) %>%
  select(`Locus name`,`Physical Map Position`)
phy$`Locus name` <- gsub("X","",phy$`Locus name`) %>%
  tolower()

# Set WOrking Directory

Cons_map <- read.table("map_Cons-qtls_map.txt", sep="\t", skip=13, row.names = 1, fill = TRUE)

left_join(Cons_map,phy,by=c("V2"="Locus name"))
