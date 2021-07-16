library(readr)
library(tidyr)
library(dplyr)
library(GRanges)
# Load Genetic Data

phy <- read_tsv("Ramu et al.- Positions.txt", col_names = TRUE) %>%
  filter(Chromosome==1)
phy <- slice(phy,-grep("SB",phy$`Locus name`)) %>%
  select(`Locus name`,`Physical Map Position`)
phy$`Locus name` <- gsub("X","",phy$`Locus name`) %>%
  tolower()

# Set WOrking Directory, Granges list from all the consensus maps

files <- list.files()
i <- grep("map",files)
j <- 1
Consensus <- data.frame()
for(i in i){
  Chromosome <- strsplit(strsplit(files[i], split = " ")[[1]][2],split = ".txt")[[1]]
  f <- read.table(paste0(files[i]), sep="\t", skip=13, row.names = 1, fill = TRUE)
  f$Chromosome <- rep(Chromosome,dim(f)[1])
  Consensus <- rbind(Consensus,f)
  j=j+1
 
}
remove(files,i,j,f, Chromosome)
Consensus <- data.frame(chr=Consensus$Chromosome,start=Consensus$V3,end=Consensus$V3,feature=Consensus$V2)
GR <- makeGRangesFromDataFrame(Consensus, keep.extra.columns=TRUE)

# Get markers inside or flanking the Meta-QTL

hits <- subjectHits(findOverlaps(gr,GR))
GR[(min(hits)-1):(max(hits)+1)]$feature