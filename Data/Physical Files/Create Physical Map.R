library(readr)
library(tidyr)
library(dplyr)
phy <- read.table("Ramu et al.- Positions.txt", sep="\t", header = TRUE, row.names = 1) 
phy <- filter(phy, Chromosome==1) %>%
  slice(-grep("SB",phy$Locus.name)) %>%
  select(c("Locus.name","Physical.Map.Position"))
phy[,1] <- gsub("X","",phy[,1])
write.table(phy,"phy_r.txt",sep = "\t",quote = FALSE, row.names = TRUE)

#read genetic files
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Data/Genetic Maps/Chromosome A")
files <- list.files()
i <- grep("map",files)
j <- 1
genetic <- list()
for(i in i){
  genetic[[j]] <- read.table(paste0(files[i]), sep="\t", skip=13, row.names = 1, fill = TRUE)
  j=j+1
}
remove(files,i,j)

# Modified maps from full consensus
cons[[5]] <- read_tsv("mac2010a consensus.txt", col_names = TRUE) %>%
  filter(map_name==1) %>%
  select(c("feature_name","feature_aliases","feature_start")) %>%
  semi_join(tokeep, by=c("feature_name"="tokeep")) %>%
  arrange("Mean")
write.table(cons,"sanches.txtt",sep = "\t",quote = FALSE, row.names = TRUE)


nms <- combn( names(genetic) , 2 , FUN = paste0 , collapse = "" , simplify = FALSE )
# Make the combinations of list elements
ll <- combn( l , 2 , simplify = FALSE )
# Intersect the list elements
out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
# Output with names
setNames( out , nms )