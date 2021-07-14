library(readr)
library(tidyr)
library(dplyr)

#Genetic_Files Data Loading, make sure map files are in the working directory.
files <- list.files()
i <- grep("map",files)
j <- 1
genetic <- list()
for(i in i){
  genetic[[j]] <- read.table(paste0(files[i]), sep="\t", skip=13, row.names = 1, fill = TRUE)
  j=j+1
}
remove(files,i,j)

# Consensus Map Loading
cons <- list()
cons[[3]] <- read_tsv("patt2003.txt", col_names = TRUE) %>%
  filter(map_name=="F") %>%
  select(c("feature_name","feature_start")) %>%
  # semi_join(tokeep, by=c("feature_name"="tokeep")) %>%
  arrange(feature_start)
cons[[3]]$feature_name <- tolower(cons[[3]]$feature_name)
cons[[3]]$feature_name <- gsub("[[:punct:]]","",cons[[3]]$feature_name)
write.table(cons[[1]],"klein3.txt",sep = "\t",quote = FALSE, row.names = TRUE)
