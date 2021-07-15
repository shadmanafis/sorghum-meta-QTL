library(readr)
library(tidyr)
library(dplyr)
cons[[1]] <- read_tsv("kl2004 consensus.txt" ,col_names = TRUE) %>%
  filter(map_name=="LG-01") %>%
  select(c("feature_name","feature_aliases","feature_start")) %>%
  #semi_join(all, by=c("feature_name"="all")) %>%
  arrange("Mean")

cons[[2]] <- read_tsv("mac2010a consensus.txt" ,col_names = TRUE) %>%
  filter(map_name=="1") %>%
  select(c("feature_name","feature_aliases","feature_start")) %>%
  #semi_join(all, by=c("feature_name"="all")) %>%
  arrange("Mean")

cons[[3]] <- read_tsv("mac2011.txt.tsv" ,col_names = TRUE) %>%
  filter(LG=="1") %>%
  select(c("Marker", "Mean")) %>%
  #semi_join(all, by=c("feature_name"="all")) %>%
  arrange("Mean")

cons[[4]] <- read_tsv("patt2003-1 consensus.txt" ,col_names = TRUE) %>%
  filter(map_name=="A") %>%
  select(c("feature_name","feature_aliases","feature_start")) %>%
  #semi_join(all, by=c("feature_name"="all")) %>%
  arrange("Mean")

cons[[5]]<- read_tsv("Sanches_1.txt" ,col_names = TRUE) %>%
  filter(map_name=="A") %>%
  select(c("feature_name","feature_aliases","feature_start")) %>%
  #semi_join(all, by=c("feature_name"="all")) %>%
  arrange("Mean")

write.table(cons[[1]][,c(1,3)],"consensus1.txt",sep = "\t",quote = FALSE, row.names = TRUE)

write.table(cons[[1]],"consensus1.txt",sep = "\t",quote = FALSE, row.names = TRUE)
write.table(cons[[2]][,c(1,3)],"consensus2.txt",sep = "\t",quote = FALSE, row.names = TRUE)
write.table(cons[[3]][,c(1,3)],"consensus3.txt",sep = "\t",quote = FALSE, row.names = TRUE)
write.table(cons[[4]][,c(1,3)],"consensus4.txt",sep = "\t",quote = FALSE, row.names = TRUE)


