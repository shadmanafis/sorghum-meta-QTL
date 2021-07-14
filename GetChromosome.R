library(readr)
library(tidyr)
library(dplyr)
# Create a list of marker name in respective chromosomes

L <- rep(list(NULL),10)
Ldat <- function(file_name,order){
  file <- read_tsv(file_name, col_names = TRUE)
  chr <- unique(file$map_name)
  chr <- chr[is.na(chr)==FALSE]
  for (i in 1:10){
    x <- filter(file,map_name==chr[order[i]]) %>%
      select("feature_name")
    x <- tolower(x$feature_name)
    x <- gsub("[[:punct:]]","",x)
    L[[i]] <- append(L[[i]],x)
    L[[i]] <- unique(L[[i]])
  }
  return(L)
}
L <- Ldat("kl2004.txt",1:10)
L <- Ldat("mac2011.txt",1:10)
L <- Ldat("patt2003.txt",c(3,2,1,6,8,4,10,5,7,9))

# Intersect LG markers with list, and get the chromosome with max intersections.

getchromosome <- function(x){
  m <- NULL
  for(i in 1:10){
    print(intersect(L[[i]],x))
    m <- append(m,length(intersect(L[[i]],x)))
  }
  return(which.max(m))
}

