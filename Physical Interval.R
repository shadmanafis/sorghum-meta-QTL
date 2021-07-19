library(readr)
library(tidyr)
library(dplyr)
library(GenomicRanges)

# Load Genetic Data
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Data/Physical Files")
phy <- read_tsv("Ramu_Jun.txt", col_names = TRUE) #%>%
  # filter(Chromosome==1)
phy <- dplyr::slice(phy,-grep("SB",phy$`Locus name`)) %>%
  dplyr::select(`Locus name`,`Physical Map Position`)
phy$`Locus name` <- gsub("X","",phy$`Locus name`) %>%
  tolower()
phy <- phy[!duplicated(phy$`Locus name`),]
# Set WOrking Directory, Granges list from all the consensus maps
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results/Consensus Maps")
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
Consensus <- data.frame(chr=Consensus$Chromosome,start=Consensus$V3,end = Consensus$V3+1,feature=tolower(Consensus$V2))
Consensus <- left_join(Consensus,phy, by=c("feature"="Locus name"))
Consensus <- Consensus[is.na(Consensus$`Physical Map Position`)==FALSE,]
Consensus <- Consensus[!duplicated(Consensus$feature),]

Temp <- split.data.frame(Consensus, Consensus$chr)
vec <- function(x){
  v <- 0
  for(i in 2:dim(x)[1]){
    v <- c(v,x[i,5]-x[i-1,5])
  }
  return(v)
}
Temp <- sapply(Temp, vec, simplify = `array`)
Temp <- unlist(Temp, use.names = FALSE)
Consensus$Diff <- Temp
Consensus <- Consensus[Consensus$Diff>=0,]

GR <- makeGRangesFromDataFrame(Consensus, keep.extra.columns=TRUE)

setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results")
meta <- read_tsv("Meta.txt")
gr <- makeGRangesFromDataFrame(meta)

# Get markers inside or flanking the Meta-QTL
df <- data.frame()
for(i in 1:length(ranges(gr))){
   hits <- subjectHits(findOverlaps(gr[i],GR))
   if(length(hits)!=0){
     df <- rbind(df,paste(GR[(min(hits)-1):(max(hits)+1)]$feature,collapse = ","))
   } else {
      hits <- nearest(gr[i],GR)
     df <- rbind(df,paste(GR[(min(hits)-1):(max(hits)+1)]$feature,collapse = ","))
   }
}

colnames(df) <- "markers"
values(gr) <- df
df <- data.frame(MQTL = row.names(df),markers = df$markers)%>%
  separate_rows(markers,sep = ",", convert = FALSE) %>%
  left_join(phy,by=c("markers"="Locus name"))

df <- left_join(df,data.frame(MQTL=row.names(meta),meta[,3]),by=c("MQTL"="MQTL"))
df <- right_join(df,Consensus, by=c("markers"="feature","chr"="chr"))[,1:5]
df <- df[!is.na(df$MQTL),]
df <- df[,c(1,4,2,5,3)]
df <- df[!duplicated(df),]  
colnames(df) <- c("MQTL","Chr","Marker","cM","Mb")

getinterval <- function(i){
  fit <- filter(df, df$MQTL==i)%>%
        lm(formula=Mb~cM)%>%predict(newdata = data.frame(cM=unlist(meta[i,1:2], use.names = FALSE)), interval = "confidence", level=0.95)
  fit <- c(fit[1,2],fit[dim(fit)[1],3])
  return(fit)
}
P95 <- sapply(1:32,getinterval)
P90 <- sapply(1:32,getinterval)
P75 <- sapply(1:32,getinterval)
