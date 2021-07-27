library(dplyr)
library(tidyr)

setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results/Consensus Maps")
qtl <- read.table("ALLQTL.txt",sep = "\t",skip = 1) %>% group_by(V6)
colnames(qtl)[6] <- "chr"
Sb <- dplyr::select(qtl,chr,V10,V11,V12) %>%
   transmute(chr=chr,mean=V10,sd=(V12-V11)/(2*1.96))
result <- data.frame()
output <- data.frame()

mx <- c()
for(j in unique(Consensus$chr))
{
  mx <- c(mx,filter(Consensus, chr==j)%>%dplyr::select(cM)%>%max())
}

for(j in unique(Sb$chr))
  {
  S <- filter(Sb, chr==j)
    x <- seq(0,mx[which(unique(Sb$chr)==j)])
    y <- 0
    for(i in 1:nrow(S)){
      y <- y + dnorm(x,mean = S$mean[i],sd=S$sd[i])
    }
    output <- cbind(rep(j,length(x)),x,y)
    result <- rbind(result,output)
  }
result[,1] <- as.factor(result[,1])

result$x <- as.numeric(result$x)
result$y <- as.numeric(result$y)
result$y <- result$y/10

# Import Meta QTL Locations
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results")
meta <- read.table("Meta-weight.txt", header = TRUE)
colnames(meta)[3] <- "chr"

# Data Cleaning for Physical Map
pp95<-data.frame(reduce(p95))[,1:3]
pp95$seqnames <- c("10","02","02","03","03","03","04")
pp95 <- arrange(pp95, seqnames)
pp95$regions <- 1:7
reg_lim<-data.frame(regions=rep(1:7,2),x=c(pp95$start,pp95$end))
p95_red <- makeGRangesFromDataFrame(pp95, keep.extra.columns = TRUE)
seqlevels(GR) <- c("01","02","03","04","07","10")
GR <- makeGRangesFromDataFrame(arrange(data.frame(GR),data.frame(GR)$seqnames))
all_genes <- NULL
for(i in 1:length(p95_red)){
  temp <- data.frame(subsetByOverlaps(GR,p95_red[i]))
  temp$region <- rep(i, nrow(temp))
  all_genes <- rbind(all_genes,temp)
}

reg_lim$x[14]<-tail(all_genes[4580,3])

meta$MQTL <- row.names(meta)
meta$MQTL <- as.integer(meta$MQTL)

merge <- left_join(data.frame(P95),meta, by=c("MQTL"="MQTL"))
merge$region <- c(1,2,2,2,2,3,4,4,4,4,5,6,7)
merge <- merge[,c(7,5,6,9,1,2)]

all_counts <- NULL
for(i in 1:7){
  den <- GRanges(seqnames = pp95[i,1], ranges = IRanges(start=seq(pp95[i,2],pp95[i,3],by=250000), width = 250000), region=pp95[i,4])
  den$count <- countOverlaps(den, GR)
  all_counts <- rbind(all_counts, data.frame(den))
}

all_counts <- all_counts[1:164,]
