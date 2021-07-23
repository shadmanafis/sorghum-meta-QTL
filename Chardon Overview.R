library(dplyr)
library(tidyr)

setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results/Consensus Maps")
qtl <- read.table("ALLQTL.txt",sep = "\t",skip = 1) %>% group_by(V6)
colnames(qtl)[6] <- "chr"
Sb <- select(qtl,chr,V10,V11,V12) %>%
   transmute(chr=chr,mean=V10,sd=(V12-V11)/(2*1.96))
result <- data.frame()
output <- data.frame()
for(j in unique(Sb$chr))
  {
  S <- filter(Sb, chr==j)
    x <- seq(0,190)
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
col_fun = colorRamp2(c(min(meta$Wt), mean(meta$Wt), max(meta$Wt)), c("blue", "white", "red"))