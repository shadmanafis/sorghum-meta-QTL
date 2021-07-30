library(circlize)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(RColorBrewer)
library(yarrr)
#############################   DATA PREPARATION    #######################################################
source("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Physical Interval.R")
## Consensus Map Data
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
Consensus <- data.frame(chr=as.factor(Consensus$Chromosome),cM=Consensus$V3,Marker=tolower(Consensus$V2))
Consensus1 <- filter(Consensus,row_number() %% 3 == 1) ## Select every 3rd row starting from first row
Consensus1 <- filter(Consensus1,row_number() %% 3 == 1) ## Select every 3rd row starting from first row
Consensus1 <- filter(Consensus1,row_number() %% 3 == 1) ## Select every 3rd row starting from first row


## QTL Overview Data
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results/Consensus Maps")
qtl <- read.table("ALLQTL.txt",sep = "\t",skip = 1) %>% group_by(V6)
colnames(qtl)[6] <- "chr"
Sb <- dplyr::select(qtl,chr,V10,V11,V12) %>%
  transmute(chr=chr,mean=V10,sd=(V12-V11)/(2*1.96))
result <- data.frame()
output <- data.frame()
# Limit of x axis on all chromosomes
mx <- c()
for(j in unique(Consensus$chr)){mx <- c(mx,filter(Consensus, chr==j)%>%dplyr::select(cM)%>%max())}
# Values of y every cM on x axis
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


## Meta QTL Locations
setwd("C:/Users/Shabana Usmani/sorghum-meta/sorghum-meta-QTL/Results")
meta <- read.table("Meta-weight.txt", header = TRUE)
colnames(meta)[3] <- "chr"


## Physical Map Data Prep
pp95<-data.frame(reduce(p95))[,1:3]
pp95[1,3] <- end(ranges(GR[24621]))
pp95$seqnames <- c("10","02","02","03","03","03","04")
pp95 <- arrange(pp95, seqnames)
pp95$regions <- 1:7
reg_lim<-data.frame(regions=rep(1:7,2),x=c(pp95$start,pp95$end))
p95_red <- makeGRangesFromDataFrame(pp95, keep.extra.columns = TRUE)
seqlevels(GR) <- c("01","02","03","04","07","10")
GR <- makeGRangesFromDataFrame(arrange(data.frame(GR),data.frame(GR)$seqnames))
meta$MQTL <- row.names(meta)
meta$MQTL <- as.integer(meta$MQTL)
P95[13,2] <- end(ranges(GR[24621]))

# Number of Genes in a 250KB Bins in each region
all_counts <- NULL
for(i in 1:7){
  den <- GRanges(seqnames = pp95[i,1], ranges = IRanges(start=seq(pp95[i,2],pp95[i,3],by=250000), width = 250000), region=pp95[i,4])
  den$count <- countOverlaps(den, GR)
  all_counts <- rbind(all_counts, data.frame(den))
}

# Correspondance data for nested zooming
merge <- left_join(data.frame(P95),meta, by=c("MQTL"="MQTL"))
merge$region <- c(1,2,2,2,2,3,4,4,4,4,5,6,7)
merge <- merge[,c(7,5,6,9,1,2)]
merge2 <- NULL
for(j in unique(merge$region)){
  mr <- NULL
  mr <- split(merge, merge$region)[[j]] %>% summarise(across(c(2,3,5,6),c(min,max))) %>%
        dplyr::select(start.y_1,end.y_2,start.x_1,end.x_2)
  mr <- cbind(mr,split(merge, merge$region)[[j]][1,c(1,4)])
  merge2 <- rbind(merge2,mr)
}
merge2 <- merge2[,c(5,1,2,6,3,4)]

## SETTING COLOURS
tr_col <- brewer.pal(11, "BrBG")[6]
ov_col <- brewer.pal(11, "BrBG")[10:11]
col_fun = colorRamp2(c(min(meta$Wt), mean(meta$Wt), max(meta$Wt)), c("green", "black", "red"))
col_fun2 = colorRamp2(c(min(all_counts$count), max(all_counts$count)), c("#D9F0A3", "#004529"))
ht <- brewer.pal(9, "Greens")[c(1,4)]
gn <- brewer.pal(9, "Blues")
ni = brewer.pal(9,"Greens")


#############################   THE CIRCULAR PLOT    #######################################################

f1 <- function(){
  circos.par(cell.padding= c(0,0,0,0))
  circos.initialize(Consensus$chr, x=Consensus$cM)
  circos.track(Consensus$chr,x=Consensus$cM,ylim = c(0, 1), panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    xcenter=CELL_META$xcenter
    ycenter=CELL_META$ycenter
    circos.rect(xlim[1], 0, xlim[2], 0.3, col = rand_color(10, hue = "green", luminosity = "dark"), border=NA)
    circos.segments(x, rep(0,length(x)),x, rep(0.3,length(x)), col="black")
    x = Consensus1$cM[Consensus1$chr==CELL_META$sector.index]
    circos.segments(x, rep(0,length(x)),x, rep(0.4,length(x)), col="black")
    circos.text(x, rep(0.75,length(x)),Consensus1$Marker[Consensus1$chr==CELL_META$sector.index],niceFacing = TRUE,facing = "clockwise",cex=0.8, font = 3)
    circos.text(xcenter,ycenter/3,labels = CELL_META$sector.index, col = "floralwhite", font = 2)
  },bg.border = NA)
  circos.par(track.height=0.12)
  circos.track(result$V1,x=result$x,ylim = c(0, 0.15),panel.fun = function(x, y){
    circos.axis(h="top",direction = "inside" ,labels.facing = "clockwise",labels.cex = 0.5, labels.col = "grey17", col="grey17")
    circos.lines(result[result$V1==CELL_META$sector.index,2],result[result$V1==CELL_META$sector.index,3], area = TRUE, type = 'l', col= ov_col[1],border = ov_col[2])
    for(h in seq(0, 0.15, by = 0.05)) {
      circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
    }
  }, bg.border = NA, bg.col=tr_col)
  circos.par(track.height=0.05)
  circos.track(result$V1,x=result$x,ylim = c(0,1),panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    circos.rect(meta$start[meta$chr==chr], rep(0,length(meta$start[meta$chr==chr])), meta$end[meta$chr==chr], rep(1,length(meta$start[meta$chr==chr])), col = col_fun(meta$Wt[meta$chr==chr]), border=NA)
  }, bg.border = ov_col[2], bg.col=ht[1])
}

f2 <- function(){
  circos.par(start.degree=+60)
  circos.par(track.height=0.10)
  circos.par(track.margin=c(0,0))
  circos.initialize(reg_lim$regions, x=reg_lim$x)
  circos.track(c(all_counts$region,all_counts$region),x=c(all_counts$start,all_counts$end),ylim = c(0, 1), bg.col=ni[reg_lim$regions], track.index=1 ,panel.fun = function(x, y) {
    reg = CELL_META$sector.index
    xlim = CELL_META$xlim
    circos.text(CELL_META$xcenter, y=CELL_META$ycenter, labels=paste0("SBI-",pp95$seqnames[which(pp95$regions==reg)]),niceFacing = TRUE,facing = "inside",cex=0.6, col = "grey17")
    circos.axis(h="bottom",direction = "outside" ,labels.facing = "clockwise",major.at = reg_lim$x[reg_lim$regions==reg] ,labels = paste0(round(reg_lim$x[reg_lim$regions==reg]/10^6, digits = 2)),labels.cex = 0.7, labels.col = "grey17", col="grey17")
  },bg.border = NA)
  circos.par(track.height=0.15)
  circos.track(c(all_counts$region,all_counts$region),x=c(all_counts$start,all_counts$end),ylim = c(0, 1),track.index=2 ,panel.fun = function(x, y) {
    reg = CELL_META$sector.index
    circos.rect(all_counts[all_counts$region==reg,2], 0, all_counts[all_counts$region==reg,3], 1, col = col_fun2(all_counts[all_counts$region==reg,7]), border=NA)
  },bg.border = NA)
}

circos.nested(f1,f2,merge2, connection_height = mm_h(7), adjust_start_degree = FALSE, connection_col = ni[merge2[[4]]], connection_border = NA)
