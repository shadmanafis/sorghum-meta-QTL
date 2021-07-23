library(circlize)
library(dplyr)
library(tidyr)
library(GenomicRanges)
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

circos.initialize(Consensus$chr, x=Consensus$cM)
circos.track(Consensus$chr,x=Consensus$cM,ylim = c(0, 1), panel.fun = function(x, y) {
          xlim = CELL_META$xlim
          circos.rect(xlim[1], 0, xlim[2], 0.3, col = rand_color(10, hue = "green", luminosity = "bright"), border=NA)
          circos.segments(x, rep(0,length(x)),x, rep(0.3,length(x)), col="grey")
          x = Consensus1$cM[Consensus1$chr==CELL_META$sector.index]
          circos.segments(x, rep(0,length(x)),x, rep(0.5,length(x)), col="black")
          circos.text(x, rep(1,length(x)),Consensus1$Marker[Consensus1$chr==CELL_META$sector.index],niceFacing = TRUE,facing = "clockwise",cex=0.8)
   },bg.border = NA)
circos.par(track.height=0.1)
circos.track(result$V1,x=result$x,ylim = c(0, 0.15),panel.fun = function(x, y) {
  circos.lines(result[result$V1==CELL_META$sector.index,2],result[result$V1==CELL_META$sector.index,3])
}, bg.border = NA)
circos.par(track.height=0.05)
circos.track(result$V1,x=result$x,ylim = c(0,1),panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  circos.rect(meta$start[meta$Chr==CELL_META$sector.index], rep(0,length(meta$start[meta$Chr==CELL_META$sector.index])), meta$end[meta$Chr==CELL_META$sector.index], rep(1,length(meta$start[meta$Chr==CELL_META$sector.index])), col = col_fun(meta$Wt[meta$Chr==chr]), border=NA)
})

