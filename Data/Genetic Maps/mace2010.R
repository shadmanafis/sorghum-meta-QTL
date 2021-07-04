m <- read.delim("mac2010a consensus.txt", header = TRUE, sep = "\t", colClasses = "character")
m1<- m[grep("txp*",m[,6], ignore.case=FALSE),]