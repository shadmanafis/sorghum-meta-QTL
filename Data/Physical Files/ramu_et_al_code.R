my_data  <- read.delim("Ramu et al.- Positions.txt", header = TRUE, sep = "\t", colClasses = "character")
my_data <- my_data[-grep("SB",my_data[,2], ignore.case=FALSE),]
my_data[,2] <- gsub("X*", "",my_data[,2], ignore.case = FALSE)