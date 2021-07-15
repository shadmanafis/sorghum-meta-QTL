names(gen_con) <- c("Haus","Sab","Klein","Mac","phy")

#INTERSECTIONS
nm <- combn(names(gen_con),2,FUN = paste0, collapse = "&",simplify = FALSE)
# Make the combinations of list elements
cmb <- combn(gen_con,2,simplify = FALSE)
# Intersect the list elements
int <- unlist(lapply(cmb, function(x) intersect( x[[1]] , x[[2]])))
out <- lapply(cmb, function(x) length( intersect( x[[1]] , x[[2]])))
# Output with names
setNames(out,nm)


for(i in 1:9){
  genetic[[i]]$V2 <- gsub("[[:punct:]]","",genetic[[i]]$V2)
  # genetic[[i]]$V2 <- gsub("psb0","psb",genetic[[i]]$V2)
  # genetic[[i]]$V2 <- gsub("umc0","umc",genetic[[i]]$V2)
  # genetic[[i]]$V2 <- gsub("csu0","csu",genetic[[i]]$V2)
  genetic[[i]]$V2 <- tolower(genetic[[i]]$V2)
  arrange(genetic[[i]],"V3")
}

write.table(genetic[[i]],paste0(i,"genetic.txt"),sep = "\t",quote = FALSE, row.names = TRUE)