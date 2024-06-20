library(vegan)
data_lxy <- read.table("clipboard", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy, index = "shannon")
Simpson <- diversity(data_lxy, index = "simpson")
div<-data.frame(Shannon.Wiener,Simpson)
write.csv(div,"Data.csv")


