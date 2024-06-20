library(vegan)
setwd("path/to/your/directory")  #设置工作路径
data_lxy_26 <- read.table("alpha data_-26d.txt", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy_26, index = "shannon")
Simpson <- diversity(data_lxy_26, index = "simpson")
div<-data.frame(Shannon.Wiener,Simpson)
write.csv(div,"Data_26.csv")


data_lxy_4 <- read.table("alpha data_4d.txt", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy_4, index = "shannon")
Simpson <- diversity(data_lxy_4, index = "simpson")
div<-data.frame(Shannon.Wiener,Simpson)
write.csv(div,"Data_4.csv")

