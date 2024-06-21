library(vegan)
setwd("path/to/your/directory")  #设置工作路径
data_lxy_26 <- read.table("alpha data_-26d.txt", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy_26, index = "shannon")
Simpson <- diversity(data_lxy_26, index = "simpson")
div<-data.frame(Shannon.Wiener,Simpson)
write.csv(div,"Data_26.csv")
library(ggpubr)
library(ggplot2)
data_simpson_26 <-read.table("simpson_-26.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(data_simpson_26, x="group", y="Simpson", fill = "group", palette = "jco", add = "jitter",facet.by = "time")     
p
compare_means(Simpson~group, data=data_simpson_26, method = "wilcox.test", paired = F,  group.by = NULL)


data_lxy_4 <- read.table("alpha data_4d.txt", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy_4, index = "shannon")
Simpson <- diversity(data_lxy_4, index = "simpson")
div<-data.frame(Shannon.Wiener,Simpson)
write.csv(div,"Data_4.csv")
data_simpson_4 <-read.table("simpson_4.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(data_simpson_4, x="group", y="Simpson", fill = "group",palette = "jco", add = "jitter", facet.by = "time")   
p
compare_means(Simpson~group, data=data_simpson_4, method = "wilcox.test", paired = F,  group.by = NULL)
