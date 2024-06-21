#Figure 2C Boxplots showing significant differential Kyoto Encyclopedia of Genes and Genomes Orthologies (KOs). After the boxplots was produced by the p-value script, the significance calculated by p-value script  were marked in adobe Illustrator.
#以K03073为例
library(ggpubr)
library(ggplot2)
setwd("path/to/your/directory")  #设置工作路径
df <-read.table("kegg_4.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="Group", y="K03073", color = "Group", palette = "jco", add = "jitter")     
p
compare_means(K03073~Group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)


 
#Figure 2D Boxplots showing significant differential carbohydrate-active enzymes. After the boxplots was produced by the p-value script, the significance calculated by p-value script  were marked in adobe Illustrator.
#以GH13_18为例,
library(ggpubr)
library(ggplot2)
setwd("path/to/your/directory")  #设置工作路径
df <-read.table("cazy_4.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="Group", y="GH13_18", color = "Group", palette = "jco", add = "jitter",fill = "Group")     
p
compare_means(GH13_18~Group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)
